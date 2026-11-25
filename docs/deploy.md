# Deploying the Microbiome Platform (for IT)

A Flask website plus a background worker, installed from this repository. The website only shows
pages and queues analyses; the worker runs them, one at a time. Your reverse proxy (nginx or similar)
serves HTTPS in front of the website. Ready-made files are in `deploy/`.

## Requirements

- **Linux x86-64**, 16 GB RAM or more (large analyses use 7–8 GB), a few GB of disk for the
  environment, and space for uploads and results.
- **uv 0.13.7** ([astral.sh/uv](https://docs.astral.sh/uv/)), installed where root finds it:
  `curl -LsSf https://astral.sh/uv/0.13.7/install.sh | sudo env UV_INSTALL_DIR=/usr/local/bin sh`.
  uv installs **Python 3.13** itself, so the system Python does not matter. Internet access to PyPI
  and `download.pytorch.org` is needed during installation only (the CPU build of PyTorch; no GPU).
- **System libraries for Qt** (one of the analysis tools draws its tree through PyQt5, with no
  screen). Without them that analysis fails with `libGL.so.1: cannot open shared object file`:
  - RHEL / AlmaLinux / Rocky: `dnf install libglvnd-glx mesa-libGL fontconfig libXext libX11`
  - Debian / Ubuntu: `apt install libgl1 libfontconfig1 libxext6 libx11-6`
- `git`, and a TLS certificate for the site's host name.

## Install

The layout below matches the service files in `deploy/systemd/`: code and environment in
`/opt/microplatform`, data in `/var/lib/microplatform`, settings in `/etc/microplatform/`, run by a
dedicated account `microplatform`.

```bash
sudo useradd --system --home-dir /var/lib/microplatform --shell /sbin/nologin microplatform
sudo install -d -m 750 -o microplatform -g microplatform /var/lib/microplatform
sudo git clone <repository URL> /opt/microplatform        # branch main
cd /opt/microplatform
# Keep uv's Python inside /opt: the services cannot read anyone's home folder (ProtectHome=true).
sudo env UV_PYTHON_INSTALL_DIR=/opt/microplatform/python uv sync --locked --extra server
sudo chmod -R a+rX /opt/microplatform                     # the service account must read the code
ls /opt/microplatform/.venv/bin/gunicorn                  # must exist
```

`--locked` installs exactly the tested versions in `uv.lock`; if it refuses, stop and tell us.

> **Always `uv sync --locked --extra server`. Never a plain `uv sync`, and never `uv run …`.**
> gunicorn comes only with the `server` extra. A plain `uv sync` (and `uv run`, which syncs first)
> makes the environment match the lock *without* extras and **removes gunicorn**. The site that is
> already running keeps working, so nothing looks wrong, until the next restart or reboot, when it
> fails to start. To run a command in the environment, use `/opt/microplatform/.venv/bin/python -m …`.

## Configure

```bash
sudo install -d -m 750 -o root -g microplatform /etc/microplatform
sudo install -m 640 -o root -g microplatform deploy/microplatform.env.example /etc/microplatform/microplatform.env
python3 -c "import secrets; print(secrets.token_hex(32))"     # paste as MICROPLATFORM_SECRET_KEY
sudo nano /etc/microplatform/microplatform.env
```

Set at least:

| Setting | Value |
|---|---|
| `MICROPLATFORM_SECRET_KEY` | the fresh random value printed above (never reuse one, never commit it) |
| `MICROPLATFORM_PUBLIC_URL` | the site's address, e.g. `https://microbiome.example.ac.il` (used in results links) |
| `MICROPLATFORM_BIND` | where gunicorn listens: `127.0.0.1:<port>` on this machine (default `127.0.0.1:8000`) |
| `MICROPLATFORM_TRUST_PROXY`, `MICROPLATFORM_COOKIE_SECURE` | `1` and `1` behind an HTTPS reverse proxy |
| `MICROPLATFORM_MAX_UPLOAD_MB` | `100` (keep the proxy's body-size limit the same) |
| `MICROPLATFORM_SMTP_*` | **leave empty**: email stays off until the lab provides mail-server details |

The template explains every line. Uploads can be clinical data: keep `/var/lib/microplatform` private
to the service account and do not back it up to places the lab has not approved.

## Run

Two services, from `deploy/systemd/`:

```bash
sudo cp deploy/systemd/microplatform-web.service deploy/systemd/microplatform-worker.service /etc/systemd/system/
sudo systemctl daemon-reload
sudo systemctl enable --now microplatform-web microplatform-worker
systemctl is-active microplatform-web microplatform-worker      # active, active
curl -s -o /dev/null -w "%{http_code}\n" http://127.0.0.1:8000/  # 200 (use your port)
```

The website service runs:

```
/opt/microplatform/.venv/bin/gunicorn -c /opt/microplatform/deploy/gunicorn.conf.py microplatform.web.wsgi:app
```

and the worker `/opt/microplatform/.venv/bin/python -m microplatform.worker --instance /var/lib/microplatform`.
gunicorn binds to `MICROPLATFORM_BIND` and trusts forwarded headers only from `127.0.0.1`
(`deploy/gunicorn.conf.py`; change `forwarded_allow_ips` if the proxy is on another machine).
On SELinux systems, allow the proxy to connect to the port (`setsebool -P httpd_can_network_connect 1`).

**Reverse proxy:** `deploy/nginx/microplatform.conf` (and `microplatform.almalinux.conf`) is a
complete example: HTTP → HTTPS redirect, TLS, an allow-list for the lab network, a 100 MB body
limit, and `proxy_pass http://127.0.0.1:8000`. Replace the `<HOSTNAME>`, certificate paths, the
allowed address ranges and the port. The site has **no login**: anyone who can reach it can use it,
so restrict it to the lab network or VPN.

The worker's service waits up to 4 hours on a stop or restart, so a long analysis finishes instead
of being cut (`TimeoutStopSec=4h`); `systemctl kill microplatform-worker` stops it at once.

## Check

Open the site, click **Try the example**, run it, and download the results ZIP. Logs:
`journalctl -u microplatform-web -u microplatform-worker`.

## Update

```bash
cd /opt/microplatform && sudo git pull --ff-only
sudo env UV_PYTHON_INSTALL_DIR=/opt/microplatform/python uv sync --locked --extra server
sudo chmod -R a+rX /opt/microplatform
sudo systemctl restart microplatform-web microplatform-worker
```

Restarting the worker waits for an analysis in progress to finish.

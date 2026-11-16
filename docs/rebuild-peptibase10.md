# Rebuilding the platform on peptibase10 from scratch

For the day the server is wiped or reinstalled. The platform runs under the account `parashn`
as two systemd *user* services, with no administrator rights; everything lives in the home folder.
Rebuilding takes about 30 minutes, most of it the package download (PyTorch).

**What survives a wipe, and what does not**

| Thing | Where it comes back from |
|---|---|
| The code and the release tags | GitHub (`noaparash/microbiome-platform`, private), through the laptop |
| The Python environment | Rebuilt exactly from `uv.lock` (step 4) |
| The service definitions | `deploy/systemd-user/` in the repo (step 7) |
| The settings file | Recreated from `deploy/microplatform.user.env.example` with a **new** secret key (step 6) |
| Uploads, results, cache (`~/microplatform-instance`) | **Only from your own backup copy** (step 8). Not in Git, and never to be put there: it may hold clinical data |
| Logs | Not restored (start empty) |

The server has no access to the private GitHub repo, so the code always travels
**laptop → server**: the laptop pulls from GitHub and pushes to a bare repository on the server.

---

## 1. Reach the server from the laptop

Start the SSH tunnel (Bitvise) as usual, then in PowerShell:

```powershell
ssh -4 peptibase10 "echo ok"
```

- **"REMOTE HOST IDENTIFICATION HAS CHANGED"** (a reinstalled server has a new host key):
  `ssh-keygen -R "[localhost]:222"`, then connect again and accept the new key.
- **Asked for a password** (the wipe removed `~/.ssh/authorized_keys`): the password still works.
  Put the laptop's key back so later steps don't ask every time:

  ```powershell
  type $env:USERPROFILE\.ssh\id_ed25519_peptibase10.pub | ssh -4 peptibase10 "mkdir -p ~/.ssh && chmod 700 ~/.ssh && cat >> ~/.ssh/authorized_keys && chmod 600 ~/.ssh/authorized_keys && restorecon -R ~/.ssh"
  ```

  (`restorecon` matters: on this server the home folder may be unlabelled for SELinux, and then
  the SSH server ignores the key.)

Always use `ssh -4` / `scp -4`: the tunnel listens on IPv4 only.

## 2. Install uv (on the server)

The same way it was installed before (Astral's official installer, into `~/.local/bin`):

```bash
curl -LsSf https://astral.sh/uv/0.13.7/install.sh | sh
~/.local/bin/uv --version        # uv 0.13.7
```

If that download is blocked, the same uv comes from PyPI: `python3 -m pip install --user uv==0.13.7`.
uv later downloads its own Python 3.13, so the system Python does not matter.

## 3. The code

On the server, an empty repository for the laptop to push to:

```bash
mkdir -p ~/git && git init --bare ~/git/microbiome-platform.git
```

On the laptop (in the project folder; the remote `peptibase10` already points there):

```powershell
git pull origin main
$env:GIT_SSH_COMMAND = "ssh -4"; git push peptibase10 main --tags
```

On the server, the running copy:

```bash
git clone ~/git/microbiome-platform.git ~/microbiome-platform
```

## 4. The Python environment (on the server)

```bash
cd ~/microbiome-platform && ~/.local/bin/uv sync --locked --extra server
```

`--locked` installs exactly the tested versions in `uv.lock` (CPU-only PyTorch). If it refuses
because the lock is out of date, stop: something was changed without updating the lock.

## 5. Private folders

```bash
mkdir -m 700 -p ~/microplatform-instance ~/microplatform-logs ~/.config/microplatform
```

## 6. The settings file (never in Git)

```bash
install -m 600 ~/microbiome-platform/deploy/microplatform.user.env.example ~/.config/microplatform/microplatform.env
F=~/.config/microplatform/microplatform.env
sed -i "s/^MICROPLATFORM_SECRET_KEY=.*/MICROPLATFORM_SECRET_KEY=$(python3 -c 'import secrets; print(secrets.token_hex(32))')/" $F
sed -i "s#/home/<you>/#$HOME/#g" $F
grep -c "SECRET_KEY=." $F         # 1: the key is set
```

The example already holds this server's choices (address `127.0.0.1:8000`, 16 threads per
analysis, no HTTPS yet). A new secret key is fine: it only signs short-lived browser cookies.
Email stays off until the SMTP lines are filled in (see `deploy/microplatform.user.env.example`).

## 7. The services

```bash
mkdir -p ~/.config/systemd/user && cp ~/microbiome-platform/deploy/systemd-user/*.service ~/.config/systemd/user/
systemctl --user daemon-reload && systemctl --user enable --now microplatform-web microplatform-worker
loginctl enable-linger parashn     # keeps them running without an open session, and starts them at boot
restorecon -R /home/parashn        # SELinux labels for the home folder (worked without sudo before)
```

If `loginctl enable-linger` is refused after the update, ask the administrator to run
`sudo loginctl enable-linger parashn` (see `docs/admin-peptibase10.md`). Until then the services
stop when your last session on the server closes.

## 8. Restore the data (only if you kept a copy)

Stop the services, copy the saved `microplatform-instance` folder back, keep it private, start again:

```bash
systemctl --user stop microplatform-web microplatform-worker
# copy the backup into ~/microplatform-instance (same layout: runs/, cache/, ...)
chmod -R go-rwx ~/microplatform-instance
systemctl --user start microplatform-web microplatform-worker
```

Without a copy the site starts empty; old results links then show "not found".

## 9. Check that it works

On the server:

```bash
systemctl --user is-active microplatform-web microplatform-worker      # active, active
curl -s -o /dev/null -w "%{http_code}\n" http://127.0.0.1:8000/          # 200
```

From the laptop: `ssh -N -L 8000:127.0.0.1:8000 parashn@peptibase10.math.biu.ac.il`, open
http://127.0.0.1:8000, click **Try the example**, run it, and download the ZIP.

**LOCATE's exact numbers depend on the processor.** If the update moved the platform to different
hardware, check before offering LOCATE:

```bash
cd ~/microbiome-platform && ~/.local/bin/uv run pytest tests/test_reproduction.py -k locate -q
```

If that fails (or skips with "re-validate first"), switch LOCATE off until it is re-validated:
add `metabolites` to `MICROPLATFORM_DISABLED_GOALS` in the settings file and restart the services.

---

## Before a planned update: what to save

1. Make sure GitHub has every commit and tag on the server (`git ls-remote origin` on the laptop
   against `git -C ~/git/microbiome-platform.git for-each-ref` on the server).
2. Copy `~/microplatform-instance` to a place the lab has approved for research data (not GitHub,
   not a personal cloud drive), if any run in it is worth keeping. The cache folder inside it can
   be left out: it is rebuilt automatically.

#!/bin/bash
# Deploy the laptop's current main to both lab servers (peptibase10, then its mirror peptibase9).
#
# Run from the repository on the laptop, in Git Bash:   bash deploy/deploy-both.sh
#
# 1. Refuses unless you are on main with no uncommitted changes to tracked files.
# 2. Pushes main to the GitHub backup (origin) and to each server's bare repository.
# 3. On each server: skips it if an analysis is queued or running (a restart would wait for it);
#    otherwise git pull --ff-only, uv sync --locked --extra server (never a plain uv sync: it
#    removes gunicorn), checks that gunicorn is there, restarts both user services, and checks
#    that the site answers on the port in its own settings file (MICROPLATFORM_BIND).
# 4. Checks that both servers now run the same commit as the laptop.
#
# Uses the ssh aliases peptibase10 and peptibase9 from ~/.ssh/config (always over IPv4).
set -euo pipefail
SERVERS=(peptibase10 peptibase9)

cd "$(git rev-parse --show-toplevel)"
[ "$(git rev-parse --abbrev-ref HEAD)" = main ] || { echo "Not on main: stopping."; exit 1; }
git diff --quiet && git diff --cached --quiet || { echo "Uncommitted changes to tracked files: commit them first."; exit 1; }
want=$(git rev-parse HEAD)
echo "Deploying $(git log --oneline -1)"

git push -q origin main
for s in "${SERVERS[@]}"; do git push -q "$s" main; done
echo "Pushed main to origin ${SERVERS[*]}"

server_step() {  # runs on the server
  ssh -4 -o BatchMode=yes "$1" 'bash -s' <<'EOF'
set -eu
cd ~/microbiome-platform
busy=$(grep -l -E '"status": "(queued|running)"' ~/microplatform-instance/runs/*/meta.json 2>/dev/null | wc -l)
if [ "$busy" -gt 0 ]; then echo "  SKIPPED: $busy analysis queued or running; deploy again when it ends."; exit 3; fi
git pull -q --ff-only
~/.local/bin/uv sync -q --locked --extra server
[ -x .venv/bin/gunicorn ] || { echo "  gunicorn missing after uv sync: NOT restarting."; exit 4; }
systemctl --user restart microplatform-web microplatform-worker
bind=$(sed -n 's/^MICROPLATFORM_BIND=//p' ~/.config/microplatform/microplatform.env); bind=${bind:-127.0.0.1:8000}
for i in $(seq 1 30); do code=$(curl -s -o /dev/null -w '%{http_code}' "http://$bind/" || true); [ "$code" = 200 ] && break; sleep 1; done
echo "  services: $(systemctl --user is-active microplatform-web microplatform-worker | tr '\n' ' ')| site http://$bind/: $code"
echo "  running: $(git log --oneline -1)"
[ "$code" = 200 ]
EOF
}

status=0
for s in "${SERVERS[@]}"; do
  echo "== $s"
  server_step "$s" 2> >(grep -v 'bad CNAME' >&2) || { echo "  $s: deploy step FAILED (exit $?)"; status=1; }
done

echo "== check"
for s in "${SERVERS[@]}"; do
  got=$(ssh -4 -o BatchMode=yes "$s" 'git -C ~/microbiome-platform rev-parse HEAD' 2>/dev/null)
  if [ "$got" = "$want" ]; then echo "  $s: ${got:0:7} = laptop"; else echo "  $s: ${got:0:7} DIFFERS from laptop ${want:0:7}"; status=1; fi
done
exit $status

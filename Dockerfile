FROM python:3.13-slim-bookworm

COPY --from=ghcr.io/astral-sh/uv:0.13.7 /uv /usr/local/bin/uv

WORKDIR /app

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
       libgl1 libfontconfig1 libxext6 libx11-6 libglib2.0-0 libgomp1 \
    && rm -rf /var/lib/apt/lists/*

COPY . .

RUN uv sync --locked --extra server --no-dev --no-editable --no-cache

EXPOSE 5000

RUN uv add gunicorn

CMD ["/app/.venv/bin/gunicorn", "-c", "deploy/gunicorn.conf.py", "--bind", "0.0.0.0:5000", "microplatform.web.wsgi:app"]

FROM python:3.11-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1

WORKDIR /app

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    gcc \
    g++ \
    curl \
    ca-certificates \
    zlib1g \
    zlib1g-dev \
 && rm -rf /var/lib/apt/lists/*

COPY requirements*.txt ./ 
RUN if [ -f requirements.txt ]; then pip install -r requirements.txt; fi

COPY . .

RUN chmod +x ./run_server.sh

ENV PYTHONPATH=/app

EXPOSE 8000 8501

ENTRYPOINT ["./run_server.sh"]

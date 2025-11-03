FROM python:3.11-slim

ENV PYTHONUNBUFFERED=1 

WORKDIR /app

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    gcc \
    g++ \
    curl \
    ca-certificates \
    zlib1g \
    zlib1g-dev \
    libgomp1 \
    libfreetype6 \
    libpng16-16 \
 && rm -rf /var/lib/apt/lists/*

COPY requirements.txt ./

RUN python -m pip install --upgrade pip && \
    pip install --no-cache-dir -r requirements.txt

COPY . .

RUN chmod +x ./run_server.sh
ENV PYTHONPATH=/app

EXPOSE 8000 8501

ENTRYPOINT ["./run_server.sh"]

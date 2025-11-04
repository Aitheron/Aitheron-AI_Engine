FROM python:3.11-slim AS builder
ENV PIP_NO_CACHE_DIR=1 PIP_DISABLE_PIP_VERSION_CHECK=1 PIP_PREFER_BINARY=1
WORKDIR /wheels
RUN apt-get update && apt-get install -y --no-install-recommends gcc g++ zlib1g-dev && rm -rf /var/lib/apt/lists/*
COPY requirements.txt .
RUN python -m pip install --upgrade pip
RUN pip wheel --no-cache-dir --no-deps --index-url https://download.pytorch.org/whl/cpu torch==2.9.0 -w /wheels
RUN sed -i '/^torch==/d' requirements.txt
RUN pip wheel --no-cache-dir -r requirements.txt -w /wheels

FROM python:3.11-slim
ENV PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1 PIP_NO_CACHE_DIR=1 PIP_DISABLE_PIP_VERSION_CHECK=1
WORKDIR /app
RUN apt-get update && apt-get install -y --no-install-recommends zlib1g libgomp1 libfreetype6 libpng16-16 && rm -rf /var/lib/apt/lists/*
COPY --from=builder /wheels /tmp/wheels
RUN python -m pip install --upgrade pip && pip install --no-cache-dir /tmp/wheels/*.whl --no-compile && rm -rf /tmp/wheels
COPY . .
RUN chmod +x ./run_server.sh
ENV PYTHONPATH=/app
EXPOSE 8000 8501
ENTRYPOINT ["./run_server.sh"]

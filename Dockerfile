FROM python:3.12-slim

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
      build-essential \
      gcc \
      python3-dev \
      libffi-dev \
      libssl-dev \
      rustc \
      cargo \
      libxrender1 \
      libxext6 \
      libsm6 \
      libfontconfig1 \
      libgl1-mesa-glx && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /app

COPY requirements-prod.txt .

RUN pip install --upgrade pip setuptools wheel \
 && pip install --no-cache-dir -r requirements-prod.txt

COPY . .

EXPOSE 5000
CMD ["python", "run.py"]

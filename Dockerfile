FROM python:3.11-slim

LABEL maintainer="snailQH"
LABEL description="ab1tools: comprehensive toolkit for AB1 chromatogram files"
LABEL version="2.0.0"

# Install system dependencies for pysam
RUN apt-get update && apt-get install -y --no-install-recommends \
    gcc \
    libc6-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Install Python dependencies
COPY pyproject.toml .
COPY ab1tools/ ab1tools/
RUN pip install --no-cache-dir .

# Verify installation
RUN ab1tools --help

# Default working directory for data
WORKDIR /data

ENTRYPOINT ["ab1tools"]
CMD ["--help"]

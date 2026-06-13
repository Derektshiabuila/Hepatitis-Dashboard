# Use Debian base image with multi-architecture support
FROM debian:bookworm-slim

# Prevent interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies (including Wine, Xvfb for headless execution)
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    curl \
    git \
    ca-certificates \
    unzip \
    gnupg \
    procps \
    xvfb \
    xauth \
    docker.io \
    && dpkg --add-architecture i386 \
    && apt-get update \
    && apt-get install -y --no-install-recommends \
        wine \
        wine32:i386 \
        wine64 \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install Miniconda
ENV CONDA_DIR=/opt/conda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    /bin/bash /tmp/miniconda.sh -b -p $CONDA_DIR && \
    rm /tmp/miniconda.sh

# Put conda on PATH
ENV PATH=$CONDA_DIR/bin:$PATH

# Create the working directory
WORKDIR /app

# Copy the environment file and create the phylo Conda environment
COPY envs/phylo.yaml /app/envs/phylo.yaml
RUN conda config --remove channels defaults || true && \
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main || true && \
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r || true && \
    conda env create -f /app/envs/phylo.yaml && \
    conda clean -afy

# Set path to use the conda environment binaries by default
ENV PATH=/opt/conda/envs/phylo/bin:$PATH

# Install Python dependencies from requirements.txt
COPY requirements.txt /app/requirements.txt
RUN pip install --no-cache-dir -r requirements.txt

# Initialize Wine in headless mode
ENV WINEPREFIX=/root/.wine
ENV WINEDEBUG=-all
ENV WINEDLLOVERRIDES="mscoree,mshtml,winhttp=d"
RUN wineboot --init

# Attempt to download and silent-install RDP5 inside Wine
# Falls back gracefully to the pre-loaded scripts/RDP5CL.exe if UCT servers are unreachable
RUN mkdir -p /root/.wine/drive_c/RDP5 && \
    wget -q http://www.cs.uct.ac.za/~darren/RDP5/RDP5_Setup.zip -O /tmp/rdp5.zip || true && \
    if [ -f /tmp/rdp5.zip ]; then \
        unzip -q /tmp/rdp5.zip -d /tmp/rdp5_extracted && \
        find /tmp/rdp5_extracted -name "*Setup*.exe" -exec xvfb-run wine {} /VERYSILENT /SUPPRESSMSGBOXES /DIR="C:\RDP5" \; || true; \
    fi && \
    rm -rf /tmp/rdp5.zip /tmp/rdp5_extracted

# Copy the rest of the application files
COPY . /app/

# Set environment variables for RDP5 runner script
# Points to either the silent-installed path or falls back to /app/scripts/RDP5CL.exe
ENV RDP5_EXE=/root/.wine/drive_c/RDP5/RDP5CL.exe
ENV RDP5_WINE=wine

# Expose the dashboard port
EXPOSE 8051

# Start the dashboard using Gunicorn
CMD ["gunicorn", "Full_Hepatitis_page:server", "--bind", "0.0.0.0:8051", "--workers", "2", "--timeout", "120"]

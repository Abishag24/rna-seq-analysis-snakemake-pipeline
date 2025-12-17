FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies
RUN apt-get update && apt-get install -y \
    wget \
    curl \
    unzip \
    git \
    build-essential \
    python3 \
    python3-pip \
    r-base \
    default-jre \
    vim-common \
    && rm -rf /var/lib/apt/lists/*

# Install bioinformatics tools
RUN apt-get update && apt-get install -y \
    fastqc \
    samtools \
    subread \
    && rm -rf /var/lib/apt/lists/*

# Install STAR
RUN wget https://github.com/alexdobin/STAR/archive/2.7.11b.tar.gz && \
    tar -xzf 2.7.11b.tar.gz && \
    cd STAR-2.7.11b/source && make STAR && \
    cp STAR /usr/local/bin && \
    cd / && rm -rf STAR-2.7.11b*

# Install fastp
RUN wget http://opengene.org/fastp/fastp && \
    chmod +x fastp && \
    mv fastp /usr/local/bin/

# Install Snakemake
RUN pip3 install snakemake

# Set working directory
WORKDIR /pipeline

# Copy workflow files only
COPY Snakefile .
COPY workflow/ workflow/
COPY scripts/ scripts/
COPY envs/ envs/
COPY README.md .

CMD ["bash"]

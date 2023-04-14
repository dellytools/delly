# use the ubuntu base image
FROM ubuntu:22.04

MAINTAINER Tobias Rausch rausch@embl.de

# install required packages
RUN apt-get update && apt-get install -y \
    autoconf \
    build-essential \
    cmake \
    g++ \
    gfortran \
    git \
    libcurl4-gnutls-dev \
    hdf5-tools \
    libboost-date-time-dev \
    libboost-program-options-dev \
    libboost-system-dev \
    libboost-filesystem-dev \
    libboost-iostreams-dev \
    libbz2-dev \
    libdeflate-dev \
    libhdf5-dev \
    libncurses-dev \
    liblzma-dev \
    pkg-config \
    zlib1g-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# set environment
ENV BOOST_ROOT /usr

# install delly
RUN cd /opt \
    && git clone --recursive https://github.com/dellytools/delly.git \
    && cd /opt/delly/ \
    && make STATIC=1 PARALLEL=1 all \
    && make install


# Multi-stage build
FROM alpine:latest
RUN mkdir -p /opt/delly/bin
WORKDIR /opt/delly/bin
COPY --from=0 /opt/delly/bin/delly .

# Workdir
WORKDIR /root/

# Add Delly to PATH
ENV PATH="/opt/delly/bin:${PATH}"

# Set OpenMP threads
ENV OMP_NUM_THREADS 2

# by default /bin/sh is executed
CMD ["/bin/sh"]

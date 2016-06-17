# use the ubuntu base image
FROM ubuntu:14.04

MAINTAINER Tobias Rausch rausch@embl.de

# install required packages
RUN apt-get update && apt-get install -y \
    ant \
    build-essential \
    cmake \
    g++ \
    gfortran \
    git \
    hdf5-tools \
    libboost-date-time-dev \
    libboost-program-options-dev \
    libboost-system-dev \
    libboost-filesystem-dev \
    libboost-iostreams-dev \
    libbz2-dev \
    libhdf5-dev \
    libncurses-dev \
    python \
    python-dev \
    python-pip \
    zlib1g-dev \
    && apt-get clean

# set environment
ENV BOOST_ROOT /usr
ENV SEQTK_ROOT /opt/htslib

# install delly
RUN cd /opt \
    && git clone https://github.com/samtools/htslib.git \
    && cd /opt/htslib \
    && make \
    && make lib-static \
    && make install
RUN cd /opt \
    && git clone https://github.com/samtools/bcftools.git \
    && cd /opt/bcftools \
    && make all \
    && make install
RUN cd /opt \
    && git clone https://github.com/samtools/samtools.git \
    && cd /opt/samtools \
    && make all \
    && make install
RUN cd /opt \
    && git clone https://github.com/tobiasrausch/delly.git \
    && cd /opt/delly/ \
    && touch .bcftools .htslib .boost \
    && make all \
    && install -p /opt/delly/src/delly /usr/local/bin/ \
    && install -p /opt/delly/src/stats /usr/local/bin/ \
    && install -p /opt/delly/src/cov /usr/local/bin/

# install python variant filtering dependencies
#RUN pip install -r delly/complexVariants/requirements.txt


# Delly Dockerfile
#
# Version 0.0.5
# build with: sudo docker build --rm --no-cache -t delly - < Dockerfile
# run with: sudo docker run -i -t delly /bin/bash

# use the ubuntu base image
FROM ubuntu

MAINTAINER Tobias Rausch rausch@embl.de

# update package repository
RUN apt-get update

# install g++, git, zlib, cmake, ant, python
RUN apt-get install -y build-essential g++ git cmake zlib1g-dev ant libbz2-dev libboost-date-time-dev libboost-program-options-dev libboost-system-dev libboost-filesystem-dev libboost-iostreams-dev python python-dev python-pip gfortran libncurses-dev hdf5-tools libhdf5-dev

# set environment
ENV BOOST_ROOT /usr
ENV SEQTK_ROOT /htslib

# install delly
RUN git clone https://github.com/samtools/htslib.git
RUN cd htslib && make && make lib-static && cd ..
RUN git clone https://github.com/samtools/bcftools.git
RUN cd bcftools && make all && cd ..
RUN git clone https://github.com/samtools/samtools.git
RUN cd samtools && make all && cd ..
RUN git clone https://github.com/tobiasrausch/delly.git
RUN cd delly/ && touch .bcftools .htslib .boost && make all && cd ..

# install python variant filtering dependencies
RUN pip install -r delly/complexVariants/requirements.txt




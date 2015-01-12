# Delly Dockerfile
#
# Version 0.0.3
# build with: sudo docker build --rm --no-cache -t delly - < Dockerfile
# run with: sudo docker run -i -t delly /bin/bash

# use the ubuntu base image
# FROM ubuntu:13.10
FROM ubuntu

MAINTAINER Tobias Rausch rausch@embl.de

# update package repository
RUN apt-get update

# install g++, git, zlib, cmake, boost, ...
RUN apt-get install -y build-essential g++ git cmake zlib1g-dev libboost-date-time-dev libboost-program-options-dev libboost-system-dev libboost-filesystem-dev libboost-iostreams-dev

# set environment
ENV BOOST_ROOT /usr

# install delly
RUN git clone --recursive https://github.com/tobiasrausch/delly.git
RUN cd delly/ ; make all

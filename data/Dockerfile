FROM ubuntu:18.04

# Install system dependencies
RUN apt-get update
RUN apt-get install -qq -y git curl wget nano pkg-config libfreetype6-dev libpng-dev zlib1g-dev gfortran apt-utils

# Install python 3
RUN apt-get install -qq -y software-properties-common && \
    add-apt-repository ppa:deadsnakes/ppa && \
    apt-get update && \
    apt-get install -qq -y python3.6 python3-pip python3-dev

# Download python dependencies
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get install -qq -y libproj-dev proj-bin proj-data libtcmalloc-minimal4 libhdf5-serial-dev python3-grib python3-tk
RUN pip3 install pyproj urllib3 'numpy>=1.9.0' pandas matplotlib tables

RUN mkdir data

# TODO: move out of the build phase
RUN wget -O ssi71_location.h5 https://github.com/KMarshland/valbal-proc-data/blob/master/ssi71_location.h5?raw=true

ADD atmotools.py atmotools.py


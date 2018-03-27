FROM debian:testing
MAINTAINER Nicolas Richart <nicolas.richart@epfl.ch>

# Make sure the package repository is up to date.
RUN apt-get -qq update && apt-get -qq -y install \
    g++ gfortran  cmake \
    libmumps-seq-dev libscotch-dev \
    libboost-dev libopenblas-dev \
    python3 python3-numpy python3-scipy \
    swig3.0 gmsh \
    && rm -rf /var/lib/apt/lists/*

# apt-get on one line due to https://docs.docker.com/develop/develop-images/dockerfile_best-practices/#run

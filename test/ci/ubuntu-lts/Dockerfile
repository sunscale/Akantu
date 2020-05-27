FROM ubuntu:xenial
MAINTAINER Nicolas Richart <nicolas.richart@epfl.ch>

ENV TZ=Europe/Zurich
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# Make sure the package repository is up to date.
RUN apt-get -qq update && apt-get -qq -y install \
    g++ gfortran  cmake ccache \
    openmpi-bin libmumps-dev libscotch-dev \
    libboost-dev libopenblas-dev libcgal-dev \
    python3 python3-dev \
    python3-numpy python3-scipy python3-mpi4py\
    python3-yaml \
    python3-pytest \
    gmsh curl git xsltproc \
    && rm -rf /var/lib/apt/lists/*

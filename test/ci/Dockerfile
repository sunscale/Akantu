FROM debian:testing
MAINTAINER Nicolas Richart <nicolas.richart@epfl.ch>

# Make sure the package repository is up to date.
RUN apt-get -qq update && apt-get -qq -y install \
    g++ gfortran  cmake \
    libmumps-seq-dev libscotch-dev \
    libboost-dev libopenblas-dev \
    python3 python3-dev \
    python3-numpy python3-scipy python3-mpi4py\
    python3-phabricator python3-click python3-yaml \
    python3-pytest \
    swig3.0 gmsh curl flake8 \
    git clang-format xsltproc jq \
    php-cli php-curl php-xml \
    && rm -rf /var/lib/apt/lists/*

# apt-get on one line due to https://docs.docker.com/develop/develop-images/dockerfile_best-practices/#run
RUN git clone https://github.com/phacility/libphutil.git /libphutil
RUN git clone https://github.com/phacility/arcanist.git /arcanist

ENV PATH="$PATH:/arcanist/bin/"
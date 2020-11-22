FROM debian:testing
MAINTAINER Nicolas Richart <nicolas.richart@epfl.ch>

# library dependencies
RUN apt-get -qq update && apt-get -qq -y install \
    g++ gfortran clang cmake \
    openmpi-bin libmumps-dev libscotch-dev \
    libboost-dev libopenblas-dev \
    python3 python3-dev python3-numpy python3-scipy python3-mpi4py \
    && rm -rf /var/lib/apt/lists/*

# for documentation
RUN apt-get -qq update && apt-get -qq -y install \
    python3-sphinx python3-sphinx-rtd-theme \
    python3-sphinxcontrib.bibtex \
    python3-breathe doxygen graphviz \
    && rm -rf /var/lib/apt/lists/*

# for ci
RUN apt-get -qq update && apt-get -qq -y install \
    gmsh python3-pytest \
    ccache clang-format python3-flake8 clang-tidy \
    curl git xsltproc jq \
    gcovr binutils \
    && rm -rf /var/lib/apt/lists/*

COPY .openmpi /root/.openmpi

# for debug
RUN apt-get -qq update && apt-get -qq -y install \
    gdb valgrind \
    && rm -rf /var/lib/apt/lists/*

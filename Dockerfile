FROM debian:testing
MAINTAINER Nicolas Richart <nicolas.richart@epfl.ch>

# Make sure the package repository is up to date.
RUN apt-get -qq update
RUN apt-get -qq -y upgrade

RUN apt-get -qq -y install g++ gfortran libmumps-seq-dev libscotch-dev libboost-dev libopenblas-dev gmsh cmake

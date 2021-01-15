FROM quay.io/pypa/manylinux2010_x86_64
MAINTAINER Guillaume Anciaux <guillaume.anciaux@epfl.ch>

RUN git clone https://github.com/spack/spack.git
RUN yum install -y xz
RUN mkdir -p /root/.spack/linux/
COPY compilers.yaml /root/.spack/linux/
RUN spack/bin/spack install perl
# RUN /spack/bin/spack install akantu@master external_solvers=mumps ~mpi ^mumps~mpi


FROM quay.io/pypa/manylinux2010_x86_64
MAINTAINER Guillaume Anciaux <guillaume.anciaux@epfl.ch>

RUN git clone https://github.com/spack/spack.git
RUN yum install -y xz
RUN mkdir -p /root/.spack/linux/
COPY compilers.yaml /root/.spack/linux/
RUN spack/bin/spack install perl
RUN /spack/bin/spack install boost
RUN /spack/bin/spack install cmake
RUN /spack/bin/spack install openblas
RUN /spack/bin/spack install mumps@5.2.0 ~mpi
RUN /spack/bin/spack install akantu@master external_solvers=mumps ~mpi ^mumps@5.2.0~mpi
RUN yum install -y environment-modules
RUN git clone https://gitlab.com/akantu/pypi_akantu.git
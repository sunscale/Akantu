MUMPS_TYPE = @MUMPS_TYPE@
PLAT     = @MUMPS_PREFIX@
LIBEXT   = @CMAKE_STATIC_LIBRARY_SUFFIX@
SHLIBEXT = @CMAKE_SHARED_LIBRARY_SUFFIX@
OUTC    = -o 
OUTF    = -o 
ISCOTCH    = -I@SCOTCH_INCLUDE_DIR@
LSCOTCH    = @MUMPS_SCOTCH_LIBRARIES@
LPORDDIR = $(topdir)/PORD/lib/
IPORD    = -I$(topdir)/PORD/include/
LPORD    = -L$(LPORDDIR) -lpord$(PLAT)
ORDERINGSF  = -Dpord -Dscotch
ORDERINGSC  = $(ORDERINGSF)
LORDERINGS = $(LMETIS) $(LPORD) $(LSCOTCH)
IORDERINGSF = $(ISCOTCH)
IORDERINGSC = $(IMETIS) $(IPORD) $(ISCOTCH)
RM      = '@CMAKE_COMMAND@' -E remove -f
MKDIR   = '@CMAKE_COMMAND@' -E make_directory
CP      = cp -af

ifeq ($(MUMPS_TYPE),seq)
# CC : C compiler
CC      = '@CMAKE_C_COMPILER@'
# FC : Fortran 90 compiler
FC      = '@CMAKE_Fortran_COMPILER@'
# FL : Fortran linker
FL      = '@CMAKE_Fortran_COMPILER@'
else
# CC : C compiler
CC      = '@MPI_C_COMPILER@'
# FC : Fortran 90 compiler
FC      = '@MPI_Fortran_COMPILER@'
# FL : Fortran linker
FL      = '@MPI_Fortran_COMPILER@'
endif

# AR : Archive object in a library
AR      = '@CMAKE_AR@' vr 
# RANLIB : generate index of an archive file
RANLIB  = '@CMAKE_RANLIB@'
# SCALAP should define the SCALAPACK and  BLACS libraries.
SCALAP  = '@SCALAPACK_LIBRARIES@'

INCPAR  =
LIBPAR  = $(SCALAP)

# The parallel version is not concerned by the next two lines.
# They are related to the sequential library provided by MUMPS,
# to use instead of ScaLAPACK and MPI.
INCSEQ  = -I$(topdir)/libseq
LIBSEQ  = -L$(topdir)/libseq -lmpiseq$(PLAT)
LIBBLAS = @MUMPS_BLAS_LIBRARIES@
LIBOTHERS = @AKANTU_MUMPS_PTHREAD@

# FORTRAN/C COMPATIBILITY:
#  Use:
#    -DAdd_ if your Fortran compiler adds an underscore at the end
#              of symbols,
#     -DAdd__ if your Fortran compiler adds 2 underscores,
#
#     -DUPPER if your Fortran compiler uses uppercase symbols
#
#     leave empty if your Fortran compiler does not change the symbols.
#
CDEFS = @AKANTU_MUMPS_CDEFS@
#COMPILER OPTIONS
OPTF    = -O -w -fPIC
OPTC    = -O -I. -fPIC
OPTL    = -O -fPIC
#Sequential:
ifeq ($(MUMPS_TYPE),seq)
INCS = $(INCSEQ)
LIBS = $(LIBSEQ)
LIBSEQNEEDED = libseqneeded
endif

#Parallel:
ifeq ($(MUMPS_TYPE),par)
INCS = $(INCPAR)
LIBS = $(LIBPAR)
LIBSEQNEEDED =
endif

prefix  = ${INSTALL_DIR}
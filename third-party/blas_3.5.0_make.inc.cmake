####################################################################
#  BLAS make include file.                                         #
#  March 2007                                                      #
####################################################################
#
SHELL = @CMAKE_SH@
#
#  The machine (platform) identifier to append to the library names
#
PLAT = _@CMAKE_SYSTEM_NAME@
#  
#  Modify the FORTRAN and OPTS definitions to refer to the
#  compiler and desired compiler options for your machine.  NOOPT
#  refers to the compiler options desired when NO OPTIMIZATION is
#  selected.  Define LOADER and LOADOPTS to refer to the loader and 
#  desired load options for your machine.
#
FORTRAN  = '@CMAKE_Fortran_COMPILER@'
OPTS     = -O3
DRVOPTS  = $(OPTS)
NOOPT    =
LOADER   = '@CMAKE_Fortran_COMPILER@'
LOADOPTS =
#
#  The archiver and the flag(s) to use when building archive (library)
#  If you system has no ranlib, set RANLIB = echo.
#
ARCH     = '@CMAKE_AR@'
ARCHFLAGS= cr
#
#  The location and name of the Reference BLAS library.
#
BLASLIB      = blas$(PLAT).a

RANLIB   = '@CMAKE_Fortran_COMPILER@' @CMAKE_SHARED_LIBRARY_Fortran_FLAGS@ -shared \
             -o @CMAKE_SHARED_LIBRARY_PREFIX@blas@CMAKE_SHARED_LIBRARY_SUFFIX@ \
             -Wl,--whole-archive $(BLASLIB) -Wl,--no-whole-archive; \
             '@CMAKE_COMMAND@' -E echo

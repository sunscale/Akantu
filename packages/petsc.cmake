#===============================================================================
# @file   petsc.cmake
#
# @author Alejandro M. Aragón <alejandro.aragon@epfl.ch>
# @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Mon Nov 21 2011
# @date last modification: Tue Jan 19 2016
#
# @brief  package description for PETSc support
#
# @section LICENSE
#
# Copyright (©)  2010-2012, 2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de
# Lausanne)  Laboratory (LSMS  -  Laboratoire de  Simulation  en Mécanique  des
# Solides)
#
# Akantu is free  software: you can redistribute it and/or  modify it under the
# terms  of the  GNU Lesser  General Public  License as  published by  the Free
# Software Foundation, either version 3 of the License, or (at your option) any
# later version.
#
# Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
# details.
#
# You should  have received  a copy  of the GNU  Lesser General  Public License
# along with Akantu. If not, see <http://www.gnu.org/licenses/>.
#
#===============================================================================
package_declare(PETSc EXTERNAL
  DESCRIPTION "Add PETSc support in akantu"
  EXTRA_PACKAGE_OPTIONS ARGS "VERSION;3.5"
  DEPENDS parallel)

package_declare_sources(petsc
  model/common/dof_manager/dof_manager_petsc.cc
  model/common/dof_manager/dof_manager_petsc.hh
  model/common/non_linear_solver/non_linear_solver_petsc.cc
  model/common/non_linear_solver/non_linear_solver_petsc.hh
  solver/petsc_wrapper.hh
  solver/solver_petsc.cc
  solver/solver_petsc.hh
  solver/solver_vector_petsc.cc
  solver/solver_vector_petsc.hh
  solver/sparse_matrix_petsc.cc
  solver/sparse_matrix_petsc.hh
  )

package_declare_extra_files_to_package(PETSc
  PROJECT cmake/Modules/FindPETSc.cmake
  )

package_declare_documentation(PETSc
  "This package enables PETSc as a solver in Akantu"
  ""
  "Under Ubuntu (14.04 LTS) the installation can be performed using the commands:"
  "\\begin{command}"
  "  > sudo apt-get install libpetsc3.4.2-dev"
  "\\end{command}"
  ""
)

package_get_option_name(PETSc _opt_name)
if(${_opt_name})
  include(CheckTypeSize)

  package_get_include_dir(PETSc _petsc_include_dir)
  if(_petsc_include_dir)
      package_get_include_dir(MPI _mpi_include_dir)
    set(CMAKE_EXTRA_INCLUDE_FILES petscsys.h)
    set(CMAKE_REQUIRED_INCLUDES ${_petsc_include_dir} ${_mpi_include_dir})
    check_type_size("PetscInt" PETSC_INT_SIZE)

    if(PETSC_INT_SIZE AND NOT PETSC_INT_SIZE EQUAL AKANTU_INTEGER_SIZE)
      message(SEND_ERROR "This version ofma PETSc cannot be used, it is compiled with the wrong size for PetscInt.")
    endif()

    check_type_size("PetscReal" PETSC_REAL_SIZE)
    if(PETSC_REAL_SIZE AND NOT PETSC_REAL_SIZE EQUAL AKANTU_FLOAT_SIZE)
      message(SEND_ERROR "This version of PETSc cannot be used, it is compiled with the wrong size for PetscInt.")
    endif()
  endif()
endif()


package_set_package_system_dependency(PETSc deb libpetsc3.4.2)
package_set_package_system_dependency(PETSc deb-src libpetsc3.4.2-dev)

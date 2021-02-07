/**
 * @file   aka_config.hh.in
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Sun Sep 26 2010
 * @date last modification: Thu Jan 25 2018
 *
 * @brief  Compilation time configuration of Akantu
 *
 *
 * Copyright (©)  2010-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne) Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the terms  of the  GNU Lesser  General Public  License as published by  the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#ifndef AKANTU_AKA_CONFIG_HH_
#define AKANTU_AKA_CONFIG_HH_

#define AKANTU_VERSION_MAJOR 4
#define AKANTU_VERSION_MINOR 0
#define AKANTU_VERSION_PATCH 0
#define AKANTU_VERSION (AKANTU_VERSION_MAJOR * 10000 \
                        + AKANTU_VERSION_MINOR * 100 \
                        + AKANTU_VERSION_PATCH)


namespace akantu {
using Real = double;
using Int = int;
using UInt = unsigned int;
} // akantu

#define AKANTU_INTEGER_SIZE 4
#define AKANTU_FLOAT_SIZE 8

/* #undef AKANTU_HAS_STRDUP */

/* #undef AKANTU_USE_BLAS */
#define AKANTU_USE_LAPACK

#define AKANTU_PARALLEL
#define AKANTU_USE_MPI

#define AKANTU_USE_SCOTCH
/* #undef AKANTU_USE_PTSCOTCH */
/* #undef AKANTU_SCOTCH_NO_EXTERN */

#define AKANTU_IMPLICIT
#define AKANTU_USE_MUMPS
/* #undef AKANTU_USE_PETSC */

#define AKANTU_USE_IOHELPER
/* #undef AKANTU_USE_QVIEW */
/* #undef AKANTU_USE_BLACKDYNAMITE */

#define AKANTU_USE_PYBIND11

/* #undef AKANTU_USE_OBSOLETE_GETTIMEOFDAY */

/* #undef AKANTU_EXTRA_MATERIALS */
/* #undef AKANTU_STUDENTS_EXTRA_PACKAGE */
#define AKANTU_DAMAGE_NON_LOCAL

#define AKANTU_SOLID_MECHANICS
#define AKANTU_STRUCTURAL_MECHANICS
#define AKANTU_HEAT_TRANSFER
#define AKANTU_PYTHON_INTERFACE

#define AKANTU_COHESIVE_ELEMENT
/* #undef AKANTU_PARALLEL_COHESIVE_ELEMENT */

/* #undef AKANTU_IGFEM */

/* #undef AKANTU_USE_CGAL */
/* #undef AKANTU_EMBEDDED */

// Debug tools
/* #undef AKANTU_NDEBUG */
/* #undef AKANTU_DEBUG_TOOLS */
#define READLINK_COMMAND /bin/readlink
#define ADDR2LINE_COMMAND /usr/bin/addr2line

#endif /* AKANTU_AKA_CONFIG_HH_ */

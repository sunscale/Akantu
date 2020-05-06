/**
 * @file   aka_warning_restore.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Mon Nov 14 2016
 * @date last modification: Fri Dec 02 2016
 *
 * @brief  file to include to reactivate the previously deactivatied warnings
 *
 *
 * Copyright (©) 2016-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

#if defined(__INTEL_COMPILER)
//#pragma warning ( disable : 383 )
#elif defined(__clang__) // test clang to be sure that when we test for gnu it
// is only gnu
#pragma clang diagnostic pop
#elif defined(__GNUG__)
#if GCC_VERSION > 40600
#pragma GCC diagnostic pop
#else
#if defined(AKANTU_WARNING_IGNORE_UNUSED_PARAMETER)
#pragma GCC diagnostic warning "-Wunused-parameter"
#endif
#endif
#endif

#undef AKANTU_WARNING_IGNORE_UNUSED_PARAMETER

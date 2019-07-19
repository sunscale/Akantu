/**
 * @file   aka_voigthelper.cc
 *
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 * @author Till Junge <till.junge@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Dec 20 2013
 * @date last modification: Mon Jun 19 2017
 *
 * @brief  Voigt indices
 *
 * @section LICENSE
 *
 * Copyright (©) 2014-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

/* -------------------------------------------------------------------------- */
#include "aka_voigthelper.hh"
#include "aka_common.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
/* clang-format off */
template <> const UInt VoigtHelper<1>::mat[][1] = {{0}};
template <> const UInt VoigtHelper<2>::mat[][2] = {{0, 2},
                                                   {3, 1}};
template <> const UInt VoigtHelper<3>::mat[][3] = {{0, 5, 4},
                                                   {8, 1, 3},
                                                   {7, 6, 2}};
template <> const UInt VoigtHelper<1>::vec[][2] = {{0, 0}};
template <> const UInt VoigtHelper<2>::vec[][2] = {{0, 0},
                                                   {1, 1},
                                                   {0, 1},
                                                   {1, 0}};
template <> const UInt VoigtHelper<3>::vec[][2] = {{0, 0},
                                                   {1, 1},
                                                   {2, 2},
                                                   {1, 2},
                                                   {0, 2},
                                                   {0, 1},
                                                   {2, 1},
                                                   {2, 0},
                                                   {1, 0}};
template <> const Real VoigtHelper<1>::factors[] = {1.};
template <> const Real VoigtHelper<2>::factors[] = {1., 1., 2.};
template <> const Real VoigtHelper<3>::factors[] = {1., 1., 1.,
                                                    2., 2., 2.};

/* clang-format on */
} // namespace akantu

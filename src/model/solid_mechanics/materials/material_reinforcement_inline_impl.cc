/**
 * @file   material_reinforcement_inline_impl.cc
 *
 * @author Lucas Frérot <lucas.frerot@epfl.ch>
 *
 * @date creation: Thu Mar 12 2015
 * @date last modification: Thu Mar 12 2015
 *
 * @brief  Reinforcement material
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as  published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */

template<UInt d>
inline void MaterialReinforcement<d>::computeDirectingCosinesOnQuad(const Matrix<Real> & nodes,
                                                                 Matrix<Real> & cosines) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_ASSERT(nodes.cols() == 2, "Higher order reinforcement elements not implemented");

  const Vector<Real> & a = nodes(0), b = nodes(1);

  cosines.clear();

  Real sq_length = 0.;
  for (UInt i = 0 ; i < d ; i++) {
    sq_length += Math::pow<2, Real>(b(i) - a(i));
  }

  // Fill the first row of cosine matrix
  for (UInt i = 0 ; i < d ; i++) {
    cosines(0, i) = Math::pow<2, Real>(b(i) - a(i)) / sq_length;
  }

  AKANTU_DEBUG_OUT();
}


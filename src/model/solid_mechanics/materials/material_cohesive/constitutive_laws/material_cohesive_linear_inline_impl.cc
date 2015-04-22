/**
 * @file   material_cohesive_linear_inline_impl.cc
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 * @date   Tue Apr 21 14:49:56 2015
 *
 * @brief  Inline functions of the MaterialCohesiveLinear
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

template<UInt dim>
inline Real MaterialCohesiveLinear<dim>::computeEffectiveNorm(const Matrix<Real> & stress,
							      const Vector<Real> & normal,
							      const Vector<Real> & tangent,
							      Vector<Real> & normal_traction) const {
  normal_traction.mul<false>(stress, normal);

  Real normal_contrib = normal_traction.dot(normal);

  /// in 3D tangential components must be summed
  Real tangent_contrib = 0;

  if (dim == 2) {
    Real tangent_contrib_tmp = normal_traction.dot(tangent);
    tangent_contrib += tangent_contrib_tmp * tangent_contrib_tmp;
  }
  else if (dim == 3) {
    for (UInt s = 0; s < dim - 1; ++s) {
      const Vector<Real> tangent_v(tangent.storage() + s * dim, dim);
      Real tangent_contrib_tmp = normal_traction.dot(tangent_v);
      tangent_contrib += tangent_contrib_tmp * tangent_contrib_tmp;
    }
  }

  tangent_contrib = std::sqrt(tangent_contrib);
  normal_contrib = std::max(0., normal_contrib);

  return std::sqrt(normal_contrib * normal_contrib + tangent_contrib * tangent_contrib * beta2_inv);
}

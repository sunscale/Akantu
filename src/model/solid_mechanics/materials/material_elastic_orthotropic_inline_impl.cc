/**
 * @file   material_elastic_orthotropic_inline_impl.cc
 *
 * @author Enrico Milanese <enrico.milanese@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Thu Feb 15 2018
 * @date last modification: Thu Feb 15 2018
 *
 * @brief Implementation of the inline functions of the material elastic linear
 * orthotropic
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2012, 2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de
 * Lausanne)  Laboratory (LSMS  -  Laboratoire de  Simulation  en Mécanique  des
 * Solides)
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
#include "material_elastic_orthotropic.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MATERIAL_ELASTIC_ORTHOTROPIC_INLINE_IMPL_CC__
#define __AKANTU_MATERIAL_ELASTIC_ORTHOTROPIC_INLINE_IMPL_CC__

namespace akantu {

/* -------------------------------------------------------------------------- */
template <UInt dim>
inline void MaterialElasticOrthotropic<dim>::computePotentialEnergyOnQuad(
    const Matrix<Real> & grad_u, const Matrix<Real> & sigma, Real & epot) {
  epot = .5 * sigma.doubleDot(grad_u);
}

/* -------------------------------------------------------------------------- */
inline Real vector_norm(Vector<Real> & vec) {
  Real norm = 0;
  for (UInt i = 0; i < vec.size(); ++i) {
    norm += vec(i) * vec(i);
  }
  return std::sqrt(norm);
}

} // akantu

#endif /* __AKANTU_MATERIAL_ELASTIC_ORTHOTROPIC_INLINE_IMPL_CC__ */

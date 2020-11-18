/**
 * @file   solid_mechanics_model_tmpl.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Dana Christen <dana.christen@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Tue Nov 07 2017
 *
 * @brief  template part of solid mechanics model
 *
 *
 * Copyright (©)  2010-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "material.hh"
#include "solid_mechanics_model.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_SOLID_MECHANICS_MODEL_TMPL_HH_
#define AKANTU_SOLID_MECHANICS_MODEL_TMPL_HH_

namespace akantu {

#define FWD(...) ::std::forward<decltype(__VA_ARGS__)>(__VA_ARGS__)

/* -------------------------------------------------------------------------- */
template <typename Operation>
void SolidMechanicsModel::splitByMaterial(const Array<Element> & elements,
                                          Operation && op) const {
  std::vector<Array<Element>> elements_per_mat(materials.size());
  this->splitElementByMaterial(elements, elements_per_mat);

  for (auto && mat : zip(materials, elements_per_mat)) {
    FWD(op)(FWD(*std::get<0>(mat)), FWD(std::get<1>(mat)));
  }
}

#undef FWD
/* -------------------------------------------------------------------------- */

} // namespace akantu

#endif /* AKANTU_SOLID_MECHANICS_MODEL_TMPL_HH_ */

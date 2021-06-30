/**
 * @file   structural_element_bernoulli_beam_2.hh
 *
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 * @author Sébastien Hartmann <sebastien.hartmann@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Oct 11 2017
 * @date last modification: Wed Jan 10 2018
 *
 * @brief  Specific functions for bernoulli beam 2d
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

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "structural_mechanics_model.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_STRUCTURAL_ELEMENT_BERNOULLI_BEAM_2_HH_
#define AKANTU_STRUCTURAL_ELEMENT_BERNOULLI_BEAM_2_HH_

namespace akantu {

/* -------------------------------------------------------------------------- */
template <>
void StructuralMechanicsModel::computeTangentModuli<_bernoulli_beam_2>(
    Array<Real> & tangent_moduli) {
  // auto nb_element = getFEEngine().getMesh().getNbElement(_bernoulli_beam_2);
  auto nb_quadrature_points =
      getFEEngine().getNbIntegrationPoints(_bernoulli_beam_2);
  auto tangent_size = 2;

  tangent_moduli.zero();
  auto D_it = tangent_moduli.begin(tangent_size, tangent_size);
  auto el_mat = element_material(_bernoulli_beam_2, _not_ghost).begin();

  for (auto & mat : element_material(_bernoulli_beam_2, _not_ghost)) {
    auto E = materials[mat].E;
    auto A = materials[mat].A;
    auto I = materials[mat].I;
    for (UInt q = 0; q < nb_quadrature_points; ++q, ++D_it) {
      auto & D = *D_it;
      D(0, 0) = E * A;
      D(1, 1) = E * I;
    }
  }
}

} // namespace akantu

#endif /* AKANTU_STRUCTURAL_ELEMENT_BERNOULLI_BEAM_2_HH_ */

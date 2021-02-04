/**
 * @file   structural_element_bernoulli_beam_3.hh
 *
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 * @author Sébastien Hartmann <sebastien.hartmann@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Damien Spielmann <damien.spielmann@epfl.ch>
 *
 * @date creation: Wed Oct 11 2017
 * @date last modification: Tue Feb 20 2018
 *
 * @brief  Specific functions for bernoulli beam 3d
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
#ifndef AKANTU_STRUCTURAL_ELEMENT_BERNOULLI_BEAM_3_HH_
#define AKANTU_STRUCTURAL_ELEMENT_BERNOULLI_BEAM_3_HH_

#include "structural_mechanics_model.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
template <>
void StructuralMechanicsModel::computeTangentModuli<_bernoulli_beam_3>(
    Array<Real> & tangent_moduli) {
  UInt nb_element = getFEEngine().getMesh().getNbElement(_bernoulli_beam_3);
  UInt nb_quadrature_points =
      getFEEngine().getNbIntegrationPoints(_bernoulli_beam_3);
  UInt tangent_size = 4;

  tangent_moduli.zero();
  Array<Real>::matrix_iterator D_it =
      tangent_moduli.begin(tangent_size, tangent_size);

  for (UInt e = 0; e < nb_element; ++e) {
    UInt mat = element_material(_bernoulli_beam_3, _not_ghost)(e);
    Real E = materials[mat].E;
    Real A = materials[mat].A;
    Real Iz = materials[mat].Iz;
    Real Iy = materials[mat].Iy;
    Real GJ = materials[mat].GJ;
    for (UInt q = 0; q < nb_quadrature_points; ++q, ++D_it) {
      Matrix<Real> & D = *D_it;
      D(0, 0) = E * A;
      D(1, 1) = E * Iz;
      D(2, 2) = E * Iy;
      D(3, 3) = GJ;
    }
  }
}

} // namespace akantu

#endif /* AKANTU_STRUCTURAL_ELEMENT_BERNOULLI_BEAM_3_HH_ */

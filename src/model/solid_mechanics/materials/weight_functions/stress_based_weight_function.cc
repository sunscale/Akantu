/**
 * @file   stress_based_weight_function.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Mon Aug 24 2015
 * @date last modification: Tue Feb 20 2018
 *
 * @brief  implementation of the stres based weight function classes
 *
 *
 * Copyright (©) 2015-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "stress_based_weight_function.hh"
namespace akantu {
/* -------------------------------------------------------------------------- */
StressBasedWeightFunction::StressBasedWeightFunction(NonLocalManager & manager)
    : BaseWeightFunction(manager, "stress_based")
// stress_diag("stress_diag", material), selected_stress_diag(NULL),
// stress_base("stress_base", material), selected_stress_base(NULL),
// characteristic_size("lc", material),  selected_characteristic_size(NULL)
{

  // this->registerParam("ft", this->ft, 0., _pat_parsable, "Tensile strength");
  // stress_diag.initialize(spatial_dimension);
  // stress_base.initialize(spatial_dimension * spatial_dimension);
  // characteristic_size.initialize(1);
}

/* -------------------------------------------------------------------------- */
/// During intialization the characteristic sizes for all quadrature
/// points are computed
void StressBasedWeightFunction::init() {
  // const Mesh & mesh = this->material.getModel().getFEEngine().getMesh();
  // for (UInt g = _not_ghost; g <= _ghost; ++g) {
  //   GhostType gt = GhostType(g);
  //   Mesh::type_iterator it = mesh.firstType(spatial_dimension, gt);
  //   Mesh::type_iterator last_type = mesh.lastType(spatial_dimension, gt);
  //   for(; it != last_type; ++it) {
  //     UInt nb_quadrature_points =
  // 	this->material.getModel().getFEEngine().getNbQuadraturePoints(*it, gt);
  //     const Array<UInt> & element_filter =
  //     this->material.getElementFilter(*it, gt);
  //     UInt nb_element = element_filter.size();

  //     Array<Real> ones(nb_element*nb_quadrature_points, 1, 1.);
  //     Array<Real> & lc = characteristic_size(*it, gt);
  //     this->material.getModel().getFEEngine().integrateOnQuadraturePoints(ones,
  // 								     lc,
  // 								     1,
  // 								     *it,
  // 								     gt,
  // 								     element_filter);

  //     for (UInt q = 0;  q < nb_quadrature_points * nb_element; q++) {
  // 	lc(q) = pow(lc(q), 1./ Real(spatial_dimension));
  //     }
  //   }
  // }
}

/* -------------------------------------------------------------------------- */
/// computation of principals stresses and principal directions
void StressBasedWeightFunction::updatePrincipalStress(__attribute__((unused))
                                                      GhostType ghost_type) {
  //   AKANTU_DEBUG_IN();

  //   const Mesh & mesh = this->material.getModel().getFEEngine().getMesh();

  //   Mesh::type_iterator it = mesh.firstType(spatial_dimension, ghost_type);
  //   Mesh::type_iterator last_type = mesh.lastType(spatial_dimension,
  //   ghost_type);
  //   for(; it != last_type; ++it) {
  //     Array<Real>::const_matrix_iterator sigma =
  //       this->material.getStress(*it, ghost_type).begin(spatial_dimension,
  //       spatial_dimension);
  //     auto eigenvalues =
  //       stress_diag(*it, ghost_type).begin(spatial_dimension);
  //     auto eigenvalues_end =
  //       stress_diag(*it, ghost_type).end(spatial_dimension);
  //     Array<Real>::matrix_iterator eigenvector =
  //       stress_base(*it, ghost_type).begin(spatial_dimension,
  //       spatial_dimension);

  // #ifndef __trick__
  //     auto cl = characteristic_size(*it, ghost_type).begin();
  // #endif
  //     UInt q = 0;
  //     for(;eigenvalues != eigenvalues_end; ++sigma, ++eigenvalues,
  //     ++eigenvector, ++cl, ++q) {
  //       sigma->eig(*eigenvalues, *eigenvector);
  //       *eigenvalues /= ft;
  // #ifndef __trick__
  //       // specify a lower bound for principal stress based on the size of
  //       the element
  //       for (UInt i = 0; i < spatial_dimension; ++i) {
  //         (*eigenvalues)(i) = std::max(*cl / this->R, (*eigenvalues)(i));
  //       }
  // #endif
  //     }
  //   }
  //   AKANTU_DEBUG_OUT();
}

} // namespace akantu

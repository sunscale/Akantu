/**
 * @file   non_local_neighborhood.cc
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @date   Sat Sep 26 18:43:30 2015
 *
 * @brief  Implementation of class non-local neighborhood
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
#include "non_local_neighborhood.hh"
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
NonLocalNeighborhood::NonLocalNeighborhood(SolidMechanicsModel & model, Real radius, const ID & type)  :
  NonLocalNeighborhoodBase(model, radius) {

  AKANTU_DEBUG_IN();


  for(UInt gt = _not_ghost; gt <= _ghost; ++gt) {
    GhostType ghost_type = (GhostType) gt;
    pair_weight[ghost_type] = NULL;
  }

  // if (weight_func_type == "base_wf") 
  /// this->weight_function = new Dummy();
  // case damage_wf:
  //   this->weight_function = new DamagedWeightFunction(); break;
  // case remove_wf:
  //   this->weight_function = new RemoveDamagedWeightFunction(); break;
  // case stress_wf:
  //   this->weight_function = new StressBasedWeightFunction(); break


  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
NonLocalNeighborhood::~NonLocalNeighborhood() {
  AKANTU_DEBUG_IN();


  for(UInt gt = _not_ghost; gt <= _ghost; ++gt) {
    GhostType ghost_type = (GhostType) gt;
    delete pair_weight[ghost_type];
  }

  delete weight_function;


  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NonLocalNeighborhood::computeWeights() {
  AKANTU_DEBUG_IN();

  weight_function->updateInternals();

  // for(UInt gt = _not_ghost; gt <= _ghost; ++gt) {
  //   GhostType ghost_type2 = (GhostType) gt;

  //   if(!(pair_weight[ghost_type2])) {
  //     std::string ghost_id = "";
  //     if (ghost_type2 == _ghost) ghost_id = ":ghost";
  //     std::stringstream sstr; sstr << getID() << ":pair_weight:" << ghost_id;
  //     pair_weight[ghost_type2] = new Array<Real>(0, 2, sstr.str());
  //   }

  //   pair_weight[ghost_type2]->resize(pair_list[ghost_type2].size());
  //   pair_weight[ghost_type2]->clear();

  //   PairList::const_iterator first_pair = pair_list[ghost_type2].begin();
  //   PairList::const_iterator last_pair  = pair_list[ghost_type2].end();

  //   Array<Real>::vector_iterator weight_it = pair_weight[ghost_type2]->begin(2);

  //   // Compute the weights
  //   for(;first_pair != last_pair; ++first_pair, ++weight_it) {
  //     Vector<Real> & weight = *weight_it;

  //     const QuadraturePoint & lq1 = first_pair->first;
  //     const QuadraturePoint & lq2 = first_pair->second;

  //     QuadraturePoint gq1 = this->convertToGlobalPoint(lq1);
  //     QuadraturePoint gq2 = this->convertToGlobalPoint(lq2);

  //     //   const Real q2_wJ = fem.getIntegratorInterface().getJacobians(gq2.type, gq2.ghost_type)(gq2.global_num);

  //     Array<Real>::const_vector_iterator quad_coords_1 =
  // 	quadrature_points_coordinates(lq1.type, lq1.ghost_type).begin(spatial_dimension);
  //     Array<Real>::const_vector_iterator quad_coords_2 =
  // 	quadrature_points_coordinates(lq2.type, lq2.ghost_type).begin(spatial_dimension);

  //     Array<Real> & quad_volumes_1 = quadrature_points_volumes(lq1.type, lq1.ghost_type);
  //     const Array<Real> & jacobians_2 =
  // 	fem.getIntegratorInterface().getJacobians(gq2.type, gq2.ghost_type);
  //     const Real & q2_wJ = jacobians_2(gq2.global_num);
  //     // Real & q1_volume = quad_volumes(lq1.global_num);

  //     const Vector<Real> & q1_coord = quad_coords_1[lq1.global_num];
  //     // quadrature_points_coordinates(lq1.type, lq1.ghost_type).begin(spatial_dimension)[lq1.global_num];
  //     const Vector<Real> & q2_coord = quad_coords_2[lq2.global_num];
  //     // quadrature_points_coordinates(lq1.type, lq1.ghost_type).begin(spatial_dimension)[lq2.global_num];

  //     this->weight_func->selectType(lq1.type, lq1.ghost_type, lq2.type, lq2.ghost_type);

  //     // Weight function
  //     Real r = q1_coord.distance(q2_coord);
  //     Real w1 = this->weight_func->operator()(r, lq1, lq2);
  //     weight(0) = q2_wJ * w1;
  //     //     q1_volume += weight(0);
  //     quad_volumes_1(lq1.global_num) += weight(0);

  //     if(lq2.ghost_type != _ghost && lq1.global_num != lq2.global_num) {
  // 	const Array<Real> & jacobians_1 =
  // 	  fem.getIntegratorInterface().getJacobians(gq1.type, gq1.ghost_type);
  // 	Array<Real> & quad_volumes_2 =
  // 	  quadrature_points_volumes(lq2.type, lq2.ghost_type);

  // 	const Real & q1_wJ = jacobians_1(gq1.global_num);
  // 	//Real & q2_volume = quad_volumes_2(lq2.global_num);

  // 	Real w2 = this->weight_func->operator()(r, lq2, lq1);
  // 	weight(1) = q1_wJ * w2;

  // 	quad_volumes_2(lq2.global_num) += weight(1);
  // 	//q2_volume += weight(1);
  //     } else
  // 	weight(1) = 0.;
  //   }
  // }

  // //normalize the weights
  // for(UInt gt = _not_ghost; gt <= _ghost; ++gt) {
  //   GhostType ghost_type2 = (GhostType) gt;

  //   PairList::const_iterator first_pair = pair_list[ghost_type2].begin();
  //   PairList::const_iterator last_pair  = pair_list[ghost_type2].end();

  //   Array<Real>::vector_iterator weight_it = pair_weight[ghost_type2]->begin(2);

  //   // Compute the weights
  //   for(;first_pair != last_pair; ++first_pair, ++weight_it) {
  //     Vector<Real> & weight = *weight_it;

  //     const QuadraturePoint & lq1 = first_pair->first;
  //     const QuadraturePoint & lq2 = first_pair->second;

  //     Array<Real> & quad_volumes_1 = quadrature_points_volumes(lq1.type, lq1.ghost_type);
  //     Array<Real> & quad_volumes_2 = quadrature_points_volumes(lq2.type, lq2.ghost_type);

  //     Real q1_volume = quad_volumes_1(lq1.global_num);

  //     weight(0) *= 1. / q1_volume;
  //     if(ghost_type2 != _ghost) {
  // 	Real q2_volume = quad_volumes_2(lq2.global_num);
  // 	weight(1) *= 1. / q2_volume;
  //     }
  //   }
  // }

  AKANTU_DEBUG_OUT();
}

__END_AKANTU__

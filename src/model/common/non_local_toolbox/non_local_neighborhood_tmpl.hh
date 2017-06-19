/**
 * @file   non_local_neighborhood_tmpl.hh
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Mon Sep 28 2015
 * @date last modification: Wed Nov 25 2015
 *
 * @brief  Implementation of class non-local neighborhood
 *
 * @section LICENSE
 *
 * Copyright (©) 2015 EPFL (Ecole Polytechnique Fédérale de Lausanne) Laboratory
 * (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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
#include "non_local_manager.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_NON_LOCAL_NEIGHBORHOOD_TMPL_HH__
#define __AKANTU_NON_LOCAL_NEIGHBORHOOD_TMPL_HH__

namespace akantu {

/* -------------------------------------------------------------------------- */
template<class WeightFunction>
NonLocalNeighborhood<WeightFunction>::NonLocalNeighborhood(NonLocalManager & manager,
							   const ElementTypeMapReal & quad_coordinates,
							   const ID & id,
							   const MemoryID & memory_id)  :
  NonLocalNeighborhoodBase(manager.getModel(), quad_coordinates, id, memory_id),
  non_local_manager(&manager) {

  AKANTU_DEBUG_IN();


  for(UInt gt = _not_ghost; gt <= _ghost; ++gt) {
    GhostType ghost_type = (GhostType) gt;
    pair_weight[ghost_type] = NULL;
  }

  this->weight_function = new WeightFunction(this->getNonLocalManager());

  this->registerSubSection(_st_weight_function, "weight_parameter", *weight_function);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<class WeightFunction>
NonLocalNeighborhood<WeightFunction>::~NonLocalNeighborhood() {
  AKANTU_DEBUG_IN();


  for(UInt gt = _not_ghost; gt <= _ghost; ++gt) {
    GhostType ghost_type = (GhostType) gt;
    delete pair_weight[ghost_type];
  }

  delete weight_function;


  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<class WeightFunction>
void NonLocalNeighborhood<WeightFunction>::computeWeights() {
  AKANTU_DEBUG_IN();

  this->weight_function->setRadius(this->neighborhood_radius);
  Vector<Real> q1_coord(this->spatial_dimension);
  Vector<Real> q2_coord(this->spatial_dimension);
 
  UInt nb_weights_per_pair = 2; /// w1: q1->q2, w2: q2->q1
 
  /// get the elementtypemap for the neighborhood volume for each quadrature point
  ElementTypeMapReal & quadrature_points_volumes = this->non_local_manager->getVolumes();
 
  /// update the internals of the weight function if applicable (not
  /// all the weight functions have internals and do noting in that
  /// case)
  weight_function->updateInternals();

  for(UInt gt = _not_ghost; gt <= _ghost; ++gt) {
    GhostType ghost_type = (GhostType) gt;

    /// allocate the array to store the weight, if it doesn't exist already
    if(!(pair_weight[ghost_type])) {
      std::string ghost_id = "";
      if (ghost_type == _ghost) ghost_id = ":ghost";
      std::stringstream sstr; sstr << this->id <<":pair_weight:" << ghost_id;
      pair_weight[ghost_type] = new Array<Real>(0, nb_weights_per_pair, sstr.str());
    }


    /// resize the array to the correct size
    pair_weight[ghost_type]->resize(pair_list[ghost_type].size());
    /// set entries to zero
    pair_weight[ghost_type]->clear();

    /// loop over all pairs in the current pair list array and their corresponding weights
    PairList::const_iterator first_pair = pair_list[ghost_type].begin();
    PairList::const_iterator last_pair  = pair_list[ghost_type].end();
    Array<Real>::vector_iterator weight_it = pair_weight[ghost_type]->begin(nb_weights_per_pair);

    // Compute the weights
    for(;first_pair != last_pair; ++first_pair, ++weight_it) {
      Vector<Real> & weight = *weight_it;
      const IntegrationPoint & q1 = first_pair->first;
      const IntegrationPoint & q2 = first_pair->second;

      /// get the coordinates for the given pair of quads
      Array<Real>::const_vector_iterator coords_type_1_it = this->quad_coordinates(q1.type, q1.ghost_type).begin(this->spatial_dimension);
      q1_coord = coords_type_1_it[q1.global_num];
      Array<Real>::const_vector_iterator coords_type_2_it = this->quad_coordinates(q2.type, q2.ghost_type).begin(this->spatial_dimension);
      q2_coord = coords_type_2_it[q2.global_num];

      Array<Real> & quad_volumes_1 = quadrature_points_volumes(q1.type, q1.ghost_type);
      const Array<Real> & jacobians_2 =
	this->non_local_manager->getJacobians(q2.type, q2.ghost_type);
      const Real & q2_wJ = jacobians_2(q2.global_num);

      /// compute distance between the two quadrature points
      Real r = q1_coord.distance(q2_coord);

      /// compute the weight for averaging on q1 based on the distance
      Real w1 = this->weight_function->operator()(r, q1, q2);
      weight(0) = q2_wJ * w1;

      quad_volumes_1(q1.global_num) += weight(0);

      if(q2.ghost_type != _ghost && q1.global_num != q2.global_num) {
	const Array<Real> & jacobians_1 =
	  this->non_local_manager->getJacobians(q1.type, q1.ghost_type);
	Array<Real> & quad_volumes_2 =
	  quadrature_points_volumes(q2.type, q2.ghost_type);
	/// compute the weight for averaging on q2
	const Real & q1_wJ = jacobians_1(q1.global_num);
	Real w2 = this->weight_function->operator()(r, q2, q1);
	weight(1) = q1_wJ * w2;
	quad_volumes_2(q2.global_num) += weight(1);
      } else
	weight(1) = 0.;
    }
  }

  ///  normalize the weights
  for(UInt gt = _not_ghost; gt <= _ghost; ++gt) {
    GhostType ghost_type2 = (GhostType) gt;

    PairList::const_iterator first_pair = pair_list[ghost_type2].begin();
    PairList::const_iterator last_pair  = pair_list[ghost_type2].end();

    Array<Real>::vector_iterator weight_it = pair_weight[ghost_type2]->begin(nb_weights_per_pair);

    // Compute the weights
    for(;first_pair != last_pair; ++first_pair, ++weight_it) {
      Vector<Real> & weight = *weight_it;

      const IntegrationPoint & q1 = first_pair->first;
      const IntegrationPoint & q2 = first_pair->second;

      Array<Real> & quad_volumes_1 = quadrature_points_volumes(q1.type, q1.ghost_type);
      Array<Real> & quad_volumes_2 = quadrature_points_volumes(q2.type, q2.ghost_type);

      Real q1_volume = quad_volumes_1(q1.global_num);

      weight(0) *= 1. / q1_volume;
      if(ghost_type2 != _ghost) {
	Real q2_volume = quad_volumes_2(q2.global_num);
	weight(1) *= 1. / q2_volume;
      }
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<class WeightFunction>
void NonLocalNeighborhood<WeightFunction>::saveWeights(const std::string & filename) const {
  std::ofstream pout;

  std::stringstream sstr;

  StaticCommunicator & comm = StaticCommunicator::getStaticCommunicator();
  Int prank = comm.whoAmI();
  sstr << filename << "." << prank;

  pout.open(sstr.str().c_str());

  for (UInt gt = _not_ghost; gt <= _ghost; ++gt) {
    GhostType ghost_type = (GhostType) gt;
 
    AKANTU_DEBUG_ASSERT((pair_weight[ghost_type]), "the weights have not been computed yet");
   
    Array<Real> & weights = *(pair_weight[ghost_type]);
    Array<Real>::const_vector_iterator weights_it = weights.begin(2);
    for (UInt i = 0; i < weights.getSize(); ++i, ++weights_it) 
      pout << "w1: " << (*weights_it)(0) <<" w2: " << (*weights_it)(1) << std::endl;
  }
}

/* -------------------------------------------------------------------------- */
template<class WeightFunction>
void NonLocalNeighborhood<WeightFunction>::weightedAverageOnNeighbours(const ElementTypeMapReal & to_accumulate,
								       ElementTypeMapReal & accumulated,
								       UInt nb_degree_of_freedom,
								       const GhostType & ghost_type2) const {
  AKANTU_DEBUG_IN();
  std::set<ID>::iterator it = non_local_variables.find(accumulated.getName());
  ///do averaging only for variables registered in the neighborhood
  if (it != non_local_variables.end()) {
    PairList::const_iterator first_pair = pair_list[ghost_type2].begin();
    PairList::const_iterator last_pair  = pair_list[ghost_type2].end();

    Array<Real>::vector_iterator weight_it = pair_weight[ghost_type2]->begin(2);

    // Compute the weights
    for(;first_pair != last_pair; ++first_pair, ++weight_it) {
      Vector<Real> & weight = *weight_it;

      const IntegrationPoint & q1 = first_pair->first;
      const IntegrationPoint & q2 = first_pair->second;

      const Array<Real> & to_acc_1 = to_accumulate(q1.type, q1.ghost_type);
      Array<Real> & acc_1 = accumulated(q1.type, q1.ghost_type);
      const Array<Real> & to_acc_2 = to_accumulate(q2.type, q2.ghost_type);
      Array<Real> & acc_2 = accumulated(q2.type, q2.ghost_type);

      for(UInt d = 0; d < nb_degree_of_freedom; ++d) {
	acc_1(q1.global_num, d) += weight(0) * to_acc_2(q2.global_num, d);
      }

      if(ghost_type2 != _ghost) {
	for(UInt d = 0; d < nb_degree_of_freedom; ++d) {
	  acc_2(q2.global_num, d) += weight(1) * to_acc_1(q1.global_num, d);
	}
      }
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<class WeightFunction>
void NonLocalNeighborhood<WeightFunction>::updateWeights() {

  // Update the weights for the non local variable averaging  
  if(this->weight_function->getUpdateRate() &&
     (this->non_local_manager->getNbStressCalls() % this->weight_function->getUpdateRate() == 0))  {
    this->synch_registry->asynchronousSynchronize(_gst_mnl_weight);
    this->synch_registry->waitEndSynchronize(_gst_mnl_weight);
    this->computeWeights();
  }
}


/* -------------------------------------------------------------------------- */
template<class WeightFunction>
void NonLocalNeighborhood<WeightFunction>::registerNonLocalVariable(const ID & id) {
  this->non_local_variables.insert(id);
}


} // akantu

#endif /* __AKANTU_NON_LOCAL_NEIGHBORHOOD_TMPL__ */

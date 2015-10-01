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
#include "non_local_manager.hh"
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
NonLocalNeighborhood::NonLocalNeighborhood(NonLocalManager & manager, 
					   Real radius, 
					   const ID & type,
					   const ID & id,
					   const MemoryID & memory_id)  :
  NonLocalNeighborhoodBase(manager.getModel(), radius, id, memory_id),
  non_local_manager(&manager),
  type(type) {

  AKANTU_DEBUG_IN();


  for(UInt gt = _not_ghost; gt <= _ghost; ++gt) {
    GhostType ghost_type = (GhostType) gt;
    pair_weight[ghost_type] = NULL;
  }

  // if (weight_func_type == "base_wf") 
  this->weight_function = new WeightFunction();
  this->weight_function->setRadius(radius);

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

      const QuadraturePoint & q1 = first_pair->first;
      const QuadraturePoint & q2 = first_pair->second;

      Array<Real> & quad_volumes_1 = quadrature_points_volumes(q1.type, q1.ghost_type);
      const Array<Real> & jacobians_2 =
  	this->non_local_manager->getJacobians(q2.type, q2.ghost_type);
      const Real & q2_wJ = jacobians_2(q2.global_num);


      /// compute distance between the two quadrature points
      const Vector<Real> & q1_coord = q1.getPosition();
      const Vector<Real> & q2_coord = q2.getPosition();
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

      const QuadraturePoint & q1 = first_pair->first;
      const QuadraturePoint & q2 = first_pair->second;

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
void NonLocalNeighborhood::saveWeights(const std::string & filename) const {
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


__END_AKANTU__

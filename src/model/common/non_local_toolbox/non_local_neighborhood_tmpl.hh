/**
 * @file   non_local_neighborhood_tmpl.hh
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Mon Sep 28 2015
 * @date last modification: Tue Feb 20 2018
 *
 * @brief  Implementation of class non-local neighborhood
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
#include "communicator.hh"
#include "non_local_manager.hh"
#include "non_local_neighborhood.hh"
/* -------------------------------------------------------------------------- */
#include <fstream>
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_NON_LOCAL_NEIGHBORHOOD_TMPL_HH_
#define AKANTU_NON_LOCAL_NEIGHBORHOOD_TMPL_HH_

namespace akantu {

/* -------------------------------------------------------------------------- */
template <class WeightFunction>
template <class Func>
inline void NonLocalNeighborhood<WeightFunction>::foreach_weight(
    GhostType ghost_type, Func && func) {
  auto weight_it =
      pair_weight[ghost_type]->begin(pair_weight[ghost_type]->getNbComponent());

  for (auto & pair : pair_list[ghost_type]) {
    std::forward<decltype(func)>(func)(pair.first, pair.second, *weight_it);
    ++weight_it;
  }
}

/* -------------------------------------------------------------------------- */
template <class WeightFunction>
template <class Func>
inline void NonLocalNeighborhood<WeightFunction>::foreach_weight(
    GhostType ghost_type, Func && func) const {
  auto weight_it =
      pair_weight[ghost_type]->begin(pair_weight[ghost_type]->getNbComponent());

  for (auto & pair : pair_list[ghost_type]) {
    std::forward<decltype(func)>(func)(pair.first, pair.second, *weight_it);
    ++weight_it;
  }
}

/* -------------------------------------------------------------------------- */
template <class WeightFunction>
NonLocalNeighborhood<WeightFunction>::NonLocalNeighborhood(
    NonLocalManager & manager, const ElementTypeMapReal & quad_coordinates,
    const ID & id)
    : NonLocalNeighborhoodBase(manager.getModel(), quad_coordinates, id),
      non_local_manager(manager) {
  AKANTU_DEBUG_IN();

  this->weight_function = std::make_unique<WeightFunction>(manager);

  this->registerSubSection(ParserType::_weight_function, "weight_parameter",
                           *weight_function);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class WeightFunction>
NonLocalNeighborhood<WeightFunction>::~NonLocalNeighborhood() = default;

/* -------------------------------------------------------------------------- */
template <class WeightFunction>
void NonLocalNeighborhood<WeightFunction>::computeWeights() {
  AKANTU_DEBUG_IN();

  this->weight_function->setRadius(this->neighborhood_radius);
  Vector<Real> q1_coord(this->spatial_dimension);
  Vector<Real> q2_coord(this->spatial_dimension);

  UInt nb_weights_per_pair = 2; /// w1: q1->q2, w2: q2->q1

  /// get the elementtypemap for the neighborhood volume for each quadrature
  /// point
  ElementTypeMapReal & quadrature_points_volumes =
      this->non_local_manager.getVolumes();

  /// update the internals of the weight function if applicable (not
  /// all the weight functions have internals and do noting in that
  /// case)
  weight_function->updateInternals();

  for (auto ghost_type : ghost_types) {
    /// allocate the array to store the weight, if it doesn't exist already
    if (!(pair_weight[ghost_type])) {
      pair_weight[ghost_type] =
          std::make_unique<Array<Real>>(0, nb_weights_per_pair);
    }

    /// resize the array to the correct size
    pair_weight[ghost_type]->resize(pair_list[ghost_type].size());
    /// set entries to zero
    pair_weight[ghost_type]->zero();

    /// loop over all pairs in the current pair list array and their
    /// corresponding weights
    auto first_pair = pair_list[ghost_type].begin();
    auto last_pair = pair_list[ghost_type].end();
    auto weight_it = pair_weight[ghost_type]->begin(nb_weights_per_pair);

    // Compute the weights
    for (; first_pair != last_pair; ++first_pair, ++weight_it) {
      Vector<Real> & weight = *weight_it;
      const IntegrationPoint & q1 = first_pair->first;
      const IntegrationPoint & q2 = first_pair->second;

      /// get the coordinates for the given pair of quads
      auto coords_type_1_it = this->quad_coordinates(q1.type, q1.ghost_type)
                                  .begin(this->spatial_dimension);
      q1_coord = coords_type_1_it[q1.global_num];
      auto coords_type_2_it = this->quad_coordinates(q2.type, q2.ghost_type)
                                  .begin(this->spatial_dimension);
      q2_coord = coords_type_2_it[q2.global_num];

      Array<Real> & quad_volumes_1 =
          quadrature_points_volumes(q1.type, q1.ghost_type);
      const Array<Real> & jacobians_2 =
          this->non_local_manager.getJacobians(q2.type, q2.ghost_type);
      const Real & q2_wJ = jacobians_2(q2.global_num);

      /// compute distance between the two quadrature points
      Real r = q1_coord.distance(q2_coord);

      /// compute the weight for averaging on q1 based on the distance
      Real w1 = this->weight_function->operator()(r, q1, q2);
      weight(0) = q2_wJ * w1;

      quad_volumes_1(q1.global_num) += weight(0);

      if (q2.ghost_type != _ghost && q1.global_num != q2.global_num) {
        const Array<Real> & jacobians_1 =
            this->non_local_manager.getJacobians(q1.type, q1.ghost_type);
        Array<Real> & quad_volumes_2 =
            quadrature_points_volumes(q2.type, q2.ghost_type);
        /// compute the weight for averaging on q2
        const Real & q1_wJ = jacobians_1(q1.global_num);
        Real w2 = this->weight_function->operator()(r, q2, q1);
        weight(1) = q1_wJ * w2;
        quad_volumes_2(q2.global_num) += weight(1);
      } else {
        weight(1) = 0.;
      }
    }
  }

  ///  normalize the weights
  for (auto ghost_type : ghost_types) {
    foreach_weight(ghost_type, [&](const auto & q1, const auto & q2,
                                   auto & weight) {
      auto & quad_volumes_1 = quadrature_points_volumes(q1.type, q1.ghost_type);
      auto & quad_volumes_2 = quadrature_points_volumes(q2.type, q2.ghost_type);

      Real q1_volume = quad_volumes_1(q1.global_num);
      auto ghost_type2 = q2.ghost_type;
      weight(0) *= 1. / q1_volume;
      if (ghost_type2 != _ghost) {
        Real q2_volume = quad_volumes_2(q2.global_num);
        weight(1) *= 1. / q2_volume;
      }
    });
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class WeightFunction>
void NonLocalNeighborhood<WeightFunction>::saveWeights(
    const std::string & filename) const {
  std::ofstream pout;

  std::stringstream sstr;

  const Communicator & comm = model.getMesh().getCommunicator();

  Int prank = comm.whoAmI();
  sstr << filename << "." << prank;

  pout.open(sstr.str().c_str());

  for (UInt gt = _not_ghost; gt <= _ghost; ++gt) {
    auto ghost_type = (GhostType)gt;

    AKANTU_DEBUG_ASSERT((pair_weight[ghost_type]),
                        "the weights have not been computed yet");

    Array<Real> & weights = *(pair_weight[ghost_type]);
    auto weights_it = weights.begin(2);
    for (UInt i = 0; i < weights.size(); ++i, ++weights_it) {
      pout << "w1: " << (*weights_it)(0) << " w2: " << (*weights_it)(1)
           << std::endl;
    }
  }
}

/* -------------------------------------------------------------------------- */
template <class WeightFunction>
void NonLocalNeighborhood<WeightFunction>::weightedAverageOnNeighbours(
    const ElementTypeMapReal & to_accumulate, ElementTypeMapReal & accumulated,
    UInt nb_degree_of_freedom, GhostType ghost_type2) const {

  auto it = non_local_variables.find(accumulated.getName());
  // do averaging only for variables registered in the neighborhood
  if (it == non_local_variables.end()) {
    return;
  }

  foreach_weight(
      ghost_type2,
      [ghost_type2, nb_degree_of_freedom, &to_accumulate,
       &accumulated](const auto & q1, const auto & q2, auto & weight) {
        const Vector<Real> to_acc_1 =
            to_accumulate(q1.type, q1.ghost_type)
                .begin(nb_degree_of_freedom)[q1.global_num];
        const Vector<Real> to_acc_2 =
            to_accumulate(q2.type, q2.ghost_type)
                .begin(nb_degree_of_freedom)[q2.global_num];
        Vector<Real> acc_1 = accumulated(q1.type, q1.ghost_type)
                                 .begin(nb_degree_of_freedom)[q1.global_num];
        Vector<Real> acc_2 = accumulated(q2.type, q2.ghost_type)
                                 .begin(nb_degree_of_freedom)[q2.global_num];

        acc_1 += weight(0) * to_acc_2;

        if (ghost_type2 != _ghost) {
          acc_2 += weight(1) * to_acc_1;
        }
      });
}

/* -------------------------------------------------------------------------- */
template <class WeightFunction>
void NonLocalNeighborhood<WeightFunction>::updateWeights() {
  // Update the weights for the non local variable averaging
  if (this->weight_function->getUpdateRate() &&
      (this->non_local_manager.getNbStressCalls() %
           this->weight_function->getUpdateRate() ==
       0)) {
    SynchronizerRegistry::synchronize(SynchronizationTag::_mnl_weight);
    this->computeWeights();
  }
}

} // namespace akantu

#endif /* __AKANTU_NON_LOCAL_NEIGHBORHOOD_TMPL__ */

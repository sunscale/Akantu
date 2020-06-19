/**
 * @file   material_damage_iterative.cc
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 *
 *
 * @brief  Implementation of MaterialDamageIterativeNonLocal
 *
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

/* -------------------------------------------------------------------------- */
#include "material_damage_iterative_non_local.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialDamageIterativeNonLocal<
    spatial_dimension>::computeNonLocalStresses(GhostType ghost_type) {
  AKANTU_DEBUG_IN();
  /// reset normalized maximum equivalent stress
  if (ghost_type == _not_ghost)
    this->norm_max_equivalent_stress = 0;

  MaterialDamageIterativeNonLocalParent::computeNonLocalStresses(ghost_type);

  /// find global Gauss point with highest stress
  const auto & comm = this->model.getMesh().getCommunicator();
  comm.allReduce(this->norm_max_equivalent_stress, SynchronizerOperation::_max);

  AKANTU_DEBUG_OUT();
}
/* -------------------------------------------------------------------------- */
INSTANTIATE_MATERIAL(damage_iterative_non_local,
                     MaterialDamageIterativeNonLocal);

} // namespace akantu

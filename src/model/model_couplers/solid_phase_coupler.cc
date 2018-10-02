/**
 * @file   solid_phase_coupler.cc
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 * 
 * @date creation: Fri Sep 28 2018
 * @date last modification: Fri Sep 28 2018
 *
 * @brief  class for coupling of solid mechancis and phasefield model 
 *
 * @section LICENSE
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
#include "solid_phase_coupler.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

template<SolidMechanicsModel & smm, PhaseFieldModel & pfm>
SolidPhaseCoupler<smm, pfm>::SolidPhaseCoupler() {
  this->spatial_dimension = smm.getMesh().getSpatialDimension();
}

/* -------------------------------------------------------------------------- */
template<SolidMechanicsModel & smm, PhaseFieldModel & pfm>
SolidPhaseCoupler<smm, pfm>::~SolidPhaseCoupler() {
}

/* -------------------------------------------------------------------------- */
template<SolidMechanicsModel & smm, PhaseFieldModel & pfm>
void SolidPhaseCoupler<smm, pfm>::computeDamageOnQuadPoints(const GhostType & ghost_type) {
  AKANTU_DEBUG_IN();
  
  auto & fem  = pfm.getFEEngine();
  auto & mesh = pfm.getMesh();

  for (auto & type: mesh.elementTypes(this->spatial_dimension, ghost_type)) {
    auto & damage_on_qpoints_vect = smm.getMaterial(0).getArray<Real>("damage", type);
    fem.interpolateOnIntegrationPoints(pfm.getDamage(), damage_on_qpoints_vect,
				       1, type, ghost_type); 
  }
  
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<SolidMechanicsModel & smm,PhaseFieldModel & pfm>
void SolidPhaseCoupler<smm, pfm>::computeStrainOnQuadPoints(const GhostType & ghost_type) {
  AKANTU_DEBUG_IN();

  auto & fem  = smm.getFEEngine();
  auto & mesh = smm.getMesh();

  auto & strain_on_qpoints = pfm.getStrain();
  auto & gradu_on_qpoints  = smm.getMaterial(0).getGradU();
    
  for (auto & type: mesh.elementTypes(spatial_dimension, ghost_type)) {
    auto & strain_on_qpoints_vect = strain_on_qpoints(type, ghost_type);
    auto & gradu_on_qpoints_vect  = gradu_on_qpoints(type, ghost_type);
    for (auto && values:
	   zip(make_view(strain_on_qpoints_vect, this->spatial_dimension, this->spatial_dimension),
	       make_view(gradu_on_qpoints_vect,  this->spatial_dimension, this->spatial_dimension))) {
      auto & strain = std::get<0>(values);
      auto & grad_u = std::get<1>(values);
      this->gradUToEpsilon(grad_u, strain);
    }      
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<SolidMechanicsModel & smm, PhaseFieldModel & pfm>
void SolidPhaseCoupler<smm, pfm>::solve() {
  smm.solveStep();
  this->computeStrainOnQuadPoints();

  pfm.solveStep();
  this->computeDamageOnQuadPoints();
  
  smm.updateResidual();
  // check for convergence();

}

/* -------------------------------------------------------------------------- */
template<SolidMechanicsModel & smm, PhaseFieldModel & pfm>
void SolidPhaseCoupler<smm, pfm>::gradUToEpsilon(const Matrix<Real> & grad_u,
						 Matrix<Real> & epsilon) {
  for (UInt i=0; i < this->spatial_dimension; ++i) {
    for (UInt j = 0; j < this->spatial_dimension; ++j)
	epsilon(i, j) = 0.5 * (grad_u(i, j) + grad_u(j, i));
  }
}


} // akantu

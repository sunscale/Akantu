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
#include "non_linear_solver.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

template<typename SolidType, typename PhaseType>
SolidPhaseCoupler<SolidType, PhaseType>::SolidPhaseCoupler(SolidType & solid, PhaseType & phase)
  : solid(solid), phase(phase) {
  this->spatial_dimension = solid.getMesh().getSpatialDimension();
}

/* -------------------------------------------------------------------------- */
template<typename SolidType, typename PhaseType>
SolidPhaseCoupler<SolidType, PhaseType>::~SolidPhaseCoupler() {
}

/* -------------------------------------------------------------------------- */
template<typename SolidType, typename PhaseType>
void SolidPhaseCoupler<SolidType, PhaseType>::computeDamageOnQuadPoints(const GhostType & ghost_type) {
  AKANTU_DEBUG_IN();
  
  auto & fem  = phase.getFEEngine();
  auto & mesh = phase.getMesh();

  switch (spatial_dimension) {
  case 1: {
    auto & mat = static_cast<MaterialPhaseField<1> &>(solid.getMaterial(0));
    auto & damage = mat.getDamage();
    for (auto & type: mesh.elementTypes(this->spatial_dimension, ghost_type)) {
      auto & damage_on_qpoints_vect = damage(type, ghost_type);
      fem.interpolateOnIntegrationPoints(phase.getDamage(), damage_on_qpoints_vect,
				       1, type, ghost_type); 
    }
    
    break;
  }
  case 2: {
    auto & mat = static_cast<MaterialPhaseField<2> &>(solid.getMaterial(0));
    auto & damage = mat.getDamage();
  
    for (auto & type: mesh.elementTypes(this->spatial_dimension, ghost_type)) {
      auto & damage_on_qpoints_vect = damage(type, ghost_type);
      fem.interpolateOnIntegrationPoints(phase.getDamage(), damage_on_qpoints_vect,
				       1, type, ghost_type); 
    }
    break;
  }  
  default:
    auto & mat = static_cast<MaterialPhaseField<3> &>(solid.getMaterial(0));
    break;
  }
  
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<typename SolidType,typename PhaseType>
void SolidPhaseCoupler<SolidType, PhaseType>::computeStrainOnQuadPoints(const GhostType & ghost_type) {
  AKANTU_DEBUG_IN();

  auto & mesh = solid.getMesh();

  auto & strain_on_qpoints = phase.getStrain();
  auto & gradu_on_qpoints  = solid.getMaterial(0).getGradU();
    
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
template<typename SolidType, typename PhaseType>
void SolidPhaseCoupler<SolidType, PhaseType>::solve() {

  this->convergence = true;
  UInt iter = 0;
  UInt max_iter = 10;
  
  while(iter < max_iter) {

    auto u_old = solid.getDisplacement();
    auto d_old = phase.getDamage();

    std::cerr << "---- solving solid model ------  \n";
    
    solid.solveStep();
    this->computeStrainOnQuadPoints(_not_ghost);

    std::cerr << "---- solving phasefield model ------  \n";
    
    phase.solveStep();
    this->computeDamageOnQuadPoints(_not_ghost);

    auto u_new = solid.getDisplacement();
    auto d_new = phase.getDamage();

    //auto u_new = solid.getDOFManager().getDOFs("displacement");
    //auto u_old = solid.getDOFManager().getPreviousDOFs("displacement");
    
    //auto d_new = phase.getDOFManager().getDOFs("damage");
    //auto d_old = phase.getDOFManager().getPreviousDOFs("damage");
   
    this->checkConvergence(u_new, u_old, d_new, d_old);
    
    if (this->convergence) {
      break;
    }

    iter++;
  }
  
}

/* -------------------------------------------------------------------------- */
template<typename SolidType, typename PhaseType>
void SolidPhaseCoupler<SolidType, PhaseType>::gradUToEpsilon(const Matrix<Real> & grad_u,
						     Matrix<Real> & epsilon) {
  for (UInt i=0; i < this->spatial_dimension; ++i) {
    for (UInt j = 0; j < this->spatial_dimension; ++j)
	epsilon(i, j) = 0.5 * (grad_u(i, j) + grad_u(j, i));
  }
}

/* -------------------------------------------------------------------------- */
template<typename SolidType, typename PhaseType>
bool SolidPhaseCoupler<SolidType, PhaseType>::checkConvergence(Array<Real> & u_new, Array<Real> & u_old, Array<Real> & d_new, Array<Real> & d_old) {

  const Array<bool> & blocked_dofs = solid.getBlockedDOFs(); 
  UInt nb_degree_of_freedom = u_new.size();

  auto u_n_it = u_new.begin();
  auto u_o_it = u_old.begin();
  auto bld_it = blocked_dofs.begin();

  
  Real norm = 0;
  for (UInt n = 0; n < nb_degree_of_freedom; ++n, ++u_n_it, ++u_o_it, ++bld_it) {
    if ((!*bld_it)) {
      norm += (*u_n_it - *u_o_it) * (*u_n_it - *u_o_it);
    }
  }

  norm = std::sqrt(norm);

  auto d_n_it = d_new.begin();
  auto d_o_it = d_old.begin();
  nb_degree_of_freedom = d_new.size();

  Real norm2 = 0;
  for (UInt i = 0; i < nb_degree_of_freedom; ++i) {
    norm2 += (*d_n_it - *d_o_it);
  }

  norm2 =std::sqrt(norm2);

  Real error  = std::max(norm, norm2);
  
  Real tolerance = 1e-8;
  if (error < tolerance) {
    
    this->convergence = true;
  }

}

 
template class SolidPhaseCoupler<SolidMechanicsModel, PhaseFieldModel>;
  
} // akantu


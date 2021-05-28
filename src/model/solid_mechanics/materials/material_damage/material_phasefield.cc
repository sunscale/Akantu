/**
 * @file   material_phasefield.cc
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Tue Oct 02 2018
 * @date last modification: Tue Oct 02 2018
 *
 * @brief  Specialization of the material class for the phasefield material
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
#include "material_phasefield.hh"
#include "solid_mechanics_model.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
MaterialPhaseField<spatial_dimension>::MaterialPhaseField(SolidMechanicsModel & model,
							  const ID & id)
  : Parent(model, id) {

  AKANTU_DEBUG_IN();

  this->registerParam("eta", eta, Real(0.), _pat_parsable, "eta");
  this->damage.initialize(0);

  AKANTU_DEBUG_OUT();
}



/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialPhaseField<spatial_dimension>::computeStress(ElementType el_type,
							  GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  auto dam = this->damage(el_type, ghost_type).begin();
 
  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type);

  computeStressOnQuad(grad_u, sigma, *dam);

  ++dam;
  
  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialPhaseField<spatial_dimension>::computeTangentModuli(
    ElementType el_type, Array<Real> & tangent_matrix,
    GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Parent::computeTangentModuli(el_type, tangent_matrix,
			       ghost_type);

  Real * dam = this->damage(el_type, ghost_type).storage();

  MATERIAL_TANGENT_QUADRATURE_POINT_LOOP_BEGIN(tangent_matrix);
  computeTangentModuliOnQuad(tangent, *dam);

  ++dam;

  MATERIAL_TANGENT_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}
  
/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialPhaseField<spatial_dimension>::computeTangentModuliOnQuad(
    Matrix<Real> & tangent, Real & dam) {
  tangent *= (1 - dam)*(1 - dam) + eta;
}

  
INSTANTIATE_MATERIAL(phasefield, MaterialPhaseField);
  

} // akantu



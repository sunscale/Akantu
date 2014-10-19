/**
 * @file   material_linear_isotropic_hardening.cc
 *
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 * @author Daniel Pino Muñoz <daniel.pinomunoz@epfl.ch>
 * @author Ramin Aghababaei <ramin.aghababaei@epfl.ch>
 *
 * @date creation: Thu Oct 03 2013
 * @date last modification: Fri Jun 13 2014
 *
 * @brief  Specialization of the material class for isotropic finite deformation linear hardening plasticity
 *
 * @section LICENSE
 *
 * Copyright (©) 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "material_linear_isotropic_hardening.hh"
#include "solid_mechanics_model.hh"

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template<UInt dim>
MaterialLinearIsotropicHardening<dim>::MaterialLinearIsotropicHardening(SolidMechanicsModel & model,
                                                                        const ID & id) :
  Material(model, id), MaterialPlastic<dim>(model, id) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialLinearIsotropicHardening<spatial_dimension>::computeStress(ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  MaterialThermal<spatial_dimension>::computeStress(el_type, ghost_type);

  Array<Real>::iterator<> sigma_th_it =
    this->sigma_th(el_type, ghost_type).begin();

  Array<Real>::iterator<> previous_sigma_th_it =
    this->sigma_th.previous(el_type, ghost_type).begin();

  Array<Real>::matrix_iterator previous_gradu_it =
    this->gradu.previous(el_type, ghost_type).begin(spatial_dimension, spatial_dimension);

  Array<Real>::matrix_iterator previous_stress_it =
    this->stress.previous(el_type, ghost_type).begin(spatial_dimension, spatial_dimension);

  Array<Real>::matrix_iterator inelastic_strain_it =
    this->inelastic_strain(el_type, ghost_type).begin(spatial_dimension,spatial_dimension);

  Array<Real>::matrix_iterator previous_inelastic_strain_it =
    this->inelastic_strain.previous(el_type, ghost_type).begin(spatial_dimension,spatial_dimension);

  Array<Real>::iterator<> iso_hardening_it =
    this->iso_hardening(el_type, ghost_type).begin();

  Array<Real>::iterator<> previous_iso_hardening_it =
    this->iso_hardening.previous(el_type, ghost_type).begin();


  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type);
  Matrix<Real> & inelastic_strain_tensor = *inelastic_strain_it;
  Matrix<Real> & previous_inelastic_strain_tensor = *previous_inelastic_strain_it;
  Matrix<Real> & previous_grad_u = *previous_gradu_it;
  Matrix<Real> & previous_sigma = *previous_stress_it;

  computeStressOnQuad(grad_u,
                      previous_grad_u,
                      sigma,
                      previous_sigma,
                      inelastic_strain_tensor,
                      previous_inelastic_strain_tensor,
                      *iso_hardening_it,
                      *previous_iso_hardening_it,
                      *sigma_th_it,
                      *previous_sigma_th_it);
  ++sigma_th_it;
  ++inelastic_strain_it;
  ++iso_hardening_it;
  ++previous_sigma_th_it;
  ++previous_stress_it;
  ++previous_gradu_it;
  ++previous_inelastic_strain_it;
  ++previous_iso_hardening_it;

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
void MaterialLinearIsotropicHardening<spatial_dimension>::computeTangentModuli(__attribute__((unused)) const ElementType & el_type,
                                                                    Array<Real> & tangent_matrix,
                                                                    __attribute__((unused)) GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Array<Real>::const_matrix_iterator previous_gradu_it =
    this->gradu.previous(el_type, ghost_type).begin(spatial_dimension, spatial_dimension);

  Array<Real>::const_matrix_iterator previous_stress_it =
    this->stress.previous(el_type, ghost_type).begin(spatial_dimension, spatial_dimension);

  Array<Real>::const_scalar_iterator iso_hardening= this->iso_hardening(el_type, ghost_type).begin();

  MATERIAL_TANGENT_QUADRATURE_POINT_LOOP_BEGIN(tangent_matrix);

  computeTangentModuliOnQuad(tangent, grad_u, *previous_gradu_it, sigma_tensor, *previous_stress_it, *iso_hardening);

  ++previous_gradu_it;
  ++previous_stress_it;
  ++iso_hardening;

  MATERIAL_TANGENT_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

INSTANSIATE_MATERIAL(MaterialLinearIsotropicHardening);

__END_AKANTU__

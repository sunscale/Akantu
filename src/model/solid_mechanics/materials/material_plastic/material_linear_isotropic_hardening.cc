/**
 * @file   material_linear_isotropic_hardening.cc
 *
 * @author Ramin Aghababaei <ramin.aghababaei@epfl.ch>
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 * @author Benjamin Paccaud <benjamin.paccaud@epfl.ch>
 * @author Daniel Pino Muñoz <daniel.pinomunoz@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Mon Apr 07 2014
 * @date last modification: Sat Dec 02 2017
 *
 * @brief  Specialization of the material class for isotropic finite deformation
 * linear hardening plasticity
 *
 *
 * Copyright (©) 2014-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "material_linear_isotropic_hardening.hh"
#include "solid_mechanics_model.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
template <UInt dim>
MaterialLinearIsotropicHardening<dim>::MaterialLinearIsotropicHardening(
    SolidMechanicsModel & model, const ID & id)
    : MaterialPlastic<dim>(model, id) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
MaterialLinearIsotropicHardening<spatial_dimension>::
    MaterialLinearIsotropicHardening(SolidMechanicsModel & model, UInt dim,
                                     const Mesh & mesh, FEEngine & fe_engine,
                                     const ID & id)
    : MaterialPlastic<spatial_dimension>(model, dim, mesh, fe_engine, id) {}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialLinearIsotropicHardening<spatial_dimension>::computeStress(
    ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  // NOLINTNEXTLINE(bugprone-parent-virtual-call)
  MaterialThermal<spatial_dimension>::computeStress(el_type, ghost_type);

  // infinitesimal and finite deformation
  auto sigma_th_it = this->sigma_th(el_type, ghost_type).begin();

  auto previous_sigma_th_it =
      this->sigma_th.previous(el_type, ghost_type).begin();

  auto previous_gradu_it = this->gradu.previous(el_type, ghost_type)
                               .begin(spatial_dimension, spatial_dimension);

  auto previous_stress_it = this->stress.previous(el_type, ghost_type)
                                .begin(spatial_dimension, spatial_dimension);

  auto inelastic_strain_it = this->inelastic_strain(el_type, ghost_type)
                                 .begin(spatial_dimension, spatial_dimension);

  auto previous_inelastic_strain_it =
      this->inelastic_strain.previous(el_type, ghost_type)
          .begin(spatial_dimension, spatial_dimension);

  auto iso_hardening_it = this->iso_hardening(el_type, ghost_type).begin();

  auto previous_iso_hardening_it =
      this->iso_hardening.previous(el_type, ghost_type).begin();

  //
  // Finite Deformations
  //
  if (this->finite_deformation) {
    auto previous_piola_kirchhoff_2_it =
        this->piola_kirchhoff_2.previous(el_type, ghost_type)
            .begin(spatial_dimension, spatial_dimension);

    auto green_strain_it = this->green_strain(el_type, ghost_type)
                               .begin(spatial_dimension, spatial_dimension);

    MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type);

    auto & inelastic_strain_tensor = *inelastic_strain_it;
    auto & previous_inelastic_strain_tensor = *previous_inelastic_strain_it;
    auto & previous_grad_u = *previous_gradu_it;
    auto & previous_sigma = *previous_piola_kirchhoff_2_it;

    auto & green_strain = *green_strain_it;
    this->template gradUToE<spatial_dimension>(grad_u, green_strain);
    Matrix<Real> previous_green_strain(spatial_dimension, spatial_dimension);
    this->template gradUToE<spatial_dimension>(previous_grad_u,
                                                         previous_green_strain);
    Matrix<Real> F_tensor(spatial_dimension, spatial_dimension);
    this->template gradUToF<spatial_dimension>(grad_u, F_tensor);

    computeStressOnQuad(green_strain, previous_green_strain, sigma,
                        previous_sigma, inelastic_strain_tensor,
                        previous_inelastic_strain_tensor, *iso_hardening_it,
                        *previous_iso_hardening_it, *sigma_th_it,
                        *previous_sigma_th_it, F_tensor);

    ++sigma_th_it;
    ++inelastic_strain_it;
    ++iso_hardening_it;
    ++previous_sigma_th_it;
    //++previous_stress_it;
    ++previous_gradu_it;
    ++green_strain_it;
    ++previous_inelastic_strain_it;
    ++previous_iso_hardening_it;
    ++previous_piola_kirchhoff_2_it;

    MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  }
  // Infinitesimal deformations
  else {
    MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type);

    auto & inelastic_strain_tensor = *inelastic_strain_it;
    auto & previous_inelastic_strain_tensor = *previous_inelastic_strain_it;
    auto & previous_grad_u = *previous_gradu_it;
    auto & previous_sigma = *previous_stress_it;

    computeStressOnQuad(
        grad_u, previous_grad_u, sigma, previous_sigma, inelastic_strain_tensor,
        previous_inelastic_strain_tensor, *iso_hardening_it,
        *previous_iso_hardening_it, *sigma_th_it, *previous_sigma_th_it);
    ++sigma_th_it;
    ++inelastic_strain_it;
    ++iso_hardening_it;
    ++previous_sigma_th_it;
    ++previous_stress_it;
    ++previous_gradu_it;
    ++previous_inelastic_strain_it;
    ++previous_iso_hardening_it;

    MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialLinearIsotropicHardening<spatial_dimension>::computeTangentModuli(
    ElementType el_type, Array<Real> & tangent_matrix,
    GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  auto previous_gradu_it = this->gradu.previous(el_type, ghost_type)
                               .begin(spatial_dimension, spatial_dimension);

  auto previous_stress_it = this->stress.previous(el_type, ghost_type)
                                .begin(spatial_dimension, spatial_dimension);

  auto iso_hardening = this->iso_hardening(el_type, ghost_type).begin();

  MATERIAL_TANGENT_QUADRATURE_POINT_LOOP_BEGIN(tangent_matrix);

  computeTangentModuliOnQuad(tangent, grad_u, *previous_gradu_it, sigma,
                             *previous_stress_it, *iso_hardening);

  ++previous_gradu_it;
  ++previous_stress_it;
  ++iso_hardening;

  MATERIAL_TANGENT_QUADRATURE_POINT_LOOP_END;

  this->was_stiffness_assembled = true;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

INSTANTIATE_MATERIAL(plastic_linear_isotropic_hardening,
                     MaterialLinearIsotropicHardening);

} // namespace akantu

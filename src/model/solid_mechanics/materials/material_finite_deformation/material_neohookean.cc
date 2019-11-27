/**
 * @file   material_neohookean.cc
 *
 * @author Daniel Pino Muñoz <daniel.pinomunoz@epfl.ch>
 *
 * @date creation: Mon Apr 08 2013
 * @date last modification: Wed Nov 08 2017
 *
 * @brief  Specialization of the material class for finite deformation
 * neo-hookean material
 *
 * @section LICENSE
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
#include "material_neohookean.hh"
#include "solid_mechanics_model.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
MaterialNeohookean<spatial_dimension>::MaterialNeohookean(
    SolidMechanicsModel & model, const ID & id)
    : PlaneStressToolbox<spatial_dimension>(model, id) {
  AKANTU_DEBUG_IN();

  this->registerParam("E", E, Real(0.), _pat_parsable | _pat_modifiable,
                      "Young's modulus");
  this->registerParam("nu", nu, Real(0.5), _pat_parsable | _pat_modifiable,
                      "Poisson's ratio");
  this->registerParam("lambda", lambda, _pat_readable,
                      "First Lamé coefficient");
  this->registerParam("mu", mu, _pat_readable, "Second Lamé coefficient");
  this->registerParam("kapa", kpa, _pat_readable, "Bulk coefficient");

  this->finite_deformation = true;
  this->initialize_third_axis_deformation = true;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialNeohookean<spatial_dimension>::initMaterial() {
  AKANTU_DEBUG_IN();
  PlaneStressToolbox<spatial_dimension>::initMaterial();
  if (spatial_dimension == 1)
    nu = 0.;
  this->updateInternalParameters();
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <> void MaterialNeohookean<2>::initMaterial() {
  AKANTU_DEBUG_IN();
  PlaneStressToolbox<2>::initMaterial();
  this->updateInternalParameters();

  if (this->plane_stress)
    this->third_axis_deformation.setDefaultValue(1.);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialNeohookean<spatial_dimension>::updateInternalParameters() {
  lambda = nu * E / ((1 + nu) * (1 - 2 * nu));
  mu = E / (2 * (1 + nu));

  kpa = lambda + 2. / 3. * mu;
}

/* -------------------------------------------------------------------------- */
template <UInt dim>
void MaterialNeohookean<dim>::computeCauchyStressPlaneStress(
    ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  PlaneStressToolbox<dim>::computeCauchyStressPlaneStress(el_type, ghost_type);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <>
void MaterialNeohookean<2>::computeCauchyStressPlaneStress(
    ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  auto gradu_it = this->gradu(el_type, ghost_type).begin(2, 2);
  auto gradu_end = this->gradu(el_type, ghost_type).end(2, 2);
  auto piola_it = this->piola_kirchhoff_2(el_type, ghost_type).begin(2, 2);
  auto stress_it = this->stress(el_type, ghost_type).begin(2, 2);
  auto c33_it = this->third_axis_deformation(el_type, ghost_type).begin();

  for (; gradu_it != gradu_end; ++gradu_it, ++piola_it, ++stress_it, ++c33_it) {
    Matrix<Real> & grad_u = *gradu_it;
    Matrix<Real> & piola = *piola_it;
    Matrix<Real> & sigma = *stress_it;

    StoCauchy<2>(gradUToF<2>(grad_u), piola, sigma, *c33_it);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt dim>
void MaterialNeohookean<dim>::computeStress(ElementType el_type,
                                            GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type);
  computeStressOnQuad(grad_u, sigma);
  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <>
void MaterialNeohookean<2>::computeStress(ElementType el_type,
                                          GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  if (this->plane_stress) {
    PlaneStressToolbox<2>::computeStress(el_type, ghost_type);

    auto c33_it = this->third_axis_deformation(el_type, ghost_type).begin();

    MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type);
    computeStressOnQuad(grad_u, sigma, *c33_it);
    ++c33_it;
    MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;
  } else {

    MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type);
    computeStressOnQuad(grad_u, sigma);
    MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;
  }
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt dim>
void MaterialNeohookean<dim>::computeThirdAxisDeformation(
    ElementType /*el_type*/, GhostType /*ghost_type*/) {}

/* -------------------------------------------------------------------------- */
template <>
void MaterialNeohookean<2>::computeThirdAxisDeformation(ElementType el_type,
                                                        GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_ASSERT(this->plane_stress, "The third component of the strain "
                                          "can only be computed for 2D "
                                          "problems in Plane Stress!!");

  Array<Real>::scalar_iterator c33_it =
      this->third_axis_deformation(el_type, ghost_type).begin();

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type);
  computeThirdAxisDeformationOnQuad(grad_u, *c33_it);
  ++c33_it;
  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialNeohookean<spatial_dimension>::computePotentialEnergy(
    ElementType el_type) {
  AKANTU_DEBUG_IN();

  Material::computePotentialEnergy(el_type);

  Array<Real>::scalar_iterator epot = this->potential_energy(el_type).begin();

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, _not_ghost);

  computePotentialEnergyOnQuad(grad_u, *epot);
  ++epot;

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialNeohookean<spatial_dimension>::computeTangentModuli(
    __attribute__((unused)) const ElementType & el_type,
    Array<Real> & tangent_matrix,
    __attribute__((unused)) GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  MATERIAL_TANGENT_QUADRATURE_POINT_LOOP_BEGIN(tangent_matrix);
  computeTangentModuliOnQuad(tangent, grad_u);
  MATERIAL_TANGENT_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <>
void MaterialNeohookean<2>::computeTangentModuli(__attribute__((unused))
                                                 const ElementType & el_type,
                                                 Array<Real> & tangent_matrix,
                                                 __attribute__((unused))
                                                 GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  if (this->plane_stress) {
    PlaneStressToolbox<2>::computeStress(el_type, ghost_type);

    Array<Real>::const_scalar_iterator c33_it =
        this->third_axis_deformation(el_type, ghost_type).begin();

    MATERIAL_TANGENT_QUADRATURE_POINT_LOOP_BEGIN(tangent_matrix);
    computeTangentModuliOnQuad(tangent, grad_u, *c33_it);
    ++c33_it;
    MATERIAL_TANGENT_QUADRATURE_POINT_LOOP_END;

  } else {

    MATERIAL_TANGENT_QUADRATURE_POINT_LOOP_BEGIN(tangent_matrix);
    computeTangentModuliOnQuad(tangent, grad_u);
    MATERIAL_TANGENT_QUADRATURE_POINT_LOOP_END;
  }
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
Real MaterialNeohookean<spatial_dimension>::getPushWaveSpeed(
    __attribute__((unused)) const Element & element) const {
  return sqrt((this->lambda + 2 * this->mu) / this->rho);
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
Real MaterialNeohookean<spatial_dimension>::getShearWaveSpeed(
    __attribute__((unused)) const Element & element) const {
  return sqrt(this->mu / this->rho);
}

/* -------------------------------------------------------------------------- */

INSTANTIATE_MATERIAL(neohookean, MaterialNeohookean);

} // namespace akantu

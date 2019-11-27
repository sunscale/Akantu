/**
 * @file   material_elastic.cc
 *
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 * @author Daniel Pino Muñoz <daniel.pinomunoz@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Mon Jan 29 2018
 *
 * @brief  Specialization of the material class for the elastic material
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
#include "material_elastic.hh"
#include "solid_mechanics_model.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
template <UInt dim>
MaterialElastic<dim>::MaterialElastic(SolidMechanicsModel & model,
                                      const ID & id)
    : Parent(model, id), was_stiffness_assembled(false) {
  AKANTU_DEBUG_IN();
  this->initialize();
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt dim>
MaterialElastic<dim>::MaterialElastic(SolidMechanicsModel & model,
                                      __attribute__((unused)) UInt a_dim,
                                      const Mesh & mesh, FEEngine & fe_engine,
                                      const ID & id)
    : Parent(model, dim, mesh, fe_engine, id), was_stiffness_assembled(false) {
  AKANTU_DEBUG_IN();
  this->initialize();
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt dim> void MaterialElastic<dim>::initialize() {
  this->registerParam("lambda", lambda, _pat_readable,
                      "First Lamé coefficient");
  this->registerParam("mu", mu, _pat_readable, "Second Lamé coefficient");
  this->registerParam("kapa", kpa, _pat_readable, "Bulk coefficient");
}

/* -------------------------------------------------------------------------- */
template <UInt dim> void MaterialElastic<dim>::initMaterial() {
  AKANTU_DEBUG_IN();
  Parent::initMaterial();

  if (dim == 1)
    this->nu = 0.;

  this->updateInternalParameters();
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt dim> void MaterialElastic<dim>::updateInternalParameters() {
  MaterialThermal<dim>::updateInternalParameters();

  this->lambda = this->nu * this->E / ((1 + this->nu) * (1 - 2 * this->nu));
  this->mu = this->E / (2 * (1 + this->nu));

  this->kpa = this->lambda + 2. / 3. * this->mu;

  this->was_stiffness_assembled = false;
}

/* -------------------------------------------------------------------------- */
template <> void MaterialElastic<2>::updateInternalParameters() {
  MaterialThermal<2>::updateInternalParameters();

  this->lambda = this->nu * this->E / ((1 + this->nu) * (1 - 2 * this->nu));
  this->mu = this->E / (2 * (1 + this->nu));

  if (this->plane_stress)
    this->lambda = this->nu * this->E / ((1 + this->nu) * (1 - this->nu));

  this->kpa = this->lambda + 2. / 3. * this->mu;

  this->was_stiffness_assembled = false;
}

/* -------------------------------------------------------------------------- */
template <UInt dim>
void MaterialElastic<dim>::computeStress(ElementType el_type,
                                         GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Parent::computeStress(el_type, ghost_type);

  Array<Real>::const_scalar_iterator sigma_th_it =
      this->sigma_th(el_type, ghost_type).begin();

  if (!this->finite_deformation) {
    MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type);
    const Real & sigma_th = *sigma_th_it;
    this->computeStressOnQuad(grad_u, sigma, sigma_th);
    ++sigma_th_it;
    MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;
  } else {
    /// finite gradus
    Matrix<Real> E(dim, dim);

    MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type);

    /// compute E
    this->template gradUToE<dim>(grad_u, E);

    const Real & sigma_th = *sigma_th_it;

    /// compute second Piola-Kirchhoff stress tensor
    this->computeStressOnQuad(E, sigma, sigma_th);

    ++sigma_th_it;
    MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt dim>
void MaterialElastic<dim>::computeTangentModuli(const ElementType & el_type,
                                                Array<Real> & tangent_matrix,
                                                GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  MATERIAL_TANGENT_QUADRATURE_POINT_LOOP_BEGIN(tangent_matrix);
  this->computeTangentModuliOnQuad(tangent);
  MATERIAL_TANGENT_QUADRATURE_POINT_LOOP_END;

  this->was_stiffness_assembled = true;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt dim>
Real MaterialElastic<dim>::getPushWaveSpeed(const Element &) const {
  return sqrt((lambda + 2 * mu) / this->rho);
}

/* -------------------------------------------------------------------------- */
template <UInt dim>
Real MaterialElastic<dim>::getShearWaveSpeed(const Element &) const {
  return sqrt(mu / this->rho);
}

/* -------------------------------------------------------------------------- */
template <UInt dim>
void MaterialElastic<dim>::computePotentialEnergy(ElementType el_type) {
  AKANTU_DEBUG_IN();

  // MaterialThermal<dim>::computePotentialEnergy(ElementType)
  // needs to be implemented
  // MaterialThermal<dim>::computePotentialEnergy(el_type);

  auto epot = this->potential_energy(el_type, _not_ghost).begin();

  if (!this->finite_deformation) {
    MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, _not_ghost);

    this->computePotentialEnergyOnQuad(grad_u, sigma, *epot);
    ++epot;

    MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;
  } else {

    MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, _not_ghost);
    auto E = this->template gradUToE<dim>(grad_u);

    this->computePotentialEnergyOnQuad(E, sigma, *epot);
    ++epot;

    MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt dim>
void MaterialElastic<dim>::computePotentialEnergyByElement(
    ElementType type, UInt index, Vector<Real> & epot_on_quad_points) {
  auto gradu_it = this->gradu(type).begin(dim, dim);
  auto gradu_end = this->gradu(type).begin(dim, dim);
  auto stress_it = this->stress(type).begin(dim, dim);

  if (this->finite_deformation)
    stress_it = this->piola_kirchhoff_2(type).begin(dim, dim);

  UInt nb_quadrature_points = this->fem.getNbIntegrationPoints(type);

  gradu_it += index * nb_quadrature_points;
  gradu_end += (index + 1) * nb_quadrature_points;
  stress_it += index * nb_quadrature_points;

  Real * epot_quad = epot_on_quad_points.storage();

  Matrix<Real> grad_u(dim, dim);

  if (this->finite_deformation) {
    for (; gradu_it != gradu_end; ++gradu_it, ++stress_it, ++epot_quad) {
      auto E = this->template gradUToE<dim>(*gradu_it);
      this->computePotentialEnergyOnQuad(E, *stress_it, *epot_quad);
    }
  } else {
    for (; gradu_it != gradu_end; ++gradu_it, ++stress_it, ++epot_quad) {
      this->computePotentialEnergyOnQuad(*gradu_it, *stress_it, *epot_quad);
    }
  }
}

/* -------------------------------------------------------------------------- */
template <>
Real MaterialElastic<1>::getPushWaveSpeed(const Element & /*element*/) const {
  return std::sqrt(this->E / this->rho);
}

template <>
Real MaterialElastic<1>::getShearWaveSpeed(const Element & /*element*/) const {
  AKANTU_EXCEPTION("There is no shear wave speed in 1D");
}
/* -------------------------------------------------------------------------- */

INSTANTIATE_MATERIAL(elastic, MaterialElastic);

} // namespace akantu

/**
 * @file   material_viscoelastic_maxwell.hh
 *
 * @author Emil Gallyamov <emil.gallyamov@epfl.ch>
 *
 * @date creation: Tue May 08 2018
 * @date last modification: Tue May 08 2018
 *
 * @brief  Material Visco-elastic, based on Maxwell chain,
 * see
 * [] R. de Borst and A.H. van den Boogaard "Finite-element modeling of
 * deformation and cracking in early-age concrete", J.Eng.Mech., 1994
 * as well as
 * [] Manual of DIANA FEA Theory manual v.10.2 Section 37.6
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
#include "material_viscoelastic_maxwell.hh"
#include "solid_mechanics_model.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
MaterialViscoelasticMaxwell<spatial_dimension>::MaterialViscoelasticMaxwell(
    SolidMechanicsModel & model, const ID & id)
    : MaterialElastic<spatial_dimension>(model, id),
      C(voigt_h::size, voigt_h::size), sigma_v("sigma_v", *this),
    dissipated_energy("dissipated_energy", *this) {

  AKANTU_DEBUG_IN();

  this->registerParam("Einf", Einf, Real(1.), _pat_parsable | _pat_modifiable,
                      "Stiffness of the elastic element");
  this->registerParam("previous_dt", previous_dt, Real(0.), _pat_readable,
                      "Time step of previous solveStep");
  this->registerParam("Eta", Eta, _pat_parsable | _pat_modifiable,
                      "Viscosity of a Maxwell element");
  this->registerParam("Ev", Ev, _pat_parsable | _pat_modifiable,
                      "Stiffness of a Maxwell element");
  this->use_previous_stress = true;
  this->use_previous_gradu = true;
  this->use_previous_stress_thermal = true;

  this->dissipated_energy.initialize(1);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialViscoelasticMaxwell<spatial_dimension>::initMaterial() {
  AKANTU_DEBUG_IN();

  // this->E = Einf + Ev1;
  this->E = std::min(this->Einf, this->Ev(0));
  MaterialElastic<spatial_dimension>::initMaterial();
  AKANTU_DEBUG_ASSERT(this->Eta.size() == this->Ev.size(), "Eta and Ev have different dimensions! Please correct.");

  UInt stress_size = spatial_dimension * spatial_dimension;
  this->sigma_v.initialize(stress_size * this->Ev.size());

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialViscoelasticMaxwell<
    spatial_dimension>::updateInternalParameters() {
  MaterialElastic<spatial_dimension>::updateInternalParameters();

  Real pre_mult = 1 / (1 + this->nu) / (1 - 2 * this->nu);
  UInt n = voigt_h::size;
  Real Miiii = pre_mult * (1 - this->nu);
  Real Miijj = pre_mult * this->nu;
  Real Mijij = pre_mult * 0.5 * (1 - 2 * this->nu);

  if (spatial_dimension == 1) {
    C(0, 0) = 1;
  } else {
    C(0, 0) = Miiii;
  }
  if (spatial_dimension >= 2) {
    C(1, 1) = Miiii;
    C(0, 1) = Miijj;
    C(1, 0) = Miijj;
    C(n - 1, n - 1) = Mijij;
  }

  if (spatial_dimension == 3) {
    C(2, 2) = Miiii;
    C(0, 2) = Miijj;
    C(1, 2) = Miijj;
    C(2, 0) = Miijj;
    C(2, 1) = Miijj;
    C(3, 3) = Mijij;
    C(4, 4) = Mijij;
  }
}

/* -------------------------------------------------------------------------- */
template <> void MaterialViscoelasticMaxwell<2>::updateInternalParameters() {
  MaterialElastic<2>::updateInternalParameters();

  Real pre_mult = 1 / (1 + this->nu) / (1 - 2 * this->nu);
  UInt n = voigt_h::size;
  Real Miiii = pre_mult * (1 - this->nu);
  Real Miijj = pre_mult * this->nu;
  Real Mijij = pre_mult * 0.5 * (1 - 2 * this->nu);

  C(0, 0) = Miiii;
  C(1, 1) = Miiii;
  C(0, 1) = Miijj;
  C(1, 0) = Miijj;
  C(n - 1, n - 1) = Mijij;
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialViscoelasticMaxwell<spatial_dimension>::setToSteadyState(
    ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  /* Array<Real> & stress_dev_vect = stress_dev(el_type, ghost_type);
  Array<Real> & history_int_vect = history_integral(el_type, ghost_type);

  Array<Real>::matrix_iterator stress_d =
      stress_dev_vect.begin(spatial_dimension, spatial_dimension);
  Array<Real>::matrix_iterator history_int =
      history_int_vect.begin(spatial_dimension, spatial_dimension);

  /// Loop on all quadrature points
  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type);

  Matrix<Real> & dev_s = *stress_d;
  Matrix<Real> & h = *history_int;

  /// Compute the first invariant of strain
  Real Theta = grad_u.trace();

  for (UInt i = 0; i < spatial_dimension; ++i)
    for (UInt j = 0; j < spatial_dimension; ++j) {
      dev_s(i, j) = 2 * this->mu * (.5 * (grad_u(i, j) + grad_u(j, i)) -
                                    1. / 3. * Theta * (i == j));
      h(i, j) = 0.;
    }

  /// Save the deviator of stress
  ++stress_d;
  ++history_int;

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;
  */
  AKANTU_DEBUG_OUT();
}
/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialViscoelasticMaxwell<spatial_dimension>::computeStress(
    ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  MaterialThermal<spatial_dimension>::computeStress(el_type, ghost_type);

  auto sigma_th_it = this->sigma_th(el_type, ghost_type).begin();

  auto previous_sigma_th_it =
      this->sigma_th.previous(el_type, ghost_type).begin();

  auto previous_gradu_it = this->gradu.previous(el_type, ghost_type)
                               .begin(spatial_dimension, spatial_dimension);

  auto previous_stress_it = this->stress.previous(el_type, ghost_type)
                                .begin(spatial_dimension, spatial_dimension);

  auto sigma_v_it = this->sigma_v(el_type, ghost_type)
      .begin(spatial_dimension, spatial_dimension, this->Eta.size());

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type);

  auto & previous_grad_u = *previous_gradu_it;
  auto & previous_sigma = *previous_stress_it;

  computeStressOnQuad(grad_u, previous_grad_u, sigma, previous_sigma,
                      *sigma_v_it, *sigma_th_it, *previous_sigma_th_it);
  ++sigma_th_it;
  ++previous_sigma_th_it;
  ++previous_stress_it;
  ++previous_gradu_it;
  ++sigma_v_it;

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}


/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialViscoelasticMaxwell<spatial_dimension>::computeStressOnQuad(
    const Matrix<Real> & grad_u, const Matrix<Real> & previous_grad_u,
    Matrix<Real> & sigma, const Matrix<Real> & previous_sigma,
    Tensor3<Real> & sigma_v, const Real & sigma_th,
    const Real & previous_sigma_th) {

  Matrix<Real> grad_delta_u(grad_u);
  grad_delta_u -= previous_grad_u;
  Real delta_sigma_th = sigma_th - previous_sigma_th;

  // Wikipedia convention:
  // 2*eps_ij (i!=j) = voigt_eps_I
  // http://en.wikipedia.org/wiki/Voigt_notation
  Vector<Real> voigt_current_strain(voigt_h::size);
  Vector<Real> voigt_previous_strain(voigt_h::size);
  Vector<Real> voigt_stress(voigt_h::size);
  Vector<Real> voigt_sigma_v(voigt_h::size);

  for (UInt I = 0; I < voigt_h::size; ++I) {
    Real voigt_factor = voigt_h::factors[I];
    UInt i = voigt_h::vec[I][0];
    UInt j = voigt_h::vec[I][1];

    voigt_current_strain(I) =
        voigt_factor * (grad_u(i, j) + grad_u(j, i)) / 2.;
    voigt_previous_strain(I) =
        voigt_factor * (previous_grad_u(i, j) + previous_grad_u(j, i)) / 2.;
  }

  voigt_stress =
      this->Einf * this->C * voigt_current_strain;
  Real dt = this->model.getTimeStep();

  for (UInt k = 0; k < Eta.size(); ++k) {
    Real lambda = this->Eta(k) / this->Ev(k);
    Real exp_dt_lambda = exp(-dt / lambda);
    Real E_additional = (1 - exp_dt_lambda) * this->Ev(k) * lambda / dt;

    for (UInt I = 0; I < voigt_h::size; ++I) {
      UInt i = voigt_h::vec[I][0];
      UInt j = voigt_h::vec[I][1];

      voigt_sigma_v(I) = sigma_v(i, j, k);
    }

    voigt_stress += E_additional * this->C *
        (voigt_current_strain - voigt_previous_strain) +
        exp_dt_lambda * voigt_sigma_v;
  }

  for (UInt I = 0; I < voigt_h::size; ++I) {
    UInt i = voigt_h::vec[I][0];
    UInt j = voigt_h::vec[I][1];

    sigma(i, j) = sigma(j, i) =
        voigt_stress(I) + (i == j) * sigma_th;
  }
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialViscoelasticMaxwell<spatial_dimension>::afterSolveStep() {

  Material::afterSolveStep();

  for (auto & el_type : this->element_filter.elementTypes(
           _all_dimensions, _not_ghost, _ek_not_defined)) {
    auto previous_gradu_it = this->gradu.previous(el_type, _not_ghost)
                                 .begin(spatial_dimension, spatial_dimension);

    auto sigma_v_it = this->sigma_v(el_type, _not_ghost)
        .begin(spatial_dimension, spatial_dimension, this->Eta.size());

    MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, _not_ghost);

    auto & previous_grad_u = *previous_gradu_it;

    updateIntVarOnQuad(grad_u, previous_grad_u, *sigma_v_it);

    ++previous_gradu_it;
    ++sigma_v_it;

    MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;
    this->updateDissipatedEnergy(el_type, _not_ghost);
  }
}
/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialViscoelasticMaxwell<spatial_dimension>::updateIntVarOnQuad(
    const Matrix<Real> & grad_u, const Matrix<Real> & previous_grad_u,
    Tensor3<Real> & sigma_v) {

  Matrix<Real> grad_delta_u(grad_u);
  grad_delta_u -= previous_grad_u;

  Real dt = this->model.getTimeStep();
  Vector<Real> voigt_delta_strain(voigt_h::size);
  for (UInt I = 0; I < voigt_h::size; ++I) {
      Real voigt_factor = voigt_h::factors[I];
      UInt i = voigt_h::vec[I][0];
      UInt j = voigt_h::vec[I][1];

      voigt_delta_strain(I) =
          voigt_factor * (grad_delta_u(i, j) + grad_delta_u(j, i)) / 2.;
  }


  for (UInt k = 0; k < Eta.size(); ++k) {
    Real lambda = this->Eta(k) / this->Ev(k);
    Real exp_dt_lambda = exp(-dt / lambda);
    Real E_ef_v = (1 - exp_dt_lambda) * this->Ev(k) * lambda / dt;

    Vector<Real> voigt_sigma_v(voigt_h::size);

    for (UInt I = 0; I < voigt_h::size; ++I) {
      Real voigt_factor = voigt_h::factors[I];
      UInt i = voigt_h::vec[I][0];
      UInt j = voigt_h::vec[I][1];

      voigt_sigma_v(I) = sigma_v(i, j, k);
    }

    voigt_sigma_v =
        exp_dt_lambda * voigt_sigma_v + E_ef_v * this->C * voigt_delta_strain;

    for (UInt I = 0; I < voigt_h::size; ++I) {
      UInt i = voigt_h::vec[I][0];
      UInt j = voigt_h::vec[I][1];

      sigma_v(i, j, k) = sigma_v(j, i, k) = voigt_sigma_v(I);
    }
  }
}
/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialViscoelasticMaxwell<spatial_dimension>::computeTangentModuli(
    const ElementType & el_type, Array<Real> & tangent_matrix,
    GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Real dt = this->model.getTimeStep();
  Real E_ef = this->Einf;

  for (UInt k = 0; k < Eta.size(); ++k) {
    Real lambda = this->Eta(k) / this->Ev(k);
    Real exp_dt_lambda = exp(-dt / lambda);
    E_ef += (1 - exp_dt_lambda) * this->Ev(k) * lambda / dt;
  }

  this->previous_dt = dt;

  MATERIAL_TANGENT_QUADRATURE_POINT_LOOP_BEGIN(tangent_matrix);
  this->computeTangentModuliOnQuad(tangent);
  MATERIAL_TANGENT_QUADRATURE_POINT_LOOP_END;

  tangent_matrix *= E_ef;

  this->was_stiffness_assembled = true;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialViscoelasticMaxwell<spatial_dimension>::computeTangentModuliOnQuad(
    Matrix<Real> & tangent) {

  tangent.copy(C);
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialViscoelasticMaxwell<spatial_dimension>::savePreviousState() {

  for (auto & el_type : this->element_filter.elementTypes(
           _all_dimensions, _not_ghost, _ek_not_defined)) {

    auto sigma_th_it = this->sigma_th(el_type, _not_ghost).begin();

    auto previous_sigma_th_it =
        this->sigma_th.previous(el_type, _not_ghost).begin();

    auto previous_gradu_it = this->gradu.previous(el_type, _not_ghost)
                                 .begin(spatial_dimension, spatial_dimension);

    auto previous_sigma_it = this->stress.previous(el_type, _not_ghost)
                                 .begin(spatial_dimension, spatial_dimension);

    auto sigma_v_it = this->sigma_v(el_type, _not_ghost)
        .begin(spatial_dimension, spatial_dimension, this->Eta.size());

    MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, _not_ghost);
    auto & previous_grad_u = *previous_gradu_it;
    auto & previous_sigma = *previous_sigma_it;

    previous_grad_u.copy(*gradu_it);
    previous_sigma.copy(*stress_it);
    *previous_sigma_th_it = *sigma_th_it;

    ++previous_gradu_it, ++previous_sigma_it, ++previous_sigma_th_it,
        ++sigma_v_it, ++sigma_th_it;
    MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;
  }
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialViscoelasticMaxwell<spatial_dimension>::updateIntVariables() {

  for (auto & el_type : this->element_filter.elementTypes(
           _all_dimensions, _not_ghost, _ek_not_defined)) {

    auto previous_gradu_it = this->gradu.previous(el_type, _not_ghost)
                                 .begin(spatial_dimension, spatial_dimension);
    auto previous_sigma_it = this->stress.previous(el_type, _not_ghost)
                                 .begin(spatial_dimension, spatial_dimension);

    auto sigma_v_it = this->sigma_v(el_type, _not_ghost)
        .begin(spatial_dimension, spatial_dimension, this->Eta.size());

    MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, _not_ghost);
    auto & previous_grad_u = *previous_gradu_it;

    updateIntVarOnQuad(grad_u, previous_grad_u, *sigma_v_it);

    ++previous_gradu_it, ++sigma_v_it;
    MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;
  }
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void MaterialViscoelasticMaxwell<spatial_dimension>::updateDissipatedEnergy(
    ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  // if(ghost_type == _ghost) return 0.;

  Real * dis_energy = dissipated_energy(el_type, ghost_type).storage();

  /*  Array<Real> & sigma_vec = stress_dev(el_type, ghost_type);
    Array<Real> & history_int_vect = history_integral(el_type, ghost_type);

    Array<Real>::matrix_iterator stress_d =
        stress_dev_vect.begin(spatial_dimension, spatial_dimension);
    Array<Real>::matrix_iterator history_int =
        history_int_vect.begin(spatial_dimension, spatial_dimension);

    Matrix<Real> q(spatial_dimension, spatial_dimension);
    Matrix<Real> q_rate(spatial_dimension, spatial_dimension);
    Matrix<Real> epsilon_d(spatial_dimension, spatial_dimension);
    Matrix<Real> epsilon_v(spatial_dimension, spatial_dimension);

    Real dt = this->model.getTimeStep();

    Real gamma_v = Ev / this->E;
    Real alpha = 1. / (2. * this->mu * gamma_v);

    /// Loop on all quadrature points
    MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type);

    Matrix<Real> & dev_s = *stress_d;
    Matrix<Real> & h = *history_int;

    /// Compute the first invariant of strain
    this->template gradUToEpsilon<spatial_dimension>(grad_u, epsilon_d);

    Real Theta = epsilon_d.trace();
    epsilon_v.eye(Theta / Real(3.));
    epsilon_d -= epsilon_v;

    q.copy(dev_s);
    q -= h;
    q *= gamma_v;

    q_rate.copy(dev_s);
    q_rate *= gamma_v;
    q_rate -= q;
    q_rate /= tau;

    for (UInt i = 0; i < spatial_dimension; ++i)
      for (UInt j = 0; j < spatial_dimension; ++j)
        *dis_energy += (epsilon_d(i, j) - alpha * q(i, j)) * q_rate(i, j) * dt;

    /// Save the deviator of stress
    ++stress_d;
    ++history_int;
    ++dis_energy;

    MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;
  */
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
Real MaterialViscoelasticMaxwell<spatial_dimension>::getDissipatedEnergy()
    const {
  AKANTU_DEBUG_IN();

  Real de = 0.;

  /// integrate the dissipated energy for each type of elements
  for (auto & type :
       this->element_filter.elementTypes(spatial_dimension, _not_ghost)) {
    de +=
        this->fem.integrate(dissipated_energy(type, _not_ghost), type,
                            _not_ghost, this->element_filter(type, _not_ghost));
  }

  AKANTU_DEBUG_OUT();
  return de;
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
Real MaterialViscoelasticMaxwell<spatial_dimension>::getDissipatedEnergy(
    ElementType type, UInt index) const {
  AKANTU_DEBUG_IN();

  UInt nb_quadrature_points = this->fem.getNbIntegrationPoints(type);
  auto it =
      this->dissipated_energy(type, _not_ghost).begin(nb_quadrature_points);
  UInt gindex = (this->element_filter(type, _not_ghost))(index);

  AKANTU_DEBUG_OUT();
  return this->fem.integrate(it[index], type, gindex);
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
Real MaterialViscoelasticMaxwell<spatial_dimension>::getEnergy(
    const std::string & type) {
  if (type == "dissipated")
    return getDissipatedEnergy();
  else if (type == "viscoelastic_maxwell")
    return getDissipatedEnergy();
  else
    return MaterialElastic<spatial_dimension>::getEnergy(type);
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
Real MaterialViscoelasticMaxwell<spatial_dimension>::getEnergy(
    const std::string & energy_id, ElementType type, UInt index) {
  if (energy_id == "dissipated")
    return getDissipatedEnergy(type, index);
  else if (energy_id == "viscoelastic_maxwell")
    return getDissipatedEnergy(type, index);
  else
    return MaterialElastic<spatial_dimension>::getEnergy(energy_id, type,
                                                         index);
}

/* -------------------------------------------------------------------------- */

INSTANTIATE_MATERIAL(viscoelastic_maxwell, MaterialViscoelasticMaxwell);

} // namespace akantu

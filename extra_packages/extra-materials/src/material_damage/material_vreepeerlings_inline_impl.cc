/**
 * @file   material_vreepeerlings_inline_impl.cc
 *
 * @author Cyprien Wolff <cyprien.wolff@epfl.ch>
 *
 *
 * @brief  Specialization of the material class for the VreePeerlings material
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */

template <UInt spatial_dimension, template <UInt> class MatParent>
inline void
MaterialVreePeerlings<spatial_dimension, MatParent>::computeStressOnQuad(
    Matrix<Real> & grad_u, Matrix<Real> & sigma, Real & dam, Real & Equistrain,
    Real & Equistrain_rate, Real & Kapaq, Real dt,
    Matrix<Real> & strain_rate_vrplgs, Real & FullDam_ValStrain,
    Real & FullDam_ValStrain_rate, Real & Nb_damage) {

  Real I1 = 0.;     /// trace initialization of the strain tensor
  Real J2 = 0.;     /// J2 [J2=1/6(3 strain:strain - I1*I1)] initialization
  Real I1rate = 0.; /// trace initialization of the strain rate tensor
  Real J2rate =
      0.; /// J2 [J2=1/3(3 strain:strainrate - I1*I1rate)] initialization

  if (spatial_dimension == 1) {

    I1 = grad_u.trace();
    J2 = .5 * grad_u(0, 0) * grad_u(0, 0) - I1 * I1 / 6.;

    I1rate = strain_rate_vrplgs.trace();
    J2rate = grad_u(0, 0) * strain_rate_vrplgs(0, 0) - I1 * I1rate / 6.;
  } else {
    // if(this->plane_stress) {

    //   Real tmp = this->nu/(this->nu - 1);
    //   tmp *= tmp;

    //   I1 = (grad_u(0,0) + grad_u(1,1))*(1 - 2*this->nu)/(1 - this->nu);
    //   J2 = .5*(grad_u(0,0)*grad_u(0,0) +
    //            grad_u(1,1)*grad_u(1,1) +
    //            tmp*(grad_u(0,0) + grad_u(1,1))*(grad_u(0,0) + grad_u(1,1)) +
    //            .5*(grad_u(0,1) + grad_u(1,0))*(grad_u(0,1) + grad_u(1,0))) -
    //     I1*I1/6.;

    //   I1rate = (strain_rate_vrplgs(0,0) + strain_rate_vrplgs(1,1))*(1 -
    //   2*this->nu)/(1 - this->nu);
    //   J2rate = (grad_u(0,0)*strain_rate_vrplgs(0,0) +
    //             grad_u(1,1)*strain_rate_vrplgs(1,1) +
    //             tmp*(grad_u(0,0) + grad_u(1,1))*(strain_rate_vrplgs(0,0) +
    //             strain_rate_vrplgs(1,1)) +
    //             (grad_u(0,1)*strain_rate_vrplgs(1,0)) +
    //             (grad_u(0,1)*strain_rate_vrplgs(1,0))) -
    //     I1*I1rate/3.;

    // }
    // else {

    I1 = grad_u.trace();
    for (UInt i = 0; i < spatial_dimension; ++i)
      for (UInt j = i; j < spatial_dimension; ++j)
        J2 += 0.5 * ((i == j) * grad_u(i, i) * grad_u(i, i) +
                     0.5 * (1 - (i == j)) *
                         ((grad_u(i, j) + grad_u(j, i)) *
                          (grad_u(i, j) + grad_u(j, i))));

    J2 -= I1 * I1 / 6.;

    I1rate = strain_rate_vrplgs.trace();
    bool is3D = spatial_dimension == 3;
    J2rate = (grad_u(0, 0) * strain_rate_vrplgs(0, 0) +
              grad_u(1, 1) * strain_rate_vrplgs(1, 1) +
              is3D * grad_u(2, 2) * strain_rate_vrplgs(2, 2) +
              (grad_u(0, 1) * strain_rate_vrplgs(1, 0)) +
              (grad_u(0, 1) * strain_rate_vrplgs(1, 0)) +
              is3D * (grad_u(1, 2) * strain_rate_vrplgs(2, 1)) +
              is3D * (grad_u(1, 2) * strain_rate_vrplgs(2, 1)) +
              is3D * (grad_u(2, 0) * strain_rate_vrplgs(0, 2)) +
              is3D * (grad_u(2, 0) * strain_rate_vrplgs(0, 2))) -
             I1 * I1rate / 3.;

    //    }
  }

  Real tmp_1 = (Kct - 1) / (1 - 2 * this->nu);
  Real tmp_2 = (12 * Kct) / ((1 + this->nu) * (1 + this->nu));

  Equistrain = tmp_1 * I1 / (2 * Kct) +
               1. / (2 * Kct) * std::sqrt(tmp_1 * tmp_1 * I1 * I1 + tmp_2 * J2);

  if (I1 * I1rate > 0 || J2rate > 0) {
    Equistrain_rate = tmp_1 * std::abs(I1rate) / (2 * Kct) +
                      1. / (4 * Kct) *
                          (2 * tmp_1 * tmp_1 * I1 * I1rate + tmp_2 * J2rate) /
                          std::sqrt(tmp_1 * tmp_1 * I1 * I1 + tmp_2 * J2);
  } else {
    AKANTU_ERROR("This instruction here was commented but has to be checked");
    // Equistrain_rate = Equistrain_rate;
  }

  if (!this->is_non_local) {
    computeDamageAndStressOnQuad(sigma, dam, Equistrain, Equistrain_rate, Kapaq,
                                 dt, FullDam_ValStrain, FullDam_ValStrain_rate,
                                 Nb_damage);
  }
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension, template <UInt> class MatParent>
inline void MaterialVreePeerlings<spatial_dimension, MatParent>::
    computeDamageAndStressOnQuad(Matrix<Real> & sigma, Real & dam,
                                 Real & Equistrain, Real & Equistrain_rate,
                                 Real & Kapaq, Real dt,
                                 Real & FullDam_ValStrain,
                                 Real & FullDam_ValStrain_rate,
                                 Real & Nb_damage) {

  //---------------------------------------------------------------------------------------
  // rate-dependence  model
  //
  //---------------------------------------------------------------------------------------
  Real Kapao = Crate *
                   (1. + (std::pow(std::abs(Equistrain_rate) * Arate, Brate))) *
                   (1. - Kapao_randomness) +
               Kapao_randomness * Kapaq;
  Real Kapac_99 = Kapac * .99;
  Real Kapaodyn = 0;
  if (Kapao > Kapaoi) {
    if (Kapao < Kapac_99) {
      Kapaodyn = Kapao;
    } else {
      Kapaodyn = Kapac_99;
    }
  } else {
    Kapaodyn = Kapaoi;
  }

  Real Fd1 = Equistrain - Kapaoi;
  if (Fd1 > 0) {
    Real dam1 = 1. - Kapaodyn / Equistrain *
                         ((Kapac - Equistrain) / (Kapac - Kapaodyn));
    if (dam1 > dam) {
      dam = std::min(dam1, 1.);
      Nb_damage = Nb_damage + 1.;

      if (dam >= 1.) {
        FullDam_ValStrain = Equistrain;
        FullDam_ValStrain_rate = Equistrain_rate;
      }
    }
  }

  //---------------------------------------------------------------------------------------
  // delayed damage (see Marions thesis page 68)
  //---------------------------------------------------------------------------------------
  //   Real viscosity = 10.;
  //   Real damRateInfini = 10000000.;
  //   Real gequi = 1. - Kapaq/Equistrain * ((Kapac-Equistrain)/(Kapac -
  //   Kapaq));
  //   if (gequi - dam > 0){
  //
  //     Real damRate = damRateInfini * (1. - std::exp(-1*viscosity * (gequi -
  //     dam)));
  //     if (damrate > 0){
  //
  //	if (dam < 1.){
  //      dam = dam + damRate*dt;
  //	} else {
  //      dam = 1.;
  //	}
  //     }
  //   }
  //
  //---------------------------------------------------------------------------------------
  sigma *= 1 - dam;
}
/* -------------------------------------------------------------------------- */

/**
 * @file   material_damage_iterative_inline_impl.hh
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 *
 * @date   Sun Mar  8 12:54:30 2015
 *
 * @brief  Implementation of inline functions of the material damage iterative
 *
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */
/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
inline void MaterialOrthotropicDamageIterative<spatial_dimension>::
    computeDamageAndStressOnQuad(Matrix<Real> & sigma,
                                 Matrix<Real> & one_minus_D,
                                 Matrix<Real> & sqrt_one_minus_D,
                                 Matrix<Real> & damage,
                                 Matrix<Real> & first_term,
                                 Matrix<Real> & third_term) {

  // Real dmax = *(std::max_element(damage.storage(), damage.storage() +
  // spatial_dimension*spatial_dimension) );
  Real eta_effective = 0;

  // if ( (1 - dmax*dmax)  < (1 - this->eta / spatial_dimension *
  // damage.trace()) ) {

  //   eta_effective = this->spatial_dimension * dmax * dmax / damage.trace();

  // }
  // else
  eta_effective = this->eta;

  /// hydrostatic sensitivity parameter
  // Real eta = 3.;

  /// Definition of Cauchy stress based on second order damage tensor:
  /// "Anisotropic damage modelling of biaxial behaviour and rupture
  /// of concrete strucutres", Ragueneau et al., 2008, Eq. 7
  first_term.mul<false, false>(sqrt_one_minus_D, sigma);
  first_term *= sqrt_one_minus_D;

  Real second_term = 0;
  for (UInt i = 0; i < this->spatial_dimension; ++i) {
    for (UInt j = 0; j < this->spatial_dimension; ++j)
      second_term += sigma(i, j) * one_minus_D(i, j);
  }

  second_term /= (this->spatial_dimension - damage.trace());

  // for (UInt i = 0; i < this->spatial_dimension; ++i) {
  //   for (UInt j = 0; j < this->spatial_dimension; ++j)
  //     one_minus_D(i,j) *= second_term;
  // }
  one_minus_D *= second_term;

  third_term.eye(
      1. / this->spatial_dimension * sigma.trace() *
      (1 - std::min(eta_effective / (this->spatial_dimension) * damage.trace(),
                    this->max_damage)));

  sigma.copy(first_term);
  sigma -= one_minus_D;
  sigma += third_term;
}

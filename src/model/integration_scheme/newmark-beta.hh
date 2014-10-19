/**
 * @file   newmark-beta.hh
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Tue Oct 05 2010
 * @date last modification: Thu Jun 05 2014
 *
 * @brief  implementation of the  newmark-@f$\beta@f$ integration  scheme.  This
 * implementation is taken from Méthodes  numériques en mécanique des solides by
 * Alain Curnier \note{ISBN: 2-88074-247-1}
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

#ifndef __AKANTU_NEWMARK_BETA_HH__
#define __AKANTU_NEWMARK_BETA_HH__

/* -------------------------------------------------------------------------- */
#include "integration_scheme_2nd_order.hh"

/* -------------------------------------------------------------------------- */


__BEGIN_AKANTU__

/**
 * The three differentiate equations (dynamic and cinematic) are :
 *  @f{eqnarray*}{
 *   M \ddot{u}_{n+1} + C \dot{u}_{n+1} + K u_{n+1} &=& q_{n+1} \\
 *   u_{n+1} &=& u_{n} + (1 - \alpha) \Delta t \dot{u}_{n} + \alpha \Delta t \dot{u}_{n+1} + (1/2 - \alpha) \Delta t^2 \ddot{u}_n \\
 *   \dot{u}_{n+1} &=& \dot{u}_{n} + (1 - \beta) \Delta t \ddot{u}_{n} + \beta \Delta t \ddot{u}_{n+1}
 *  @f}
 *
 * Predictor:
 *  @f{eqnarray*}{
 *  u^{0}_{n+1}        &=& u_{n} +  \Delta t \dot{u}_n + \frac{\Delta t^2}{2} \ddot{u}_n \\
 *  \dot{u}^{0}_{n+1}  &=& \dot{u}_{n} +  \Delta t \ddot{u}_{n} \\
 *  \ddot{u}^{0}_{n+1} &=& \ddot{u}_{n}
 *  @f}
 *
 * Solve :
 *  @f[ (c M + d C + e K^i_{n+1}) w = = q_{n+1} - f^i_{n+1} - C \dot{u}^i_{n+1} - M \ddot{u}^i_{n+1} @f]
 *
 * Corrector :
 *  @f{eqnarray*}{
 *  \ddot{u}^{i+1}_{n+1} &=& \ddot{u}^{i}_{n+1} + c w \\
 *  \dot{u}^{i+1}_{n+1} &=& \dot{u}^{i}_{n+1} + d w \\
 *  u^{i+1}_{n+1} &=& u^{i}_{n+1} + e w
 *  @f}
 *
 * c, d and e are parameters depending on the method used to solve the equations @n
 * For acceleration : @f$ w = \delta \ddot{u}, e = \alpha \beta \Delta t^2, d = \beta \Delta t,    c = 1 @f$ @n
 * For velocity :     @f$ w = \delta \dot{u},  e = 1/\beta \Delta t,        d = 1,                 c = \alpha \Delta t @f$ @n
 * For displacement : @f$ w = \delta u,        e = 1,                       d = 1/\alpha \Delta t, c = 1/\alpha \beta \Delta t^2 @f$
 */

class NewmarkBeta : public IntegrationScheme2ndOrder {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  NewmarkBeta(Real alpha, Real beta) : beta(beta), alpha(alpha), k(0.), h(0.) {};

  ~NewmarkBeta(){};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  inline void integrationSchemePred(Real delta_t,
				    Array<Real> & u,
				    Array<Real> & u_dot,
				    Array<Real> & u_dot_dot,
				    Array<bool> & blocked_dofs) const;

  inline void integrationSchemeCorrAccel(Real delta_t,
					 Array<Real> & u,
					 Array<Real> & u_dot,
					 Array<Real> & u_dot_dot,
					 Array<bool> & blocked_dofs,
					 Array<Real> & delta) const;

  inline void integrationSchemeCorrVeloc(Real delta_t,
					 Array<Real> & u,
					 Array<Real> & u_dot,
					 Array<Real> & u_dot_dot,
					 Array<bool> & blocked_dofs,
					 Array<Real> & delta) const;

  inline void integrationSchemeCorrDispl(Real delta_t,
					 Array<Real> & u,
					 Array<Real> & u_dot,
					 Array<Real> & u_dot_dot,
					 Array<bool> & blocked_dofs,
					 Array<Real> & delta) const;

public:
  template<IntegrationSchemeCorrectorType type>
  Real getAccelerationCoefficient(Real delta_t) const;

  template<IntegrationSchemeCorrectorType type>
  Real getVelocityCoefficient(Real delta_t) const;

  template<IntegrationSchemeCorrectorType type>
  Real getDisplacementCoefficient(Real delta_t) const;

private:
  template<IntegrationSchemeCorrectorType type>
  void integrationSchemeCorr(Real delta_t,
			     Array<Real> & u,
			     Array<Real> & u_dot,
			     Array<Real> & u_dot_dot,
			     Array<bool> & blocked_dofs,
			     Array<Real> & delta) const;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  AKANTU_GET_MACRO(Beta, beta, Real);
  AKANTU_GET_MACRO(Alpha, alpha, Real);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// the @f$\beta@f$ parameter
  const Real beta;

  /// the @f$\alpha@f$ parameter
  const Real alpha;

  const Real k;
  const Real h;
};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#if defined (AKANTU_INCLUDE_INLINE_IMPL)
#  include "newmark-beta_inline_impl.cc"
#endif


/**
 * central difference method (explicit)
 * undamped stability condition :
 * @f$ \Delta t = \alpha \Delta t_{crit} = \frac{2}{\omega_{max}} \leq \min_{e} \frac{l_e}{c_e}
 *
 */
class CentralDifference : public NewmarkBeta {
public:
  CentralDifference() : NewmarkBeta(0., .5) {};
};
//#include "integration_scheme/central_difference.hh"

/// undamped trapezoidal rule (implicit)
class TrapezoidalRule2 : public NewmarkBeta {
public:
  TrapezoidalRule2() : NewmarkBeta(.5, .5) { };
};


__END_AKANTU__

#endif /* __AKANTU_NEWMARK_BETA_HH__ */

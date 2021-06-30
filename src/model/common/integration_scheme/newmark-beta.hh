/**
 * @file   newmark-beta.hh
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Tue Oct 05 2010
 * @date last modification: Wed Jan 31 2018
 *
 * @brief  implementation of the  newmark-@f$\beta@f$ integration  scheme.  This
 * implementation is taken from Méthodes  numériques en mécanique des solides by
 * Alain Curnier \note{ISBN: 2-88074-247-1}
 *
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
#include "integration_scheme_2nd_order.hh"

/* -------------------------------------------------------------------------- */

#ifndef AKANTU_NEWMARK_BETA_HH_
#define AKANTU_NEWMARK_BETA_HH_

/* -------------------------------------------------------------------------- */

namespace akantu {

/**
 * The three differentiate equations (dynamic and cinematic) are :
 *  \f{eqnarray*}{
 *   M \ddot{u}_{n+1} + C \dot{u}_{n+1} + K u_{n+1} &=& q_{n+1} \\
 *   u_{n+1} &=& u_{n} + (1 - \alpha) \Delta t \dot{u}_{n} + \alpha \Delta t
 *\dot{u}_{n+1} + (1/2 - \alpha) \Delta t^2 \ddot{u}_n \\
 *   \dot{u}_{n+1} &=& \dot{u}_{n} + (1 - \beta) \Delta t \ddot{u}_{n} + \beta
 *\Delta t \ddot{u}_{n+1}
 *  \f}
 *
 * Predictor:
 *  \f{eqnarray*}{
 *  u^{0}_{n+1}        &=& u_{n} +  \Delta t \dot{u}_n + \frac{\Delta t^2}{2}
 *\ddot{u}_n \\
 *  \dot{u}^{0}_{n+1}  &=& \dot{u}_{n} +  \Delta t \ddot{u}_{n} \\
 *  \ddot{u}^{0}_{n+1} &=& \ddot{u}_{n}
 *  \f}
 *
 * Solve :
 *  \f[ (c M + d C + e K^i_{n+1}) w = = q_{n+1} - f^i_{n+1} - C \dot{u}^i_{n+1}
 *- M \ddot{u}^i_{n+1} \f]
 *
 * Corrector :
 *  \f{eqnarray*}{
 *  \ddot{u}^{i+1}_{n+1} &=& \ddot{u}^{i}_{n+1} + c w \\
 *  \dot{u}^{i+1}_{n+1} &=& \dot{u}^{i}_{n+1} + d w \\
 *  u^{i+1}_{n+1} &=& u^{i}_{n+1} + e w
 *  \f}
 *
 * c, d and e are parameters depending on the method used to solve the equations
 *\n
 * For acceleration : \f$ w = \delta \ddot{u}, e = \alpha \beta \Delta t^2, d =
 *\beta \Delta t,    c = 1 \f$ \n
 * For velocity :     \f$ w = \delta \dot{u},  e = 1/\beta \Delta t,        d =
 *1,                 c = \alpha \Delta t \f$ \n
 * For displacement : \f$ w = \delta u,        e = 1,                       d =
 *1/\alpha \Delta t, c = 1/\alpha \beta \Delta t^2 \f$
 */

class NewmarkBeta : public IntegrationScheme2ndOrder {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  NewmarkBeta(DOFManager & dof_manager, const ID & dof_id, Real alpha = 0.,
              Real beta = 0.);

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  void predictor(Real delta_t, Array<Real> & u, Array<Real> & u_dot,
                 Array<Real> & u_dot_dot,
                 const Array<bool> & blocked_dofs) const override;

  void corrector(const SolutionType & type, Real delta_t, Array<Real> & u,
                 Array<Real> & u_dot, Array<Real> & u_dot_dot,
                 const Array<bool> & blocked_dofs,
                 const Array<Real> & delta) const override;

  void assembleJacobian(const SolutionType & type, Real delta_t) override;

public:
  Real getAccelerationCoefficient(const SolutionType & type,
                                  Real delta_t) const override;

  Real getVelocityCoefficient(const SolutionType & type,
                              Real delta_t) const override;

  Real getDisplacementCoefficient(const SolutionType & type,
                                  Real delta_t) const override;

private:
  template <SolutionType type>
  void allCorrector(Real delta_t, Array<Real> & u, Array<Real> & u_dot,
                    Array<Real> & u_dot_dot, const Array<bool> & blocked_dofs,
                    const Array<Real> & delta) const;

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
  /// the \f$\beta\f$ parameter
  Real beta;

  /// the \f$\alpha\f$ parameter
  Real alpha;

  Real k;
  Real h;

  /// last release of M matrix
  UInt m_release;

  /// last release of K matrix
  UInt k_release;

  /// last release of C matrix
  UInt c_release;
};

/**
 * central difference method (explicit)
 * undamped stability condition :
 * \f$ \Delta t = \alpha \Delta t_{crit} = \frac{2}{\omega_{max}} \leq \min_{e}
 *\frac{l_e}{c_e}\f$
 *
 */
class CentralDifference : public NewmarkBeta {
public:
  CentralDifference(DOFManager & dof_manager, const ID & dof_id)
      : NewmarkBeta(dof_manager, dof_id, 0., 1. / 2.){};

  std::vector<std::string> getNeededMatrixList() override { return {"M", "C"}; }
};

//#include "integration_scheme/central_difference.hh"

/// undamped trapezoidal rule (implicit)
class TrapezoidalRule2 : public NewmarkBeta {
public:
  TrapezoidalRule2(DOFManager & dof_manager, const ID & dof_id)
      : NewmarkBeta(dof_manager, dof_id, 1. / 2., 1. / 2.){};
};

/// Fox-Goodwin rule (implicit)
class FoxGoodwin : public NewmarkBeta {
public:
  FoxGoodwin(DOFManager & dof_manager, const ID & dof_id)
      : NewmarkBeta(dof_manager, dof_id, 1. / 6., 1. / 2.){};
};

/// Linear acceleration (implicit)
class LinearAceleration : public NewmarkBeta {
public:
  LinearAceleration(DOFManager & dof_manager, const ID & dof_id)
      : NewmarkBeta(dof_manager, dof_id, 1. / 3., 1. / 2.){};
};

/* -------------------------------------------------------------------------- */

} // namespace akantu

#endif /* AKANTU_NEWMARK_BETA_HH_ */

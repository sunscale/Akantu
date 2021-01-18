/**
 * @file   generalized_trapezoidal.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Mon Jul 04 2011
 * @date last modification: Wed Jan 31 2018
 *
 * @brief  Generalized Trapezoidal  Method.  This implementation  is taken  from
 * Méthodes  numériques en mécanique  des solides  by Alain  Curnier \note{ISBN:
 * 2-88074-247-1}
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

#ifndef AKANTU_GENERALIZED_TRAPEZOIDAL_HH_
#define AKANTU_GENERALIZED_TRAPEZOIDAL_HH_

#include "integration_scheme_1st_order.hh"

namespace akantu {

/**
 * The two differentiate equation (thermal and kinematic) are :
 * \f{eqnarray*}{
 *  C\dot{u}_{n+1} + Ku_{n+1} = q_{n+1}\\
 *  u_{n+1} = u_{n} + (1-\alpha) \Delta t \dot{u}_{n} + \alpha \Delta t
 *\dot{u}_{n+1}
 * \f}
 *
 * To solve it :
 * Predictor :
 * \f{eqnarray*}{
 * u^0_{n+1} &=& u_{n} + (1-\alpha) \Delta t v_{n} \\
 * \dot{u}^0_{n+1} &=& \dot{u}_{n}
 * \f}
 *
 * Solve :
 * \f[ (a C + b K^i_{n+1}) w = q_{n+1} - f^i_{n+1} - C \dot{u}^i_{n+1} \f]
 *
 * Corrector :
 * \f{eqnarray*}{
 * \dot{u}^{i+1}_{n+1} &=& \dot{u}^{i}_{n+1} + a w \\
 * u^{i+1}_{n+1} &=& u^{i}_{n+1} + b w
 * \f}
 *
 * a and b depends on the resolution method : temperature (u) or temperature
 *rate (\f$\dot{u}\f$)
 *
 * For temperature : \f$ w = \delta u, a = 1 / (\alpha \Delta t) , b = 1 \f$ @n
 * For temperature rate : \f$ w = \delta \dot{u}, a = 1, b = \alpha \Delta t \f$
 */
class GeneralizedTrapezoidal : public IntegrationScheme1stOrder {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  GeneralizedTrapezoidal(DOFManager & dof_manager, const ID & dof_id,
                         Real alpha = 0);

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  void predictor(Real delta_t, Array<Real> & u, Array<Real> & u_dot,
                 const Array<bool> & blocked_dofs) const override;

  void corrector(const SolutionType & type, Real delta_t, Array<Real> & u,
                 Array<Real> & u_dot, const Array<bool> & blocked_dofs,
                 const Array<Real> & delta) const override;

  void assembleJacobian(const SolutionType & type, Real delta_t) override;

public:
  /// the coeffichent \f$ b \f$ in the description
  Real getTemperatureCoefficient(const SolutionType & type,
                                 Real delta_t) const override;

  /// the coeffichent \f$ a \f$ in the description
  Real getTemperatureRateCoefficient(const SolutionType & type,
                                     Real delta_t) const override;

private:
  template <SolutionType type>
  void allCorrector(Real delta_t, Array<Real> & u, Array<Real> & u_dot,
                    const Array<bool> & blocked_dofs,
                    const Array<Real> & delta) const;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  AKANTU_GET_MACRO(Alpha, alpha, Real);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// the \f$\alpha\f$ parameter
  Real alpha;

  /// last release of K matrix
  UInt k_release;
};

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */

/**
 * Forward Euler (explicit) -> condition on delta_t
 */
class ForwardEuler : public GeneralizedTrapezoidal {
public:
  ForwardEuler(DOFManager & dof_manager, const ID & dof_id)
      : GeneralizedTrapezoidal(dof_manager, dof_id, 0.){};

  std::vector<std::string> getNeededMatrixList() override { return {"M"}; }
};

/**
 * Trapezoidal rule (implicit), midpoint rule or Crank-Nicolson
 */
class TrapezoidalRule1 : public GeneralizedTrapezoidal {
public:
  TrapezoidalRule1(DOFManager & dof_manager, const ID & dof_id)
      : GeneralizedTrapezoidal(dof_manager, dof_id, .5){};
};

/**
 * Backward Euler (implicit)
 */
class BackwardEuler : public GeneralizedTrapezoidal {
public:
  BackwardEuler(DOFManager & dof_manager, const ID & dof_id)
      : GeneralizedTrapezoidal(dof_manager, dof_id, 1.){};
};

/* -------------------------------------------------------------------------- */

} // namespace akantu

#endif /* AKANTU_GENERALIZED_TRAPEZOIDAL_HH_ */

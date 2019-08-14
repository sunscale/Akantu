/**
 * @file   integration_scheme_2nd_order.hh
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Wed Jan 31 2018
 *
 * @brief  Interface of the integrator of second order
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
#include "aka_array.hh"
#include "integration_scheme.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_INTEGRATION_SCHEME_2ND_ORDER_HH__
#define __AKANTU_INTEGRATION_SCHEME_2ND_ORDER_HH__

namespace akantu {
class SparseMatrix;
}

namespace akantu {

class IntegrationScheme2ndOrder : public IntegrationScheme {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  IntegrationScheme2ndOrder(DOFManager & dof_manager, const ID & dof_id)
      : IntegrationScheme(dof_manager, dof_id, 2){};

  ~IntegrationScheme2ndOrder() override = default;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// get list of needed matrices
  std::vector<std::string> getNeededMatrixList() override;

  /// generic interface of a predictor
  void predictor(Real delta_t) override;

  /// generic interface of a corrector
  void corrector(const SolutionType & type, Real delta_t) override;

  void assembleResidual(bool is_lumped) override;

protected:
  /// generic interface of a predictor of 2nd order
  virtual void predictor(Real delta_t, Array<Real> & u, Array<Real> & u_dot,
                         Array<Real> & u_dot_dot,
                         const Array<bool> & blocked_dofs) const = 0;

  /// generic interface of a corrector of 2nd order
  virtual void corrector(const SolutionType & type, Real delta_t,
                         Array<Real> & u, Array<Real> & u_dot,
                         Array<Real> & u_dot_dot,
                         const Array<bool> & blocked_dofs,
                         const Array<Real> & delta) const = 0;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
protected:
  virtual Real getAccelerationCoefficient(const SolutionType & type,
                                          Real delta_t) const = 0;

  virtual Real getVelocityCoefficient(const SolutionType & type,
                                      Real delta_t) const = 0;

  virtual Real getDisplacementCoefficient(const SolutionType & type,
                                          Real delta_t) const = 0;

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
};

} // namespace akantu

#include "newmark-beta.hh"

#endif /* __AKANTU_INTEGRATION_SCHEME_2ND_ORDER_HH__ */

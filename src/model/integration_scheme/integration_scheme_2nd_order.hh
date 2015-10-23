/**
 * @file   integration_scheme_2nd_order.hh
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Aug 04 2010
 * @date last modification: Thu Jun 05 2014
 *
 * @brief  Interface of the integrator of second order
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
#include "aka_array.hh"
#include "integration_scheme.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_INTEGRATION_SCHEME_2ND_ORDER_HH__
#define __AKANTU_INTEGRATION_SCHEME_2ND_ORDER_HH__

namespace akantu {
  class SparseMatrix;
}

__BEGIN_AKANTU__

class IntegrationScheme2ndOrder : public IntegrationScheme {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  IntegrationScheme2ndOrder(DOFManager & dof_manager) : IntegrationScheme(dof_manager, 2){};

  virtual ~IntegrationScheme2ndOrder(){};
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// generic interface of a predictor
  virtual void predictor(const ID & dof_id, Real delta_t);

  /// generic interface of a corrector
  virtual void corrector(const SolutionType & type, const ID & dof_id,
                         Real delta_t);

  /// generic interface of a predictor of 2nd order
  virtual void predictor(Real delta_t, Array<Real> & u, Array<Real> & u_dot,
                         Array<Real> & u_dot_dot,
                         const Array<bool> & blocked_dofs) const = 0;

  /// generic interface of a corrector of 2nd order
  virtual void corrector(const SolutionType & type,
                         Real delta_t, Array<Real> & u, Array<Real> & u_dot,
                         Array<Real> & u_dot_dot,
                         const Array<bool> & blocked_dofs,
                         const Array<Real> & delta) const = 0;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
protected:
  virtual Real
  getAccelerationCoefficient(const SolutionType & type,
                             Real delta_t) const;

  virtual Real
  getVelocityCoefficient(const SolutionType & type,
                         Real delta_t) const;

  virtual Real
  getDisplacementCoefficient(const SolutionType & type,
                             Real delta_t) const;

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
};

__END_AKANTU__

#include "newmark-beta.hh"

#endif /* __AKANTU_INTEGRATION_SCHEME_2ND_ORDER_HH__ */

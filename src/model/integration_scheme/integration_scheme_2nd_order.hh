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

#ifndef __AKANTU_INTEGRATION_SCHEME_2ND_ORDER_HH__
#define __AKANTU_INTEGRATION_SCHEME_2ND_ORDER_HH__

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "aka_array.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

class IntegrationScheme2ndOrder {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  enum IntegrationSchemeCorrectorType {
    _acceleration_corrector,
    _velocity_corrector,
    _displacement_corrector
  };

  virtual ~IntegrationScheme2ndOrder() {};
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  virtual void integrationSchemePred(Real delta_t,
				     Array<Real> & u,
				     Array<Real> & u_dot,
				     Array<Real> & u_dot_dot,
				     Array<bool> & blocked_dofs) const = 0;

  virtual void integrationSchemeCorrDispl(Real delta_t,
					  Array<Real> & u,
					  Array<Real> & u_dot,
					  Array<Real> & u_dot_dot,
					  Array<bool> & blocked_dofs,
					  Array<Real> & delta) const = 0;

  virtual void integrationSchemeCorrAccel(Real delta_t,
					  Array<Real> & u,
					  Array<Real> & u_dot,
					  Array<Real> & u_dot_dot,
					  Array<bool> & blocked_dofs,
					  Array<Real> & delta) const = 0;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:

};

__END_AKANTU__

#include "newmark-beta.hh"

#endif /* __AKANTU_INTEGRATION_SCHEME_2ND_ORDER_HH__ */

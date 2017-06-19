/**
 * @file   pseudo_time.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Wed Feb 17 09:46:05 2016
 *
 * @brief  Pseudo time integration scheme
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "integration_scheme.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_PSEUDO_TIME_HH__
#define __AKANTU_PSEUDO_TIME_HH__

namespace akantu {

class PseudoTime : public IntegrationScheme {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  PseudoTime(DOFManager & dof_manager, const ID & dof_id);
  virtual ~PseudoTime() {}

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// generic interface of a predictor
  virtual void predictor(Real delta_t);

  /// generic interface of a corrector
  virtual void corrector(const SolutionType & type, Real delta_t);

  /// assemble the jacobian matrix
  virtual void assembleJacobian(const SolutionType & type, Real delta_t);

  /// assemble the residual
  virtual void assembleResidual(bool is_lumped);

protected:
  /// last release of K matrix
  UInt k_release;
};

} // akantu

#endif /* __AKANTU_PSEUDO_TIME_HH__ */

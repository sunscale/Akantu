/**
 * @file   integration_scheme.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Wed Jan 31 2018
 *
 * @brief  This class is just a base class for the integration schemes
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
#include "aka_common.hh"
#include "parsable.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_INTEGRATION_SCHEME_HH__
#define __AKANTU_INTEGRATION_SCHEME_HH__

namespace akantu {
class DOFManager;
}

namespace akantu {

class IntegrationScheme : public Parsable {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  enum SolutionType {
    _not_defined = -1,
    _displacement = 0,
    _temperature = 0,
    _damage = 0,
    _velocity = 1,
    _temperature_rate = 1,
    _acceleration = 2,
  };

  IntegrationScheme(DOFManager & dof_manager, const ID & dof_id, UInt order);
  ~IntegrationScheme() override = default;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// generic interface of a predictor
  virtual void predictor(Real delta_t) = 0;

  /// generic interface of a corrector
  virtual void corrector(const SolutionType & type, Real delta_t) = 0;

  /// assemble the jacobian matrix
  virtual void assembleJacobian(const SolutionType & type, Real delta_t) = 0;

  /// assemble the residual
  virtual void assembleResidual(bool is_lumped) = 0;

  /// returns a list of needed matrices
  virtual std::vector<std::string> getNeededMatrixList() = 0;

  /// store dofs info (beginning of steps)
  virtual void store();

  /// restore dofs (solve failed)
  virtual void restore();

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// return the order of the integration scheme
  UInt getOrder() const;

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// The underlying DOFManager
  DOFManager & dof_manager;

  /// The id of the dof treated by this integration scheme.
  ID dof_id;

  /// The order of the integrator
  UInt order;

  /// last release of M matrix
  UInt m_release{UInt(-1)};

  /// stores the values at begining of solve
  std::vector<std::unique_ptr<Array<Real>>> u_store;
};

/* -------------------------------------------------------------------------- */
// std::ostream & operator<<(std::ostream & stream,
//                           const IntegrationScheme::SolutionType & type);
std::istream & operator>>(std::istream & stream,
                          IntegrationScheme::SolutionType & type);
/* -------------------------------------------------------------------------- */

} // namespace akantu

#endif /* __AKANTU_INTEGRATION_SCHEME_HH__ */

/**
 * @file   integration_scheme.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Tue Aug 18 2015
 * @date last modification: Wed Jan 31 2018
 *
 * @brief  Common interface to all interface schemes
 *
 *
 * Copyright (©) 2015-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "integration_scheme.hh"
#include "dof_manager.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
IntegrationScheme::IntegrationScheme(DOFManager & dof_manager,
                                     const ID & dof_id, UInt order)
    : Parsable(ParserType::_integration_scheme, dof_id),
      dof_manager(dof_manager), dof_id(dof_id), order(order), u_store(order + 1) {}

/* -------------------------------------------------------------------------- */
/// standard input stream operator for SolutionType
std::istream & operator>>(std::istream & stream,
                          IntegrationScheme::SolutionType & type) {
  std::string str;
  stream >> str;
  if (str == "displacement") {
    type = IntegrationScheme::_displacement;
  } else if (str == "temperature") {
    type = IntegrationScheme::_temperature;
  } else if (str == "velocity") {
    type = IntegrationScheme::_velocity;
  } else if (str == "temperature_rate") {
    type = IntegrationScheme::_temperature_rate;
  } else if (str == "acceleration") {
    type = IntegrationScheme::_acceleration;
  } else if (str == "damage") {
    type = IntegrationScheme::_damage;
  } else {
    stream.setstate(std::ios::failbit);
  }

  return stream;
}

/* -------------------------------------------------------------------------- */
void IntegrationScheme::store() {
  for (auto data : enumerate(u_store)) {
    auto o = std::get<0>(data);
    auto & u_store = std::get<1>(data);
    auto & u_o = dof_manager.getDOFsDerivatives(dof_id, o);
    if (not u_store) {
      u_store = std::make_unique<Array<Real>>(
          u_o, "integration_scheme_store:" + dof_id + ":" + std::to_string(o));
    } else {
      u_store->copy(u_o);
    }
  }
}

/* -------------------------------------------------------------------------- */
void IntegrationScheme::restore() {
  for (auto o : arange(order)) {
    auto & u_o = dof_manager.getDOFsDerivatives(dof_id, o);
    u_o.copy(*u_store[o]);
  }
}


} // namespace akantu

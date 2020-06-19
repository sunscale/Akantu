/**
 * @file   ntn_friction_tmpl.hh
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 *
 * @date creation: Tue Dec 02 2014
 * @date last modification: Fri Feb 23 2018
 *
 * @brief
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
// simtools
#include "ntn_contact.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
template <template <class> class FrictionLaw, class Regularisation>
NTNFriction<FrictionLaw, Regularisation>::NTNFriction(
    NTNBaseContact & contact, const ID & id, const MemoryID & memory_id)
    : FrictionLaw<Regularisation>(contact, id, memory_id) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <template <class> class FrictionLaw, class Regularisation>
void NTNFriction<FrictionLaw, Regularisation>::applyFrictionTraction() {
  AKANTU_DEBUG_IN();

  NTNContact & ntn_contact = dynamic_cast<NTNContact &>(this->contact);
  SolidMechanicsModel & model = ntn_contact.getModel();
  Array<Real> & residual = model.getInternalForce();
  UInt dim = model.getSpatialDimension();

  const SynchronizedArray<UInt> & masters = ntn_contact.getMasters();
  const SynchronizedArray<UInt> & slaves = ntn_contact.getSlaves();
  const SynchronizedArray<Real> & l_boundary_slaves =
      ntn_contact.getLumpedBoundarySlaves();
  const SynchronizedArray<Real> & l_boundary_masters =
      ntn_contact.getLumpedBoundaryMasters();

  UInt nb_contact_nodes = ntn_contact.getNbContactNodes();
  for (UInt n = 0; n < nb_contact_nodes; ++n) {
    UInt master = masters(n);
    UInt slave = slaves(n);

    for (UInt d = 0; d < dim; ++d) {
      residual(master, d) +=
          l_boundary_masters(n) * this->friction_traction(n, d);
      residual(slave, d) -=
          l_boundary_slaves(n) * this->friction_traction(n, d);
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <template <class> class FrictionLaw, class Regularisation>
void NTNFriction<FrictionLaw, Regularisation>::printself(std::ostream & stream,
                                                         int indent) const {
  AKANTU_DEBUG_IN();
  std::string space;
  for (Int i = 0; i < indent; i++, space += AKANTU_INDENT)
    ;

  stream << space << "NTNFriction [" << std::endl;
  FrictionLaw<Regularisation>::printself(stream, ++indent);
  stream << space << "]" << std::endl;

  AKANTU_DEBUG_OUT();
}

} // namespace akantu

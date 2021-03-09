/**
 * @file   ntrf_friction_tmpl.hh
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 *
 * @date creation: Tue Dec 02 2014
 * @date last modification: Fri Feb 23 2018
 *
 * @brief  implementation of node to rigid flat interface friction
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
//#include "ntrf_friction.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
template <template <class> class FrictionLaw, class Regularisation>
NTRFFriction<FrictionLaw, Regularisation>::NTRFFriction(
    NTNBaseContact & contact, const ID & id)
    : FrictionLaw<Regularisation>(contact, id) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <template <class> class FrictionLaw, class Regularisation>
void NTRFFriction<FrictionLaw, Regularisation>::printself(std::ostream & stream,
                                                          int indent) const {
  AKANTU_DEBUG_IN();
  std::string space;
  for (Int i = 0; i < indent; i++, space += AKANTU_INDENT)
    ;

  stream << space << "NTRFFriction [" << std::endl;
  FrictionLaw<Regularisation>::printself(stream, ++indent);
  stream << space << "]" << std::endl;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/*
void NTRFFriction::addDumpFieldToDumper(const std::string & dumper_name,
                    const std::string & field_id) {
  AKANTU_DEBUG_IN();

#ifdef AKANTU_USE_IOHELPER
  //  const SynchronizedArray<UInt> * nodal_filter =
&(this->contact.getSlaves());

  if(field_id == "is_sticking") {
    this->internalAddDumpFieldToDumper(dumper_name,
                   field_id,
                   new
DumperIOHelper::NodalField<bool>(this->is_sticking.getArray()));
  }
  else if(field_id == "frictional_strength") {
    this->internalAddDumpFieldToDumper(dumper_name,
                   field_id,
                   new
DumperIOHelper::NodalField<Real>(this->frictional_strength.getArray()));
  }
  else if(field_id == "friction_traction") {
    this->internalAddDumpFieldToDumper(dumper_name,
                   field_id,
                   new
DumperIOHelper::NodalField<Real>(this->friction_traction.getArray()));
  }
  else if(field_id == "slip") {
    this->internalAddDumpFieldToDumper(dumper_name,
                   field_id,
                   new DumperIOHelper::NodalField<Real>(this->slip.getArray()));
  }
  else {
    this->contact.addDumpFieldToDumper(dumper_name, field_id);
  }

#endif

  AKANTU_DEBUG_OUT();
}
*/

} // namespace akantu

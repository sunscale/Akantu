/**
 * @file   dumpable_inline_impl.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author David Simon Kammer <david.kammer@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Nov 13 2013
 * @date last modification: Wed Nov 08 2017
 *
 * @brief  Implementation of the Dumpable class
 *
 *
 * Copyright (©) 2014-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

#ifndef AKANTU_DUMPABLE_INLINE_IMPL_HH_
#define AKANTU_DUMPABLE_INLINE_IMPL_HH_

/* -------------------------------------------------------------------------- */
#ifdef AKANTU_USE_IOHELPER
#include "dumper_elemental_field.hh"
#include "dumper_nodal_field.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
template <class T>
inline void Dumpable::registerDumper(const std::string & dumper_name,
                                     const std::string & file_name,
                                     const bool is_default) {

  if (this->dumpers.find(dumper_name) != this->dumpers.end()) {
    AKANTU_DEBUG_INFO("Dumper " + dumper_name + "is already registered.");
  }

  std::string name = file_name;
  if (name.empty()) {
    name = dumper_name;
  }

  this->dumpers[dumper_name] = std::make_shared<T>(name);

  if (is_default) {
    this->default_dumper = dumper_name;
  }
}

/* -------------------------------------------------------------------------- */
template <typename T>
inline void Dumpable::addDumpFieldExternal(const std::string & field_id,
                                           const Array<T> & field) {
  this->addDumpFieldExternalToDumper<T>(this->default_dumper, field_id, field);
}

/* -------------------------------------------------------------------------- */
template <typename T>
inline void
Dumpable::addDumpFieldExternalToDumper(const std::string & dumper_name,
                                       const std::string & field_id,
                                       const Array<T> & field) {
  auto field_cont = std::make_shared<dumpers::NodalField<T>>(field);
  DumperIOHelper & dumper = this->getDumper(dumper_name);
  dumper.registerField(field_id, field_cont);
}

/* -------------------------------------------------------------------------- */
template <typename T>
inline void Dumpable::addDumpFieldExternal(const std::string & field_id,
                                           const ElementTypeMapArray<T> & field,
                                           UInt spatial_dimension,
                                           GhostType ghost_type,
                                           ElementKind element_kind) {
  this->addDumpFieldExternalToDumper(this->default_dumper, field_id, field,
                                     spatial_dimension, ghost_type,
                                     element_kind);
}

/* -------------------------------------------------------------------------- */
template <typename T>
inline void Dumpable::addDumpFieldExternalToDumper(
    const std::string & dumper_name, const std::string & field_id,
    const ElementTypeMapArray<T> & field, UInt spatial_dimension,
    GhostType ghost_type, ElementKind element_kind) {

  std::shared_ptr<dumpers::Field> field_cont;
#if defined(AKANTU_IGFEM)
  if (element_kind == _ek_igfem) {
    field_cont = std::make_shared<dumpers::IGFEMElementalField<T>>(
        field, spatial_dimension, ghost_type, element_kind);
  } else
#endif
    field_cont = std::make_shared<dumpers::ElementalField<T>>(
        field, spatial_dimension, ghost_type, element_kind);
  DumperIOHelper & dumper = this->getDumper(dumper_name);
  dumper.registerField(field_id, field_cont);
}

/* -------------------------------------------------------------------------- */
template <class T>
inline T & Dumpable::getDumper(const std::string & dumper_name) {
  DumperIOHelper & dumper = this->getDumper(dumper_name);

  try {
    auto & templated_dumper = aka::as_type<T>(dumper);
    return templated_dumper;
  } catch (std::bad_cast &) {
    AKANTU_EXCEPTION("Dumper " << dumper_name << " is not of type: "
                               << debug::demangle(typeid(T).name()));
  }
}

/* -------------------------------------------------------------------------- */

} // namespace akantu

#endif

#endif /* AKANTU_DUMPABLE_INLINE_IMPL_HH_ */

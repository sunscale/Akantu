/**
 * @file   group_manager_inline_impl.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Dana Christen <dana.christen@gmail.com>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Nov 13 2013
 * @date last modification: Sun Dec 03 2017
 *
 * @brief  Stores information relevent to the notion of domain boundary and
 * surfaces.
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

/* -------------------------------------------------------------------------- */
#include "dumper_field.hh"
#include "element_group.hh"
#include "element_type_map_filter.hh"
#ifdef AKANTU_USE_IOHELPER
#include "dumper_nodal_field.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
template <typename T, template <bool> class dump_type>
std::shared_ptr<dumpers::Field>
GroupManager::createElementalField(const ElementTypeMapArray<T> & field,
                                   const std::string & group_name,
                                   UInt spatial_dimension, ElementKind kind,
                                   ElementTypeMap<UInt> nb_data_per_elem) {

  const ElementTypeMapArray<T> * field_ptr = &field;
  if (field_ptr == nullptr) {
    return nullptr;
  }
  if (group_name == "all") {
    return this->createElementalField<dump_type<false>>(
        field, group_name, spatial_dimension, kind, nb_data_per_elem);
  }
  return this->createElementalFilteredField<dump_type<true>>(
      field, group_name, spatial_dimension, kind, nb_data_per_elem);
}

/* -------------------------------------------------------------------------- */
template <typename T, template <class> class T2,
          template <class, template <class> class, bool> class dump_type>
std::shared_ptr<dumpers::Field>
GroupManager::createElementalField(const ElementTypeMapArray<T> & field,
                                   const std::string & group_name,
                                   UInt spatial_dimension, ElementKind kind,
                                   ElementTypeMap<UInt> nb_data_per_elem) {

  const ElementTypeMapArray<T> * field_ptr = &field;
  if (field_ptr == nullptr) {
    return nullptr;
  }
  if (group_name == "all") {
    return this->createElementalField<dump_type<T, T2, false>>(
        field, group_name, spatial_dimension, kind, nb_data_per_elem);
  }
  return this->createElementalFilteredField<dump_type<T, T2, true>>(
      field, group_name, spatial_dimension, kind, nb_data_per_elem);
}

/* -------------------------------------------------------------------------- */
template <typename T, template <typename T2, bool filtered>
                      class dump_type> ///< type of InternalMaterialField
std::shared_ptr<dumpers::Field>
GroupManager::createElementalField(const ElementTypeMapArray<T> & field,
                                   const std::string & group_name,
                                   UInt spatial_dimension, ElementKind kind,
                                   ElementTypeMap<UInt> nb_data_per_elem) {
  const ElementTypeMapArray<T> * field_ptr = &field;

  if (field_ptr == nullptr) {
    return nullptr;
  }
  if (group_name == "all") {
    return this->createElementalField<dump_type<T, false>>(
        field, group_name, spatial_dimension, kind, nb_data_per_elem);
  }
  return this->createElementalFilteredField<dump_type<T, true>>(
      field, group_name, spatial_dimension, kind, nb_data_per_elem);
}

/* -------------------------------------------------------------------------- */
template <typename dump_type, typename field_type>
std::shared_ptr<dumpers::Field> GroupManager::createElementalField(
    const field_type & field, const std::string & group_name,
    UInt spatial_dimension, ElementKind kind,
    const ElementTypeMap<UInt> & nb_data_per_elem) {
  const field_type * field_ptr = &field;
  if (field_ptr == nullptr) {
    return nullptr;
  }
  if (group_name != "all") {
    throw;
  }

  auto dumper =
      std::make_shared<dump_type>(field, spatial_dimension, _not_ghost, kind);
  dumper->setNbDataPerElem(nb_data_per_elem);
  return dumper;
}

/* -------------------------------------------------------------------------- */
template <typename dump_type, typename field_type>
std::shared_ptr<dumpers::Field> GroupManager::createElementalFilteredField(
    const field_type & field, const std::string & group_name,
    UInt spatial_dimension, ElementKind kind,
    ElementTypeMap<UInt> nb_data_per_elem) {

  const field_type * field_ptr = &field;
  if (field_ptr == nullptr) {
    return nullptr;
  }
  if (group_name == "all") {
    throw;
  }

  using T = typename field_type::type;
  ElementGroup & group = this->getElementGroup(group_name);
  UInt dim = group.getDimension();
  if (dim != spatial_dimension) {
    throw;
  }
  const ElementTypeMapArray<UInt> & elemental_filter = group.getElements();

  auto * filtered = new ElementTypeMapArrayFilter<T>(field, elemental_filter,
                                                     nb_data_per_elem);

  auto dumper = std::make_shared<dump_type>(*filtered, dim, _not_ghost, kind);
  dumper->setNbDataPerElem(nb_data_per_elem);

  return dumper;
}

/* -------------------------------------------------------------------------- */
template <typename type, bool flag, template <class, bool> class ftype>
std::shared_ptr<dumpers::Field>
GroupManager::createNodalField(const ftype<type, flag> * field,
                               const std::string & group_name,
                               UInt padding_size) {
  return createStridedNodalField(field, group_name, 0, 0, padding_size);
}

/* -------------------------------------------------------------------------- */
template <typename type, bool flag, template <class, bool> class ftype>
std::shared_ptr<dumpers::Field>
GroupManager::createStridedNodalField(const ftype<type, flag> * field,
                                      const std::string & group_name, UInt size,
                                      UInt stride, UInt padding_size) {
  if (not field) {
    return nullptr;
  }

  if (group_name == "all") {
    using DumpType = typename dumpers::NodalField<type, false>;
    auto dumper = std::make_shared<DumpType>(*field, size, stride);
    dumper->setPadding(padding_size);
    return dumper;
  }

  ElementGroup & group = this->getElementGroup(group_name);
  const Array<UInt> * nodal_filter = &(group.getNodeGroup().getNodes());
  using DumpType = typename dumpers::NodalField<type, true>;
  auto dumper = std::make_shared<DumpType>(*field, size, stride, nodal_filter);
  dumper->setPadding(padding_size);
  return dumper;
}

/* -------------------------------------------------------------------------- */

} // namespace akantu

#endif

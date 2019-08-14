/**
 * @file   mesh_data_tmpl.hh
 *
 * @author Dana Christen <dana.christen@gmail.com>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri May 03 2013
 * @date last modification: Tue Feb 20 2018
 *
 * @brief  Stores generic data loaded from the mesh file
 *
 * @section LICENSE
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
#include "mesh_data.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MESH_DATA_TMPL_HH__
#define __AKANTU_MESH_DATA_TMPL_HH__

namespace akantu {

#define AKANTU_MESH_DATA_OSTREAM(r, name, elem)                                \
  case MeshDataTypeCode::BOOST_PP_TUPLE_ELEM(2, 0, elem): {                    \
    stream << BOOST_PP_STRINGIZE(BOOST_PP_TUPLE_ELEM(2, 1, elem));             \
    break;                                                                     \
  }

inline std::ostream & operator<<(std::ostream & stream,
                                 const MeshDataTypeCode & type_code) {
  switch (type_code) {
    BOOST_PP_SEQ_FOR_EACH(AKANTU_MESH_DATA_OSTREAM, name,
                          AKANTU_MESH_DATA_TYPES)
  default:
    stream << "(unknown type)";
  }
  return stream;
}
#undef AKANTU_MESH_DATA_OSTREAM

#define MESH_DATA_GET_TYPE(r, data, type)                                      \
  template <>                                                                  \
  inline MeshDataTypeCode                                                      \
  MeshData::getTypeCode<BOOST_PP_TUPLE_ELEM(2, 1, type)>() const {             \
    return MeshDataTypeCode::BOOST_PP_TUPLE_ELEM(2, 0, type);                  \
  }

/* -------------------------------------------------------------------------- */
// get the type of the data stored in elemental_data
template <typename T> inline MeshDataTypeCode MeshData::getTypeCode() const {
  AKANTU_ERROR("Type " << debug::demangle(typeid(T).name())
                       << " not implemented by MeshData.");
}

/* -------------------------------------------------------------------------- */
BOOST_PP_SEQ_FOR_EACH(MESH_DATA_GET_TYPE, void, AKANTU_MESH_DATA_TYPES)
#undef MESH_DATA_GET_TYPE

inline MeshDataTypeCode MeshData::getTypeCode(const ID & name,
                                              MeshDataType type) const {
  auto it = typecode_map.at(type).find(name);
  if (it == typecode_map.at(type).end())
    AKANTU_EXCEPTION("No dataset named " << name << " found.");
  return it->second;
}

/* -------------------------------------------------------------------------- */
//  Register new elemental data templated (and alloc data) with check if the
//  name is new
template <typename T>
ElementTypeMapArray<T> & MeshData::registerElementalData(const ID & name) {
  auto it = elemental_data.find(name);
  if (it == elemental_data.end()) {
    return allocElementalData<T>(name);
  } else {
    AKANTU_DEBUG_INFO("Data named " << name << " already registered.");
    return getElementalData<T>(name);
  }
}

/* -------------------------------------------------------------------------- */
//  Register new elemental data of a given MeshDataTypeCode with check if the
//  name is new
#define AKANTU_MESH_DATA_CASE_MACRO(r, name, elem)                             \
  case MeshDataTypeCode::BOOST_PP_TUPLE_ELEM(2, 0, elem): {                    \
    registerElementalData<BOOST_PP_TUPLE_ELEM(2, 1, elem)>(name);              \
    break;                                                                     \
  }

inline void MeshData::registerElementalData(const ID & name,
                                            MeshDataTypeCode type) {
  switch (type) {
    BOOST_PP_SEQ_FOR_EACH(AKANTU_MESH_DATA_CASE_MACRO, name,
                          AKANTU_MESH_DATA_TYPES)
  default:
    AKANTU_ERROR("Type " << type << "not implemented by MeshData.");
  }
}
#undef AKANTU_MESH_DATA_CASE_MACRO

/* -------------------------------------------------------------------------- */
/// Register new elemental data (and alloc data)
template <typename T>
ElementTypeMapArray<T> & MeshData::allocElementalData(const ID & name) {
  auto dataset =
      std::make_unique<ElementTypeMapArray<T>>(name, _id, _memory_id);
  auto * dataset_typed = dataset.get();
  elemental_data[name] = std::move(dataset);
  typecode_map[MeshDataType::_elemental][name] = getTypeCode<T>();
  return *dataset_typed;
}

/* -------------------------------------------------------------------------- */
//  Register new nodal data templated (and alloc data) with check if the
//  name is new
template <typename T>
Array<T> & MeshData::registerNodalData(const ID & name, UInt nb_components) {
  auto it = nodal_data.find(name);
  if (it == nodal_data.end()) {
    return allocNodalData<T>(name, nb_components);
  } else {
    AKANTU_DEBUG_INFO("Data named " << name << " already registered.");
    return getNodalData<T>(name);
  }
}

/* -------------------------------------------------------------------------- */
//  Register new elemental data of a given MeshDataTypeCode with check if the
//  name is new
#define AKANTU_MESH_NODAL_DATA_CASE_MACRO(r, name, elem)                       \
  case MeshDataTypeCode::BOOST_PP_TUPLE_ELEM(2, 0, elem): {                    \
    registerNodalData<BOOST_PP_TUPLE_ELEM(2, 1, elem)>(name, nb_components);   \
    break;                                                                     \
  }

inline void MeshData::registerNodalData(const ID & name, UInt nb_components,
                                        MeshDataTypeCode type) {
  switch (type) {
    BOOST_PP_SEQ_FOR_EACH(AKANTU_MESH_NODAL_DATA_CASE_MACRO, name,
                          AKANTU_MESH_DATA_TYPES)
  default:
    AKANTU_ERROR("Type " << type << "not implemented by MeshData.");
  }
}
#undef AKANTU_MESH_NODAL_DATA_CASE_MACRO

/* -------------------------------------------------------------------------- */
/// Register new elemental data (and alloc data)
template <typename T>
Array<T> & MeshData::allocNodalData(const ID & name, UInt nb_components) {
  auto dataset =
      std::make_unique<Array<T>>(0, nb_components, T(), _id + ":" + name);
  auto * dataset_typed = dataset.get();
  nodal_data[name] = std::move(dataset);
  typecode_map[MeshDataType::_nodal][name] = getTypeCode<T>();
  return *dataset_typed;
}

/* -------------------------------------------------------------------------- */
template <typename T>
const Array<T> & MeshData::getNodalData(const ID & name) const {
  auto it = nodal_data.find(name);
  if (it == nodal_data.end())
    AKANTU_EXCEPTION("No nodal dataset named " << name << " found.");
  return aka::as_type<Array<T>>(*(it->second.get()));
}

/* -------------------------------------------------------------------------- */
// Get an existing elemental data
template <typename T>
Array<T> & MeshData::getNodalData(const ID & name, UInt nb_components) {
  auto it = nodal_data.find(name);
  if (it == nodal_data.end())
    return allocNodalData<T>(name, nb_components);

  return aka::as_type<Array<T>>(*(it->second.get()));
}

/* -------------------------------------------------------------------------- */
template <typename T>
const ElementTypeMapArray<T> &
MeshData::getElementalData(const ID & name) const {
  auto it = elemental_data.find(name);
  if (it == elemental_data.end())
    AKANTU_EXCEPTION("No dataset named " << name << " found.");
  return aka::as_type<ElementTypeMapArray<T>>(*(it->second.get()));
}

/* -------------------------------------------------------------------------- */
// Get an existing elemental data
template <typename T>
ElementTypeMapArray<T> & MeshData::getElementalData(const ID & name) {
  auto it = elemental_data.find(name);
  if (it == elemental_data.end()) {
    return allocElementalData<T>(name);
  }

  return aka::as_type<ElementTypeMapArray<T>>(*(it->second.get()));
}

/* -------------------------------------------------------------------------- */
template <typename T>
bool MeshData::hasData(const ID & name, const ElementType & elem_type,
                       const GhostType & ghost_type) const {
  auto it = elemental_data.find(name);
  if (it == elemental_data.end())
    return false;

  auto & elem_map = aka::as_type<ElementTypeMapArray<T>>(*(it->second));
  return elem_map.exists(elem_type, ghost_type);
}

/* -------------------------------------------------------------------------- */
inline bool MeshData::hasData(const ID & name, MeshDataType type) const {
  if (type == MeshDataType::_elemental) {
    auto it = elemental_data.find(name);
    return (it != elemental_data.end());
  }

  if (type == MeshDataType::_nodal) {
    auto it = nodal_data.find(name);
    return (it != nodal_data.end());
  }

  return false;
}

/* -------------------------------------------------------------------------- */
inline bool MeshData::hasData(MeshDataType type) const {
  switch (type) {
  case MeshDataType::_elemental:
    return (not elemental_data.empty());
  case MeshDataType::_nodal:
    return (not nodal_data.empty());
  }

  return false;
}

/* -------------------------------------------------------------------------- */
template <typename T>
const Array<T> &
MeshData::getElementalDataArray(const ID & name, const ElementType & elem_type,
                                const GhostType & ghost_type) const {
  auto it = elemental_data.find(name);
  if (it == elemental_data.end()) {
    AKANTU_EXCEPTION("Data named " << name
                                   << " not registered for type: " << elem_type
                                   << " - ghost_type:" << ghost_type << "!");
  }
  return aka::as_type<ElementTypeMapArray<T>>(*(it->second))(elem_type,
                                                             ghost_type);
}

template <typename T>
Array<T> & MeshData::getElementalDataArray(const ID & name,
                                           const ElementType & elem_type,
                                           const GhostType & ghost_type) {
  auto it = elemental_data.find(name);
  if (it == elemental_data.end()) {
    AKANTU_EXCEPTION("Data named " << name
                                   << " not registered for type: " << elem_type
                                   << " - ghost_type:" << ghost_type << "!");
  }
  return aka::as_type<ElementTypeMapArray<T>>(*(it->second.get()))(elem_type,
                                                                   ghost_type);
}

/* -------------------------------------------------------------------------- */
// Get an elemental data array, if it does not exist: allocate it
template <typename T>
Array<T> & MeshData::getElementalDataArrayAlloc(const ID & name,
                                                const ElementType & elem_type,
                                                const GhostType & ghost_type,
                                                UInt nb_component) {
  auto it = elemental_data.find(name);
  ElementTypeMapArray<T> * dataset;
  if (it == elemental_data.end()) {
    dataset = &allocElementalData<T>(name);
  } else {
    dataset = dynamic_cast<ElementTypeMapArray<T> *>(it->second.get());
  }
  AKANTU_DEBUG_ASSERT(
      getTypeCode<T>() ==
          typecode_map.at(MeshDataType::_elemental).find(name)->second,
      "Function getElementalDataArrayAlloc called with the wrong type!");
  if (!(dataset->exists(elem_type, ghost_type))) {
    dataset->alloc(0, nb_component, elem_type, ghost_type);
  }
  return (*dataset)(elem_type, ghost_type);
}

/* -------------------------------------------------------------------------- */
#define AKANTU_MESH_DATA_CASE_MACRO(r, name, elem)                             \
  case MeshDataTypeCode::BOOST_PP_TUPLE_ELEM(2, 0, elem): {                    \
    nb_comp = getNbComponentTemplated<BOOST_PP_TUPLE_ELEM(2, 1, elem)>(        \
        name, el_type, ghost_type);                                            \
    break;                                                                     \
  }

inline UInt MeshData::getNbComponent(const ID & name,
                                     const ElementType & el_type,
                                     const GhostType & ghost_type) const {
  auto it = typecode_map.at(MeshDataType::_elemental).find(name);
  UInt nb_comp(0);
  if (it == typecode_map.at(MeshDataType::_elemental).end()) {
    AKANTU_EXCEPTION("Could not determine the type held in dataset "
                     << name << " for type: " << el_type
                     << " - ghost_type:" << ghost_type << ".");
  }
  MeshDataTypeCode type = it->second;
  switch (type) {
    BOOST_PP_SEQ_FOR_EACH(AKANTU_MESH_DATA_CASE_MACRO, name,
                          AKANTU_MESH_DATA_TYPES)
  default:
    AKANTU_ERROR(
        "Could not call the correct instance of getNbComponentTemplated.");
    break;
  }
  return nb_comp;
}
#undef AKANTU_MESH_DATA_CASE_MACRO

/* -------------------------------------------------------------------------- */
template <typename T>
inline UInt
MeshData::getNbComponentTemplated(const ID & name, const ElementType & el_type,
                                  const GhostType & ghost_type) const {
  return getElementalDataArray<T>(name, el_type, ghost_type).getNbComponent();
}

/* -------------------------------------------------------------------------- */
inline UInt MeshData::getNbComponent(const ID & name) const {
  auto it = nodal_data.find(name);
  if (it == nodal_data.end()) {
    AKANTU_EXCEPTION("No nodal dataset registered with the name" << name
                                                                 << ".");
  }

  return it->second->getNbComponent();
}

/* -------------------------------------------------------------------------- */
// get the names of the data stored in elemental_data
#define AKANTU_MESH_DATA_CASE_MACRO(r, name, elem)                             \
  case MeshDataTypeCode::BOOST_PP_TUPLE_ELEM(2, 0, elem): {                    \
    ElementTypeMapArray<BOOST_PP_TUPLE_ELEM(2, 1, elem)> * dataset;            \
    dataset =                                                                  \
        dynamic_cast<ElementTypeMapArray<BOOST_PP_TUPLE_ELEM(2, 1, elem)> *>(  \
            it->second.get());                                                 \
    exists = dataset->exists(el_type, ghost_type);                             \
    break;                                                                     \
  }

inline auto MeshData::getTagNames(const ElementType & el_type,
                                  const GhostType & ghost_type) const {
  std::vector<std::string> tags;
  bool exists(false);

  auto it = elemental_data.begin();
  auto it_end = elemental_data.end();
  for (; it != it_end; ++it) {
    MeshDataTypeCode type = getTypeCode(it->first);
    switch (type) {
      BOOST_PP_SEQ_FOR_EACH(AKANTU_MESH_DATA_CASE_MACRO, ,
                            AKANTU_MESH_DATA_TYPES)
    default:
      AKANTU_ERROR("Could not determine the proper type to (dynamic-)cast.");
      break;
    }
    if (exists) {
      tags.push_back(it->first);
    }
  }

  return tags;
}
#undef AKANTU_MESH_DATA_CASE_MACRO

/* -------------------------------------------------------------------------- */
inline auto MeshData::getTagNames() const {
  std::vector<std::string> tags;
  for (auto && data : nodal_data) {
    tags.push_back(std::get<0>(data));
  }
  return tags;
}

/* -------------------------------------------------------------------------- */
} // namespace akantu

#endif /* __AKANTU_MESH_DATA_TMPL_HH__ */

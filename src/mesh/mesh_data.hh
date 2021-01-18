/**
 * @file   mesh_data.hh
 *
 * @author Dana Christen <dana.christen@gmail.com>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri May 03 2013
 * @date last modification: Mon Dec 18 2017
 *
 * @brief  Stores generic data loaded from the mesh file
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

#ifndef AKANTU_MESH_DATA_HH_
#define AKANTU_MESH_DATA_HH_

/* -------------------------------------------------------------------------- */
#include "aka_memory.hh"
#include "element_type_map.hh"
#include <map>
#include <string>
/* -------------------------------------------------------------------------- */

namespace akantu {

#define AKANTU_MESH_DATA_TYPES                                                 \
  ((_int, Int))((_uint, UInt))((_real, Real))((_bool, bool))(                  \
      (_element, Element))((_std_string, std::string))(                        \
      (_std_vector_element, std::vector<Element>))

#define AKANTU_MESH_DATA_TUPLE_FIRST_ELEM(s, data, elem)                       \
  BOOST_PP_TUPLE_ELEM(2, 0, elem)
enum class MeshDataTypeCode : int {
  BOOST_PP_SEQ_ENUM(BOOST_PP_SEQ_TRANSFORM(AKANTU_MESH_DATA_TUPLE_FIRST_ELEM, ,
                                           AKANTU_MESH_DATA_TYPES)),
  _unknown
};

enum class MeshDataType {
  _nodal,
  _elemental,
};

class MeshData {
  /* ------------------------------------------------------------------------ */
  /* Typedefs                                                                 */
  /* ------------------------------------------------------------------------ */
private:
  using TypeCode = MeshDataTypeCode;
  using ElementalDataMap =
      std::map<std::string, std::unique_ptr<ElementTypeMapBase>>;
  using NodalDataMap = std::map<std::string, std::unique_ptr<ArrayBase>>;
  using TypeCodeMap = std::map<std::string, TypeCode>;

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  MeshData(const ID & id = "mesh_data", const ID & parent_id = "",
           const MemoryID & mem_id = 0);

  /* ------------------------------------------------------------------------ */
  /* Methods and accessors                                                    */
  /* ------------------------------------------------------------------------ */
public:
  /// tells if the given array exists
  template <typename T>
  bool hasData(const ID & data_name, ElementType elem_type,
               GhostType ghost_type = _not_ghost) const;

  /// tells if the given data exists
  bool hasData(const ID & data_name,
               MeshDataType type = MeshDataType::_elemental) const;

  bool hasData(MeshDataType type = MeshDataType::_elemental) const;

  /// get the names of the data stored in elemental_data
  inline auto getTagNames(ElementType type,
                          GhostType ghost_type = _not_ghost) const;

  /// get the names of the data stored in elemental_data
  inline auto getTagNames() const;

  /// get the type of the data stored in elemental_data
  template <typename T> TypeCode getTypeCode() const;
  inline TypeCode
  getTypeCode(const ID & name,
              MeshDataType type = MeshDataType::_elemental) const;

  /// Get an existing elemental data array
  template <typename T>
  const Array<T> &
  getElementalDataArray(const ID & data_name, ElementType elem_type,
                        GhostType ghost_type = _not_ghost) const;
  template <typename T>
  Array<T> & getElementalDataArray(const ID & data_name,
                                   ElementType elem_type,
                                   GhostType ghost_type = _not_ghost);

  /// Get an elemental data array, if it does not exist: allocate it
  template <typename T>
  Array<T> & getElementalDataArrayAlloc(
      const ID & data_name, ElementType elem_type,
      GhostType ghost_type = _not_ghost, UInt nb_component = 1);

  template <typename T>
  inline UInt getNbComponentTemplated(const ID & name,
                                      ElementType el_type,
                                      GhostType ghost_type) const;
  inline UInt getNbComponent(const ID & name, ElementType el_type,
                             GhostType ghost_type = _not_ghost) const;

  inline UInt getNbComponent(const ID & name) const;

  /// Get an existing elemental data
  template <typename T>
  const ElementTypeMapArray<T> & getElementalData(const ID & name) const;
  template <typename T>
  ElementTypeMapArray<T> & getElementalData(const ID & name);

  template <typename T>
  Array<T> & getNodalData(const ID & name, UInt nb_components = 1);
  template <typename T> const Array<T> & getNodalData(const ID & name) const;

private:
  ///  Register new elemental data (and alloc data) with check if the name is
  ///  new
  template <typename T>
  ElementTypeMapArray<T> & registerElementalData(const ID & name);
  inline void registerElementalData(const ID & name, TypeCode type);

  ///  Register new nodal data (and alloc data) with check if the name is
  ///  new
  template <typename T>
  Array<T> & registerNodalData(const ID & name, UInt nb_components = 1);
  inline void registerNodalData(const ID & name, UInt nb_components,
                                TypeCode type);

  ///  Register new elemental data (add alloc data)
  template <typename T>
  ElementTypeMapArray<T> & allocElementalData(const ID & name);

  ///  Register new nodal data (add alloc data)
  template <typename T>
  Array<T> & allocNodalData(const ID & name, UInt nb_components);

  friend class SlaveNodeInfoPerProc;

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  ID _id;
  UInt _memory_id{0};

  /// Map when elemental data is stored as ElementTypeMap
  ElementalDataMap elemental_data;

  /// Map when elemental data is stored as ElementTypeMap
  NodalDataMap nodal_data;

  /// Map when elementalType of the data stored in elemental_data
  std::map<MeshDataType, TypeCodeMap> typecode_map{
      {MeshDataType::_elemental, {}}, {MeshDataType::_nodal, {}}};
};

} // namespace akantu

#include "mesh_data_tmpl.hh"
#undef AKANTU_MESH_DATA_TUPLE_FIRST_ELEM

#endif /* AKANTU_MESH_DATA_HH_ */

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

#ifndef __AKANTU_MESH_DATA_HH__
#define __AKANTU_MESH_DATA_HH__

/* -------------------------------------------------------------------------- */
#include "aka_memory.hh"
#include "element_type_map.hh"
#include <map>
#include <string>
/* -------------------------------------------------------------------------- */

namespace akantu {

#define AKANTU_MESH_DATA_TYPES                                                 \
  ((_tc_int, Int))((_tc_uint, UInt))((_tc_real, Real))(                        \
      (_tc_element, Element))((_tc_std_string, std::string))(                  \
      (_tc_std_vector_element, std::vector<Element>))

#define AKANTU_MESH_DATA_TUPLE_FIRST_ELEM(s, data, elem)                       \
  BOOST_PP_TUPLE_ELEM(2, 0, elem)
enum MeshDataTypeCode : int {
  BOOST_PP_SEQ_ENUM(BOOST_PP_SEQ_TRANSFORM(AKANTU_MESH_DATA_TUPLE_FIRST_ELEM, ,
                                           AKANTU_MESH_DATA_TYPES)),
  _tc_unknown
};

class MeshData : public Memory {

  /* ------------------------------------------------------------------------ */
  /* Typedefs                                                                 */
  /* ------------------------------------------------------------------------ */
private:
  using TypeCode = MeshDataTypeCode;
  using ElementalDataMap = std::map<std::string, ElementTypeMapBase *>;
  using StringVector = std::vector<std::string>;
  using TypeCodeMap = std::map<std::string, TypeCode>;

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  MeshData(const ID & id = "mesh_data", const ID & parent_id = "",
           const MemoryID & memory_id = 0);
  ~MeshData() override;

  /* ------------------------------------------------------------------------ */
  /* Methods and accessors                                                    */
  /* ------------------------------------------------------------------------ */
public:
  ///  Register new elemental data (and alloc data) with check if the name is
  ///  new
  template <typename T> void registerElementalData(const ID & name);
  inline void registerElementalData(const ID & name, TypeCode type);

  /// tells if the given array exists
  template <typename T>
  bool hasDataArray(const ID & data_name, const ElementType & el_type,
                    const GhostType & ghost_type = _not_ghost) const;

  /// tells if the given data exists
  bool hasData(const ID & data_name) const;

  /// Get an existing elemental data array
  template <typename T>
  const Array<T> &
  getElementalDataArray(const ID & data_name, const ElementType & el_type,
                        const GhostType & ghost_type = _not_ghost) const;
  template <typename T>
  Array<T> & getElementalDataArray(const ID & data_name,
                                   const ElementType & el_type,
                                   const GhostType & ghost_type = _not_ghost);

  /// Get an elemental data array, if it does not exist: allocate it
  template <typename T>
  Array<T> &
  getElementalDataArrayAlloc(const ID & data_name, const ElementType & el_type,
                             const GhostType & ghost_type = _not_ghost,
                             UInt nb_component = 1);

  /// get the names of the data stored in elemental_data
  inline void getTagNames(StringVector & tags, const ElementType & type,
                          const GhostType & ghost_type = _not_ghost) const;

  /// get the type of the data stored in elemental_data
  template <typename T> TypeCode getTypeCode() const;
  inline TypeCode getTypeCode(const ID & name) const;

  template <typename T>
  inline UInt getNbComponentTemplated(const ID & name,
                                      const ElementType & el_type,
                                      const GhostType & ghost_type) const;
  inline UInt getNbComponent(const ID & name, const ElementType & el_type,
                             const GhostType & ghost_type = _not_ghost) const;

  /// Get an existing elemental data
  template <typename T>
  const ElementTypeMapArray<T> & getElementalData(const ID & name) const;
  template <typename T>
  ElementTypeMapArray<T> & getElementalData(const ID & name);

private:
  ///  Register new elemental data (add alloc data)
  template <typename T>
  ElementTypeMapArray<T> * allocElementalData(const ID & name);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// Map when elemental data is stored as ElementTypeMap
  ElementalDataMap elemental_data;
  /// Map when elementalType of the data stored in elemental_data
  TypeCodeMap typecode_map;
};

} // akantu

#include "mesh_data_tmpl.hh"
#undef AKANTU_MESH_DATA_TUPLE_FIRST_ELEM

#endif /* __AKANTU_MESH_DATA_HH__ */

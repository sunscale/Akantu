/**
 * @file   mesh_data.hh
 *
 * @author Dana Christen <dana.christen@gmail.com>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri May 03 2013
 * @date last modification: Thu Jun 05 2014
 *
 * @brief  Stores generic data loaded from the mesh file
 *
 * @section LICENSE
 *
 * Copyright (©) 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

#ifndef __AKANTU_MESH_DATA_HH__
#define __AKANTU_MESH_DATA_HH__

/* -------------------------------------------------------------------------- */
#include "element_type_map.hh"
#include "aka_memory.hh"
#include <map>
#include <string>
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

#define AKANTU_MESH_DATA_TYPES     \
          ((_tc_int,  Int))        \
          ((_tc_uint, UInt))       \
          ((_tc_real, Real))       \
          ((_tc_element, Element)) \
          ((_tc_std_string, std::string))         \
          ((_tc_std_vector_element, std::vector<Element>)) \

#define AKANTU_MESH_DATA_TUPLE_FIRST_ELEM(s, data, elem) BOOST_PP_TUPLE_ELEM(2, 0, elem)
enum MeshDataTypeCode {
  BOOST_PP_SEQ_ENUM(BOOST_PP_SEQ_TRANSFORM(AKANTU_MESH_DATA_TUPLE_FIRST_ELEM, , AKANTU_MESH_DATA_TYPES)),
  _tc_unknown
};

class MeshData : public Memory {

  /* ------------------------------------------------------------------------ */
  /* Typedefs                                                                 */
  /* ------------------------------------------------------------------------ */
private:
  typedef MeshDataTypeCode TypeCode;
  typedef std::map<std::string, ElementTypeMapBase *> ElementalDataMap;
  typedef std::vector<std::string> StringVector;
  typedef std::map<std::string, TypeCode> TypeCodeMap;

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
   MeshData(const ID & id = "mesh_data", const ID & parent_id = "", const MemoryID & memory_id = 0);
  ~MeshData();

  /* ------------------------------------------------------------------------ */
  /* Methods and accessors                                                    */
  /* ------------------------------------------------------------------------ */
public:
  template<typename T>
  void registerElementalData(const std::string & name);
  inline void registerElementalData(const std::string & name, TypeCode type);

  template<typename T>
  const Array<T> & getElementalDataArray(const std::string & data_name,
                                         const ElementType & el_type,
                                         const GhostType & ghost_type = _not_ghost) const;

  template<typename T>
  Array<T> & getElementalDataArray(const std::string & data_name,
                                   const ElementType & el_type,
                                   const GhostType & ghost_type = _not_ghost);

  template<typename T>
  Array<T> & getElementalDataArrayAlloc(const std::string & data_name,
                                        const ElementType & el_type,
                                        const GhostType & ghost_type = _not_ghost,
                                        UInt nb_component = 1);

  inline void getTagNames(StringVector & tags, const ElementType & type, const GhostType & ghost_type = _not_ghost) const;

  inline TypeCode getTypeCode(const std::string name) const;

  template<typename T>
  inline UInt getNbComponentTemplated(const std::string, const ElementType & type, const GhostType & ghost_type) const;
  inline UInt getNbComponent(const std::string name, const ElementType & el_type, const GhostType & ghost_type = _not_ghost) const;

  template<typename T>
  const ElementTypeMapArray<T> & getElementalData(const std::string & name) const;

  template<typename T>
  ElementTypeMapArray<T> & getElementalData(const std::string & name);

  template<typename T>
  TypeCode getTypeCode() const;

private:
  template<typename T>
  ElementTypeMapArray<T> * allocElementalData(const std::string & name);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  ElementalDataMap elemental_data;
  TypeCodeMap typecode_map;
};

#include "mesh_data_tmpl.hh"
#undef AKANTU_MESH_DATA_TUPLE_FIRST_ELEM

__END_AKANTU__

#endif /* __AKANTU_MESH_DATA_HH__ */



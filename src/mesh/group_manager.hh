/**
 * @file   group_manager.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Dana Christen <dana.christen@gmail.com>
 * @author David Simon Kammer <david.kammer@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Wed Nov 13 2013
 * @date last modification: Wed Feb 07 2018
 *
 * @brief  Stores information relevent to the notion of element and nodes
 * groups.
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

#ifndef AKANTU_GROUP_MANAGER_HH_
#define AKANTU_GROUP_MANAGER_HH_

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "aka_iterators.hh"
#include "element_type_map.hh"
/* -------------------------------------------------------------------------- */
#include <set>
/* -------------------------------------------------------------------------- */

namespace akantu {
class ElementGroup;
class NodeGroup;
class Mesh;
class Element;
class ElementSynchronizer;
template <bool> class CommunicationBufferTemplated;
namespace dumpers {
  class Field;
}
} // namespace akantu

namespace akantu {

/* -------------------------------------------------------------------------- */
class GroupManager {
  /* ------------------------------------------------------------------------ */
  /* Typedefs                                                                 */
  /* ------------------------------------------------------------------------ */
private:
  using ElementGroups = std::map<std::string, std::unique_ptr<ElementGroup>>;
  using NodeGroups = std::map<std::string, std::unique_ptr<NodeGroup>>;

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  GroupManager(Mesh & mesh, const ID & id = "group_manager",
               const MemoryID & mem_id = 0);
  virtual ~GroupManager();

  /* ------------------------------------------------------------------------ */
  /* Groups iterators                                                         */
  /* ------------------------------------------------------------------------ */
public:
  using node_group_iterator = NodeGroups::iterator;
  using element_group_iterator = ElementGroups::iterator;

  using const_node_group_iterator = NodeGroups::const_iterator;
  using const_element_group_iterator = ElementGroups::const_iterator;

#define AKANTU_GROUP_MANAGER_DEFINE_ITERATOR_FUNCTION(group_type, function,        \
                                                      param_in, param_out)         \
  [[deprecated(                                                                    \
      "use iterate(Element|Node)Groups or "                                        \
      "(element|node)GroupExists")]] inline BOOST_PP_CAT(BOOST_PP_CAT(const_,      \
                                                                      group_type), \
                                                         _iterator)                \
      BOOST_PP_CAT(BOOST_PP_CAT(group_type, _), function)(param_in) const {        \
    return BOOST_PP_CAT(group_type, s).function(param_out);                        \
  };                                                                               \
                                                                                   \
  [[deprecated("use iterate(Element|Node)Groups or "                               \
               "(element|node)GroupExists")]] inline BOOST_PP_CAT(group_type,      \
                                                                  _iterator)       \
      BOOST_PP_CAT(BOOST_PP_CAT(group_type, _), function)(param_in) {              \
    return BOOST_PP_CAT(group_type, s).function(param_out);                        \
  }

#define AKANTU_GROUP_MANAGER_DEFINE_ITERATOR_FUNCTION_NP(group_type, function) \
  AKANTU_GROUP_MANAGER_DEFINE_ITERATOR_FUNCTION(                               \
      group_type, function, BOOST_PP_EMPTY(), BOOST_PP_EMPTY())

  AKANTU_GROUP_MANAGER_DEFINE_ITERATOR_FUNCTION_NP(node_group, begin);
  AKANTU_GROUP_MANAGER_DEFINE_ITERATOR_FUNCTION_NP(node_group, end);
  AKANTU_GROUP_MANAGER_DEFINE_ITERATOR_FUNCTION_NP(element_group, begin);
  AKANTU_GROUP_MANAGER_DEFINE_ITERATOR_FUNCTION_NP(element_group, end);
  AKANTU_GROUP_MANAGER_DEFINE_ITERATOR_FUNCTION(element_group, find,
                                                const std::string & name, name);
  AKANTU_GROUP_MANAGER_DEFINE_ITERATOR_FUNCTION(node_group, find,
                                                const std::string & name, name);

public:
  decltype(auto) iterateNodeGroups() {
    return make_dereference_adaptor(make_values_adaptor(node_groups));
  }
  decltype(auto) iterateNodeGroups() const {
    return make_dereference_adaptor(make_values_adaptor(node_groups));
  }

  decltype(auto) iterateElementGroups() {
    return make_dereference_adaptor(make_values_adaptor(element_groups));
  }
  decltype(auto) iterateElementGroups() const {
    return make_dereference_adaptor(make_values_adaptor(element_groups));
  }

  /* ------------------------------------------------------------------------ */
  /* Clustering filter                                                        */
  /* ------------------------------------------------------------------------ */
public:
  class ClusteringFilter {
  public:
    virtual bool operator()(const Element & /*unused*/) const { return true; }
  };

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// create an empty node group
  NodeGroup & createNodeGroup(const std::string & group_name,
                              bool replace_group = false);

  /// create an element group and the associated node group
  ElementGroup & createElementGroup(const std::string & group_name,
                                    UInt dimension = _all_dimensions,
                                    bool replace_group = false);

  /* ------------------------------------------------------------------------ */
  /// renames an element group
  void renameElementGroup(const std::string & name,
                          const std::string & new_name);

  /// renames a node group
  void renameNodeGroup(const std::string & name, const std::string & new_name);

  /// copy an existing element group
  void copyElementGroup(const std::string & name, const std::string & new_name);

  /// copy an existing node group
  void copyNodeGroup(const std::string & name, const std::string & new_name);

  /* ------------------------------------------------------------------------ */

  /// create a node group from another node group but filtered
  template <typename T>
  NodeGroup & createFilteredNodeGroup(const std::string & group_name,
                                      const NodeGroup & node_group, T & filter);

  /// create an element group from another element group but filtered
  template <typename T>
  ElementGroup &
  createFilteredElementGroup(const std::string & group_name, UInt dimension,
                             const NodeGroup & node_group, T & filter);

  /// destroy a node group
  void destroyNodeGroup(const std::string & group_name);

  /// destroy an element group and the associated node group
  void destroyElementGroup(const std::string & group_name,
                           bool destroy_node_group = false);

  // /// destroy all element groups and the associated node groups
  // void destroyAllElementGroups(bool destroy_node_groups = false);

  /// create a element group using an existing node group
  ElementGroup & createElementGroup(const std::string & group_name,
                                    UInt dimension, NodeGroup & node_group);

  /// create groups based on values stored in a given mesh data
  template <typename T>
  void createGroupsFromMeshData(const std::string & dataset_name);

  /// create boundaries group by a clustering algorithm \todo extend to parallel
  UInt createBoundaryGroupFromGeometry();

  /// create element clusters for a given dimension
  UInt createClusters(UInt element_dimension, Mesh & mesh_facets,
                      const std::string & cluster_name_prefix = "cluster",
                      const ClusteringFilter & filter = ClusteringFilter());

  /// create element clusters for a given dimension
  UInt createClusters(UInt element_dimension,
                      const std::string & cluster_name_prefix = "cluster",
                      const ClusteringFilter & filter = ClusteringFilter());

private:
  /// create element clusters for a given dimension
  UInt createClusters(UInt element_dimension,
                      const std::string & cluster_name_prefix,
                      const ClusteringFilter & filter, Mesh & mesh_facets);

public:
  /// Create an ElementGroup based on a NodeGroup
  void createElementGroupFromNodeGroup(const std::string & name,
                                       const std::string & node_group,
                                       UInt dimension = _all_dimensions);

  virtual void printself(std::ostream & stream, int indent = 0) const;

  /// this function insure that the group names are present on all processors
  /// /!\ it is a SMP call
  void synchronizeGroupNames();

  /// register an elemental field to the given group name (overloading for
  /// ElementalPartionField)
  template <typename T, template <bool> class dump_type>
  std::shared_ptr<dumpers::Field> createElementalField(
      const ElementTypeMapArray<T> & field, const std::string & group_name,
      UInt spatial_dimension, ElementKind kind,
      ElementTypeMap<UInt> nb_data_per_elem = ElementTypeMap<UInt>());

  /// register an elemental field to the given group name (overloading for
  /// ElementalField)
  template <typename T, template <class> class ret_type,
            template <class, template <class> class, bool> class dump_type>
  std::shared_ptr<dumpers::Field> createElementalField(
      const ElementTypeMapArray<T> & field, const std::string & group_name,
      UInt spatial_dimension, ElementKind kind,
      ElementTypeMap<UInt> nb_data_per_elem = ElementTypeMap<UInt>());

  /// register an elemental field to the given group name (overloading for
  /// MaterialInternalField)
  template <typename T,
            /// type of InternalMaterialField
            template <typename, bool filtered> class dump_type>
  std::shared_ptr<dumpers::Field>
  createElementalField(const ElementTypeMapArray<T> & field,
                       const std::string & group_name, UInt spatial_dimension,
                       ElementKind kind,
                       ElementTypeMap<UInt> nb_data_per_elem);

  template <typename type, bool flag, template <class, bool> class ftype>
  std::shared_ptr<dumpers::Field>
  createNodalField(const ftype<type, flag> * field,
                   const std::string & group_name, UInt padding_size = 0);

  template <typename type, bool flag, template <class, bool> class ftype>
  std::shared_ptr<dumpers::Field>
  createStridedNodalField(const ftype<type, flag> * field,
                          const std::string & group_name, UInt size,
                          UInt stride, UInt padding_size);

protected:
  /// fill a buffer with all the group names
  void fillBufferWithGroupNames(
      CommunicationBufferTemplated<false> & comm_buffer) const;

  /// take a buffer and create the missing groups localy
  void checkAndAddGroups(CommunicationBufferTemplated<false> & buffer);

  /// register an elemental field to the given group name
  template <class dump_type, typename field_type>
  inline std::shared_ptr<dumpers::Field>
  createElementalField(const field_type & field, const std::string & group_name,
                       UInt spatial_dimension, ElementKind kind,
                       const ElementTypeMap<UInt> & nb_data_per_elem);

  /// register an elemental field to the given group name
  template <class dump_type, typename field_type>
  inline std::shared_ptr<dumpers::Field>
  createElementalFilteredField(const field_type & field,
                               const std::string & group_name,
                               UInt spatial_dimension, ElementKind kind,
                               ElementTypeMap<UInt> nb_data_per_elem);

  /* ------------------------------------------------------------------------ */
  /* Accessor                                                                 */
  /* ------------------------------------------------------------------------ */
public:
  // AKANTU_GET_MACRO(ElementGroups, element_groups, const ElementGroups &);

  const ElementGroup & getElementGroup(const std::string & name) const;
  const NodeGroup & getNodeGroup(const std::string & name) const;

  ElementGroup & getElementGroup(const std::string & name);
  NodeGroup & getNodeGroup(const std::string & name);

  UInt getNbElementGroups(UInt dimension = _all_dimensions) const;
  UInt getNbNodeGroups() { return node_groups.size(); };

  bool elementGroupExists(const std::string & name) {
    return element_groups.find(name) != element_groups.end();
  }

  bool nodeGroupExists(const std::string & name) {
    return node_groups.find(name) != node_groups.end();
  }

private:
  template <typename GroupsType>
  void renameGroup(GroupsType & groups, const std::string & name,
                   const std::string & new_name);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// id to create element and node groups
  ID id;
  /// memory_id to create element and node groups
  MemoryID memory_id;

  /// list of the node groups managed
  NodeGroups node_groups;

  /// list of the element groups managed
  ElementGroups element_groups;

  /// Mesh to which the element belongs
  Mesh & mesh;
};

/// standard output stream operator
inline std::ostream & operator<<(std::ostream & stream,
                                 const GroupManager & _this) {
  _this.printself(stream);
  return stream;
}

} // namespace akantu

#endif /* AKANTU_GROUP_MANAGER_HH_ */

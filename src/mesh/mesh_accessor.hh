/**
 * @file   mesh_accessor.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Tue Jun 30 2015
 * @date last modification: Tue Sep 19 2017
 *
 * @brief  this class allow to access some private member of mesh it is used for
 * IO for examples
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
#include "mesh.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_MESH_ACCESSOR_HH_
#define AKANTU_MESH_ACCESSOR_HH_

namespace akantu {
class NodeSynchronizer;
class ElementSynchronizer;
class MeshGlobalDataUpdater;
} // namespace akantu

namespace akantu {

class MeshAccessor {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  explicit MeshAccessor(Mesh & mesh) : _mesh(mesh) {}
  virtual ~MeshAccessor() = default;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// get the global number of nodes
  inline UInt getNbGlobalNodes() const { return this->_mesh.nb_global_nodes; }

  /// set the global number of nodes
  inline void setNbGlobalNodes(UInt nb_global_nodes) {
    this->_mesh.nb_global_nodes = nb_global_nodes;
  }

  /// set the mesh as being distributed
  inline void setDistributed() { this->_mesh.is_distributed = true; }

  /// get a pointer to the nodes_global_ids Array<UInt> and create it if
  /// necessary
  inline auto & getNodesGlobalIds() {
    return this->_mesh.getNodesGlobalIdsPointer();
  }

  /// get a pointer to the nodes_type Array<Int> and create it if necessary
  inline auto & getNodesFlags() { return this->_mesh.getNodesFlags(); }

  /// get a pointer to the nodes_type Array<Int> and create it if necessary
  inline void setNodePrank(UInt node, Int prank) {
    this->_mesh.nodes_prank[node] = prank;
  }

  /// get a pointer to the coordinates Array
  inline auto & getNodes() { return this->_mesh.getNodesPointer(); }

  /// get a pointer to the coordinates Array
  inline auto getNodesSharedPtr() { return this->_mesh.nodes; }

  /// get the connectivities
  inline auto & getConnectivities() { return this->_mesh.connectivities; }

  /// get the connectivity Array for the given type and create it
  /// if necessary
  inline auto & getConnectivity(ElementType type,
                                GhostType ghost_type = _not_ghost) {
    return this->_mesh.getConnectivityPointer(type, ghost_type);
  }

  /// get the connectivity for the given element
  inline decltype(auto) getConnectivity(const Element & element) {
    return this->_mesh.getConnectivityNC(element);
  }

  /// get the ghost element counter
  inline auto & getGhostsCounters(ElementType type,
                                  GhostType ghost_type = _ghost) {
    return this->_mesh.getGhostsCounters(type, ghost_type);
  }

  /// get the element_to_subelement Array for the given type and
  /// create it if necessary
  inline auto &
  getElementToSubelement(ElementType type,
                         GhostType ghost_type = _not_ghost) {
    return this->_mesh.getElementToSubelementPointer(type, ghost_type);
  }

  inline decltype(auto)
  getElementToSubelementNC(const ElementType & type,
                         const GhostType & ghost_type = _not_ghost) {
    return this->_mesh.getElementToSubelementNC(type, ghost_type);
  }

  /// get the subelement_to_element Array for the given type and
  /// create it if necessary
  inline auto &
  getSubelementToElement(ElementType type,
                         GhostType ghost_type = _not_ghost) {
    return this->_mesh.getSubelementToElementPointer(type, ghost_type);
  }

  inline decltype(auto)
  getSubelementToElementNC(const ElementType & type,
                           const GhostType & ghost_type = _not_ghost) {
    return this->_mesh.getSubelementToElementNC(type, ghost_type);
  }

  /// get  the element_to_subelement, creates it if necessary
  inline decltype(auto) getElementToSubelement() {
    return this->_mesh.getElementToSubelementNC();
  }

  /// get subelement_to_element, creates it if necessary
  inline decltype(auto) getSubelementToElement() {
    return this->_mesh.getSubelementToElementNC();
  }

  /// get a pointer to the element_to_subelement Array for element and
  /// create it if necessary
  inline decltype(auto) getElementToSubelement(const Element & element) {
    return this->_mesh.getElementToSubelementNC(element);
  }

  /// get a pointer to the subelement_to_element Array for the given element and
  /// create it if necessary
  inline decltype(auto) getSubelementToElement(const Element & element) {
    return this->_mesh.getSubelementToElementNC(element);
  }

  template <typename T>
  inline auto &
  getData(const std::string & data_name, ElementType el_type,
          GhostType ghost_type = _not_ghost, UInt nb_component = 1,
          bool size_to_nb_element = true, bool resize_with_parent = false) {
    return this->_mesh.getDataPointer<T>(data_name, el_type, ghost_type,
                                         nb_component, size_to_nb_element,
                                         resize_with_parent);
  }

  /// get the node synchonizer
  auto & getNodeSynchronizer() { return *this->_mesh.node_synchronizer; }

  /// get the element synchonizer
  auto & getElementSynchronizer() { return *this->_mesh.element_synchronizer; }

  decltype(auto) updateGlobalData(NewNodesEvent & nodes_event,
                                  NewElementsEvent & elements_event) {
    return this->_mesh.updateGlobalData(nodes_event, elements_event);
  }

  void registerGlobalDataUpdater(
      std::unique_ptr<MeshGlobalDataUpdater> && global_data_updater) {
    this->_mesh.registerGlobalDataUpdater(
        std::forward<std::unique_ptr<MeshGlobalDataUpdater>>(
            global_data_updater));
  }

  /* ------------------------------------------------------------------------ */
  void makeReady() { this->_mesh.makeReady(); }

  /* ------------------------------------------------------------------------ */
  void addPeriodicSlave(UInt slave, UInt master) {
    this->_mesh.addPeriodicSlave(slave, master);
  }

  void markMeshPeriodic() {
    for (UInt s : arange(this->_mesh.spatial_dimension)) {
      this->_mesh.is_periodic |= 1 << s;
    }
  }

  void wipePeriodicInfo() { this->_mesh.wipePeriodicInfo(); }

private:
  Mesh & _mesh;
};

} // namespace akantu

#endif /* AKANTU_MESH_ACCESSOR_HH_ */

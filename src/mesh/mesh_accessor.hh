/**
 * @file   mesh_accessor.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Tue Jun 30 2015
 *
 * @brief  this class allow to access some private member of mesh it is used for
 * IO for examples
 *
 * @section LICENSE
 *
 * Copyright (©) 2015 EPFL (Ecole Polytechnique Fédérale de Lausanne) Laboratory
 * (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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
#include "mesh.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MESH_ACCESSOR_HH__
#define __AKANTU_MESH_ACCESSOR_HH__

namespace akantu {
class NodeSynchronizer;
class ElementSynchronizer;
}  // akantu

namespace akantu {

class MeshAccessor {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  MeshAccessor(Mesh & mesh) : _mesh(mesh){};
  virtual ~MeshAccessor(){};

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// get the global number of nodes
  inline UInt getNbGlobalNodes() const { return this->_mesh.nb_global_nodes; };

  /// set the global number of nodes
  inline void setNbGlobalNodes(UInt nb_global_nodes) {
    this->_mesh.nb_global_nodes = nb_global_nodes;
  };

  /// set the mesh as being distributed
  inline void setDistributed() { this->_mesh.is_distributed = true; }

  /// get a pointer to the nodes_global_ids Array<UInt> and create it if
  /// necessary
  inline Array<UInt> & getNodesGlobalIds() {
    return *(this->_mesh.getNodesGlobalIdsPointer());
  }

  /// get a pointer to the nodes_type Array<Int> and create it if necessary
  inline Array<NodeType> & getNodesType() {
    return *(this->_mesh.getNodesTypePointer());
  }

  /// get a pointer to the coordinates Array
  inline Array<Real> & getNodes() { return *(this->_mesh.getNodesPointer()); }

  /// get a pointer to the connectivity Array for the given type and create it
  /// if necessary
  inline Array<UInt> &
  getConnectivity(const ElementType & type,
                  const GhostType & ghost_type = _not_ghost) {
    return *(this->_mesh.getConnectivityPointer(type, ghost_type));
  }

  /// get a pointer to the element_to_subelement Array for the given type and
  /// create it if necessary
  inline Array<std::vector<Element> > &
  getElementToSubelement(const ElementType & type,
                         const GhostType & ghost_type = _not_ghost) {
    return *(this->_mesh.getElementToSubelementPointer(type, ghost_type));
  }

  /// get a pointer to the subelement_to_element Array for the given type and
  /// create it if necessary
  inline Array<Element> &
  getSubelementToElement(const ElementType & type,
                         const GhostType & ghost_type = _not_ghost) {
    return *(this->_mesh.getSubelementToElementPointer(type, ghost_type));
  }

  template <typename T>
  inline Array<T> &
  getData(const std::string & data_name, const ElementType & el_type,
          const GhostType & ghost_type = _not_ghost, UInt nb_component = 1,
          bool size_to_nb_element = true, bool resize_with_parent = false) {
    return *(this->_mesh.getDataPointer<T>(data_name, el_type, ghost_type,
                                           nb_component, size_to_nb_element,
                                           resize_with_parent));
  }

  MeshData & getMeshData() { return this->_mesh.getMeshData(); }

  /// get the node synchonizer
  NodeSynchronizer & getNodeSynchronizer();

  /// get the element synchonizer
  ElementSynchronizer & getElementSynchronizer();

private:
  Mesh & _mesh;
};

} // akantu

#endif /* __AKANTU_MESH_ACCESSOR_HH__ */

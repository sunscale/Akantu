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

#ifndef __AKANTU_MESH_ACCESSOR_HH__
#define __AKANTU_MESH_ACCESSOR_HH__

__BEGIN_AKANTU__

class MeshAccessor {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  MeshAccessor() {};
  virtual ~MeshAccessor() {};

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// get the global number of nodes
  inline UInt getNbGlobalNodes(const Mesh & mesh) const { return mesh.nb_global_nodes; };

  /// set the global number of nodes
  inline void setNbGlobalNodes(Mesh & mesh, UInt nb_global_nodes) const
  { mesh.nb_global_nodes = nb_global_nodes; };

  /// get a pointer to the nodes_global_ids Array<UInt> and create it if necessary
  inline Array<UInt> & getNodesGlobalIds(Mesh & mesh) {
    return *(mesh.getNodesGlobalIdsPointer());
  }

  /// get a pointer to the nodes_type Array<Int> and create it if necessary
  inline Array<Int> & getNodesType(Mesh & mesh) {
    return *(mesh.getNodesTypePointer());
  }

  /// get a pointer to the connectivity Array for the given type and create it if necessary
  inline Array<UInt> & getConnectivity(Mesh & mesh,
				       const ElementType & type,
				       const GhostType & ghost_type = _not_ghost) {
    return *(mesh.getConnectivityPointer(type, ghost_type));
  }

  /// get a pointer to the element_to_subelement Array for the given type and create it if necessary
  inline Array< std::vector<Element> > & getElementToSubelement(Mesh & mesh,
								const ElementType & type,
								const GhostType & ghost_type = _not_ghost) {
    return *(mesh.getElementToSubelementPointer(type, ghost_type));
  }

  /// get a pointer to the subelement_to_element Array for the given type and create it if necessary
  inline Array<Element > & getSubelementToElement(Mesh & mesh,
						  const ElementType & type,
						  const GhostType & ghost_type = _not_ghost) {
    return *(mesh.getSubelementToElementPointer(type, ghost_type));
  }

  template<typename T>
  inline Array<T> & getData(Mesh & mesh,
			    const std::string & data_name,
			    const ElementType & el_type,
			    const GhostType & ghost_type = _not_ghost,
			    UInt nb_component = 1,
			    bool size_to_nb_element = true,
			    bool resize_with_parent = false) {
    return *(mesh.getDataPointer<T>(data_name,
				    el_type,
				    ghost_type,
				    nb_component,
				    size_to_nb_element,
				    resize_with_parent));
  }
};


__END_AKANTU__


#endif /* __AKANTU_MESH_ACCESSOR_HH__ */

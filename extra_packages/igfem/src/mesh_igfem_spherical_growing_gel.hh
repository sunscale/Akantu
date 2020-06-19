/**
 * @file mesh_igfem_spherical_growing_gel.hh
 *
 * @author Clement Roux-Langlois <clement.roux@epfl.ch>
 *
 * @date creation: Mon Jul 13 2015
 *
 * @brief Computation of mesh intersection with sphere(s) and growing of these
 *        spheres. This class handle the intersectors templated for every
 * element
 *        types.
 *
 *
 * Copyright (©) 2010-2015 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

//#if 0
#ifndef __AKANTU_MESH_IGFEM_SPHERICAL_GROWING_GEL_HH__
#define __AKANTU_MESH_IGFEM_SPHERICAL_GROWING_GEL_HH__

#include "aka_common.hh"
#include "mesh_sphere_intersector.hh"
#include "mesh_utils.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
/* classes for new igfem elements mesh events */
/* -------------------------------------------------------------------------- */

class NewIGFEMElementsEvent : public NewElementsEvent {
public:
  AKANTU_GET_MACRO_NOT_CONST(OldElementsList, old_elements, Array<Element> &);
  AKANTU_GET_MACRO(OldElementsList, old_elements, const Array<Element> &);

protected:
  Array<Element> old_elements;
};

class NewIGFEMNodesEvent : public NewNodesEvent {
public:
  void setNewNodePerElem(const Array<UInt> & new_node_per_elem) {
    this->new_node_per_elem = &new_node_per_elem;
  }
  void setType(ElementType new_type) { type = new_type; }
  void setGhostType(GhostType new_type) { ghost_type = new_type; }
  AKANTU_GET_MACRO(NewNodePerElem, *new_node_per_elem, const Array<UInt> &);
  AKANTU_GET_MACRO(ElementType, type, ElementType);
  AKANTU_GET_MACRO(GhostType, ghost_type, GhostType);

protected:
  ElementType type;
  GhostType ghost_type;
  const Array<UInt> * new_node_per_elem;
};

/* -------------------------------------------------------------------------- */
template <UInt dim> class MeshIgfemSphericalGrowingGel {
// definition of the element list
#define ELEMENT_LIST (_triangle_3)(_igfem_triangle_4)(_igfem_triangle_5)

// Solution 1
/* #define INTERSECTOR_DEFINITION(type)			\
   MeshSphereIntersector<2, type> intersector##type(mesh);*/

// Solution 2
#define INSTANTIATOR(_type)                                                    \
  intersectors(_type, ghost_type) =                                            \
      new MeshSphereIntersector<dim, _type>(this->mesh)

  /* --------------------------------------------------------------------------
   */
  /* Constructor/Destructor */
  /* --------------------------------------------------------------------------
   */

public:
  /// Construct from mesh
  MeshIgfemSphericalGrowingGel(Mesh & mesh);

  /// Destructor
  ~MeshIgfemSphericalGrowingGel() {
    for (ghost_type_t::iterator gt = ghost_type_t::begin();
         gt != ghost_type_t::end(); ++gt) {
      GhostType ghost_type = *gt;

      Mesh::type_iterator it = mesh.firstType(dim, ghost_type, _ek_not_defined);
      Mesh::type_iterator end = mesh.lastType(dim, ghost_type, _ek_not_defined);

      for (; it != end; ++it) {
        delete intersectors(*it, ghost_type);
      }
    }
  }

  /* --------------------------------------------------------------------------
   */
  /* Methods */
  /* --------------------------------------------------------------------------
   */
public:
  void init();

  /// Remove the additionnal nodes
  void removeAdditionalNodes();

  /// Compute the intersection points between the mesh and the query list for
  /// all element types and send the NewNodeEvent
  void computeMeshQueryListIntersectionPoint(
      const std::list<SK::Sphere_3> & query_list);

  /// increase sphere radius and build the new intersection points between the
  /// mesh and the query list for all element types and send the NewNodeEvent
  void computeMeshQueryListIntersectionPoint(
      const std::list<SK::Sphere_3> & query_list, Real inf_fact) {
    std::list<SK::Sphere_3>::const_iterator query_it = query_list.begin(),
                                            query_end = query_list.end();
    std::list<SK::Sphere_3> sphere_list;
    for (; query_it != query_end; ++query_it) {
      SK::Sphere_3 sphere(query_it->center(),
                          query_it->squared_radius() * inf_fact * inf_fact);
      sphere_list.push_back(sphere);
    }
    computeMeshQueryListIntersectionPoint(sphere_list);
  }

  /// Build the IGFEM mesh from intersection points
  void buildIGFEMMesh();

  /// Build the IGFEM mesh from spheres
  void buildIGFEMMeshFromSpheres(const std::list<SK::Sphere_3> & query_list) {
    computeMeshQueryListIntersectionPoint(query_list);
    buildIGFEMMesh();
  }

  /// Build the IGFEM mesh from spheres with increase factor
  void buildIGFEMMeshFromSpheres(const std::list<SK::Sphere_3> & query_list,
                                 Real inf_fact) {
    computeMeshQueryListIntersectionPoint(query_list, inf_fact);
    buildIGFEMMesh();
  }

  /// set the distributed synchronizer
  void setDistributedSynchronizer(DistributedSynchronizer * dist) {
    synchronizer = dist;
    buildSegmentConnectivityToNodeType();
  }

  /// update node type
  void updateNodeType(const Array<UInt> & nodes_list,
                      const Array<UInt> & new_node_per_elem, ElementType type,
                      GhostType ghost_type);

protected:
  /// Build the unordered_map needed to assign the node type to new nodes in
  /// parallel
  void buildSegmentConnectivityToNodeType();

  /* --------------------------------------------------------------------------
   */
  /* Accessors */
  /* --------------------------------------------------------------------------
   */
public:
  AKANTU_GET_MACRO(NbStandardNodes, nb_nodes_fem, UInt);
  AKANTU_GET_MACRO(NbEnrichedNodes, nb_enriched_nodes, UInt);

  /* --------------------------------------------------------------------------
   */
  /* Class Members */
  /* --------------------------------------------------------------------------
   */
protected:
  /// Mesh used to construct the primitives
  Mesh & mesh;

  /// number of fem nodes in the initial mesh
  UInt nb_nodes_fem;

  /// number of enriched nodes before intersection
  UInt nb_enriched_nodes;

  // Solution 2
  /// map of the elements types in the mesh and the corresponding intersectors
  ElementTypeMap<MeshAbstractIntersector<SK::Sphere_3> *> intersectors;

  /// Map linking pairs of nodes to a node type. The pairs of nodes
  /// contain the connectivity of the primitive segments that are
  /// intersected.
  unordered_map<std::pair<UInt, UInt>, Int>::type segment_conn_to_node_type;

  /// Pointer to the distributed synchronizer of the model
  DistributedSynchronizer * synchronizer;
};

} // namespace akantu

#include "mesh_igfem_spherical_growing_gel_tmpl.hh"

#endif // __AKANTU_MESH_IGFEM_SPHERICAL_GROWING_GEL_HH__

//#endif //

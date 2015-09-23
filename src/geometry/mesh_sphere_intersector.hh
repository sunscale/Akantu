/**
 * @file mesh_segment_intersector.hh
 *
 * @author Clement Roux-Langlois <clement.roux@epfl.ch>
 *
 * @date creation: Wed Jun 10 2015
 *
 * @brief Computation of mesh intersection with sphere(s)
 *
 * @section LICENSE
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

#ifndef __AKANTU_MESH_SPHERE_INTERSECTOR_HH__
#define __AKANTU_MESH_SPHERE_INTERSECTOR_HH__

#include "aka_common.hh"
#include "mesh_geom_intersector.hh"

#include "mesh_geom_common.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
/* class for new igfem elements mesh events                                   */
/* -------------------------------------------------------------------------- */
#if defined(AKANTU_IGFEM)
class NewIGFEMElementsEvent : public NewElementsEvent {
public:
  AKANTU_GET_MACRO_NOT_CONST(OldElementsList, old_elements, Array<Element> &);
  AKANTU_GET_MACRO(OldElementsList, old_elements, const Array<Element> &);
protected:
  Array<Element> old_elements;
};

class NewIGFEMNodesEvent : public NewNodesEvent {
public:
  void setNewNodePerElem(ElementTypeMapUInt & new_node_per_elem) {
    this->new_node_per_elem = &new_node_per_elem;
  }
  void setType(ElementType new_type) {type = new_type;}
  AKANTU_GET_MACRO(NewNodePerElem, *new_node_per_elem, const ElementTypeMapUInt &);
  AKANTU_GET_MACRO(ElementType, type, ElementType);
protected:
  ElementType type;
  ElementTypeMapUInt * new_node_per_elem;
};


#endif


/// Here, we know what kernel we have to use
typedef Spherical SK;

template<UInt dim, ElementType type>
class MeshSphereIntersector : public MeshGeomIntersector<dim, type, Line_arc<SK>, SK::Sphere_3, SK> {
  /// Parent class type
  typedef MeshGeomIntersector<dim, type, Line_arc<SK>, SK::Sphere_3, SK> parent_type;

  /// Result of intersection function type
  typedef typename IntersectionTypeHelper<TreeTypeHelper< Triangle<K>, K>, K::Segment_3>::intersection_type result_type;

  /// Pair of intersection points and element id
  typedef std::pair<SK::Circular_arc_point_3, UInt> pair_type;

public:
  /// Construct from mesh
  explicit MeshSphereIntersector(Mesh & mesh);

  /// Destructor
  virtual ~MeshSphereIntersector();

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// get the new_node_per_elem array
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(NewNodePerElem, new_node_per_elem, UInt);

public:
  /// Construct the primitive tree object
  virtual void constructData();

  /**
   * @brief Computes the intersection of the mesh with a sphere
   *
   * @param query (sphere) to compute the intersections with the mesh
   */
  virtual void computeIntersectionQuery(const SK::Sphere_3 & query);

  /// Build the IGFEM mesh
  virtual void buildResultFromQueryList(const std::list<SK::Sphere_3> & query);

  /// Remove the additionnal nodes
  void removeAdditionnalNodes();

protected:
  /// new node per element (column 0: number of new nodes, then odd is the intersection node number and even the ID of the sintersected segment)
  ElementTypeMapUInt new_node_per_elem;

  /// number of fem nodes in the initial mesh
  const UInt nb_nodes_fem;

  /// number of primitive in an element of the template type
  UInt nb_prim_by_el;

};
 
__END_AKANTU__

#include "mesh_sphere_intersector_tmpl.hh"

#endif // __AKANTU_MESH_SPHERE_INTERSECTOR_HH__

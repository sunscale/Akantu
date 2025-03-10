/**
 * @file   mesh_geom_factory_tmpl.hh
 *
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 * @author Clement Roux <clement.roux@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Fri Feb 27 2015
 * @date last modification: Wed Jan 31 2018
 *
 * @brief  Class for constructing the CGAL primitives of a mesh
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

#ifndef __AKANTU_MESH_GEOM_FACTORY_TMPL_HH__
#define __AKANTU_MESH_GEOM_FACTORY_TMPL_HH__

/* -------------------------------------------------------------------------- */

#include "aka_common.hh"
#include "mesh_geom_common.hh"
#include "mesh_geom_factory.hh"

/* -------------------------------------------------------------------------- */

namespace akantu {

template <UInt dim, ElementType type, class Primitive, class Kernel>
MeshGeomFactory<dim, type, Primitive, Kernel>::MeshGeomFactory(Mesh & mesh)
    : MeshGeomAbstract(mesh), data_tree(NULL), primitive_list() {}

template <UInt dim, ElementType type, class Primitive, class Kernel>
MeshGeomFactory<dim, type, Primitive, Kernel>::~MeshGeomFactory() {
  delete data_tree;
}

/**
 * This function loops over the elements of `type` in the mesh and creates the
 * AABB tree of geometrical primitves (`data_tree`).
 */
template <UInt dim, ElementType type, class Primitive, class Kernel>
void MeshGeomFactory<dim, type, Primitive, Kernel>::constructData(
    GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  primitive_list.clear();
  UInt nb_nodes_per_element = mesh.getNbNodesPerElement(type);

  const Array<UInt> & connectivity = mesh.getConnectivity(type, ghost_type);
  const Array<Real> & nodes = mesh.getNodes();

  UInt el_index = 0;
  Array<UInt>::const_vector_iterator it =
      connectivity.begin(nb_nodes_per_element);
  Array<UInt>::const_vector_iterator end =
      connectivity.end(nb_nodes_per_element);

  Matrix<Real> node_coordinates(dim, nb_nodes_per_element);

  // This loop builds the list of primitives
  for (; it != end; ++it, ++el_index) {
    const Vector<UInt> & el_connectivity = *it;

    for (UInt i = 0; i < nb_nodes_per_element; i++)
      for (UInt j = 0; j < dim; j++)
        node_coordinates(j, i) = nodes(el_connectivity(i), j);

    // the unique elemental id assigned to the primitive is the
    // linearized element index over ghost type
    addPrimitive(node_coordinates, el_index);
  }

  delete data_tree;

  // This condition allows the use of the mesh geom module
  // even if types are not compatible with AABB tree algorithm
  if (TreeTypeHelper<Primitive, Kernel>::is_valid)
    data_tree = new typename TreeTypeHelper<Primitive, Kernel>::tree(
        primitive_list.begin(), primitive_list.end());

  AKANTU_DEBUG_OUT();
}

template <UInt dim, ElementType type, class Primitive, class Kernel>
void MeshGeomFactory<dim, type, Primitive, Kernel>::addPrimitive(
    const Matrix<Real> & node_coordinates, UInt id) {
  this->addPrimitive(node_coordinates, id, this->primitive_list);
}

// (2D, _triangle_3) decomposed into Triangle<cgal::Cartesian>
template <>
inline void
MeshGeomFactory<2, _triangle_3, Triangle<cgal::Cartesian>, cgal::Cartesian>::
    addPrimitive(const Matrix<Real> & node_coordinates, UInt id,
                 TreeTypeHelper<Triangle<cgal::Cartesian>,
                                cgal::Cartesian>::container_type & list) {

  TreeTypeHelper<Triangle<cgal::Cartesian>, cgal::Cartesian>::point_type a(
      node_coordinates(0, 0), node_coordinates(1, 0), 0.),
      b(node_coordinates(0, 1), node_coordinates(1, 1), 0.),
      c(node_coordinates(0, 2), node_coordinates(1, 2), 0.);

  Triangle<cgal::Cartesian> t(a, b, c);
  t.setId(id);
  list.push_back(t);
}

// (2D, _triangle_6) decomposed into Triangle<cgal::Cartesian>
template <>
inline void
MeshGeomFactory<2, _triangle_6, Triangle<cgal::Cartesian>, cgal::Cartesian>::
    addPrimitive(const Matrix<Real> & node_coordinates, UInt id,
                 TreeTypeHelper<Triangle<cgal::Cartesian>,
                                cgal::Cartesian>::container_type & list) {

  TreeTypeHelper<Triangle<cgal::Cartesian>, cgal::Cartesian>::point_type a(
      node_coordinates(0, 0), node_coordinates(1, 0), 0.),
      b(node_coordinates(0, 1), node_coordinates(1, 1), 0.),
      c(node_coordinates(0, 2), node_coordinates(1, 2), 0.);

  Triangle<cgal::Cartesian> t(a, b, c);
  t.setId(id);
  list.push_back(t);
}

// (2D, _triangle_3) decomposed into Line_arc<cgal::Spherical>
template <>
inline void
MeshGeomFactory<2, _triangle_3, Line_arc<cgal::Spherical>, cgal::Spherical>::
    addPrimitive(const Matrix<Real> & node_coordinates, UInt id,
                 TreeTypeHelper<Line_arc<cgal::Spherical>,
                                cgal::Spherical>::container_type & list) {

  TreeTypeHelper<Line_arc<cgal::Spherical>, cgal::Spherical>::point_type a(
      node_coordinates(0, 0), node_coordinates(1, 0), 0.),
      b(node_coordinates(0, 1), node_coordinates(1, 1), 0.),
      c(node_coordinates(0, 2), node_coordinates(1, 2), 0.);

  /*std::cout << "elem " << id << " node 1 : x_node=" << node_coordinates(0, 0)
        << ", x_arc_node=" << a.x() << ", y_node=" << node_coordinates(1, 0)
        << ", y_arc_node=" << a.y() << std::endl ;
  std::cout << "elem " << id << " node 2 : x_node=" << node_coordinates(0, 1)
        << ", x_arc_node=" << b.x() << ", y_node=" << node_coordinates(1, 1)
        << ", y_arc_node=" << b.y() << std::endl ;
  std::cout << "elem " << id << " node 2 : x_node=" << node_coordinates(0, 2)
        << ", x_arc_node=" << c.x() << ", y_node=" << node_coordinates(1, 2)
        << ", y_arc_node=" << c.y() << std::endl ;*/

  CGAL::Line_3<cgal::Spherical> l1(a, b), l2(b, c), l3(c, a);
  Line_arc<cgal::Spherical> s1(l1, a, b), s2(l2, b, c), s3(l3, c, a);

  s1.setId(id);
  s1.setSegId(0);
  s2.setId(id);
  s2.setSegId(1);
  s3.setId(id);
  s3.setSegId(2);

  list.push_back(s1);
  list.push_back(s2);
  list.push_back(s3);
}

#if defined(AKANTU_IGFEM)

// (2D, _igfem_triangle_4) decomposed into Line_arc<cgal::Spherical>
template <>
inline void MeshGeomFactory<2, _igfem_triangle_4, Line_arc<cgal::Spherical>,
                            cgal::Spherical>::
    addPrimitive(const Matrix<Real> & node_coordinates, UInt id,
                 TreeTypeHelper<Line_arc<cgal::Spherical>,
                                cgal::Spherical>::container_type & list) {

  TreeTypeHelper<Line_arc<cgal::Spherical>, cgal::Spherical>::point_type a(
      node_coordinates(0, 0), node_coordinates(1, 0), 0.),
      b(node_coordinates(0, 1), node_coordinates(1, 1), 0.),
      c(node_coordinates(0, 2), node_coordinates(1, 2), 0.);

  CGAL::Line_3<cgal::Spherical> l1(a, b), l2(b, c), l3(c, a);
  Line_arc<cgal::Spherical> s1(l1, a, b), s2(l2, b, c), s3(l3, c, a);

  s1.setId(id);
  s1.setSegId(0);
  s2.setId(id);
  s2.setSegId(1);
  s3.setId(id);
  s3.setSegId(2);

  list.push_back(s1);
  list.push_back(s2);
  list.push_back(s3);
}

// (2D, _igfem_triangle_4) decomposed into Line_arc<cgal::Spherical>
template <>
inline void MeshGeomFactory<2, _igfem_triangle_5, Line_arc<cgal::Spherical>,
                            cgal::Spherical>::
    addPrimitive(const Matrix<Real> & node_coordinates, UInt id,
                 TreeTypeHelper<Line_arc<cgal::Spherical>,
                                cgal::Spherical>::container_type & list) {

  TreeTypeHelper<Line_arc<cgal::Spherical>, cgal::Spherical>::point_type a(
      node_coordinates(0, 0), node_coordinates(1, 0), 0.),
      b(node_coordinates(0, 1), node_coordinates(1, 1), 0.),
      c(node_coordinates(0, 2), node_coordinates(1, 2), 0.);

  CGAL::Line_3<cgal::Spherical> l1(a, b), l2(b, c), l3(c, a);
  Line_arc<cgal::Spherical> s1(l1, a, b), s2(l2, b, c), s3(l3, c, a);

  s1.setId(id);
  s1.setSegId(0);
  s2.setId(id);
  s2.setSegId(1);
  s3.setId(id);
  s3.setSegId(2);

  list.push_back(s1);
  list.push_back(s2);
  list.push_back(s3);
}

#endif

// (3D, _tetrahedron_4) decomposed into Triangle<cgal::Cartesian>
template <>
inline void
MeshGeomFactory<3, _tetrahedron_4, Triangle<cgal::Cartesian>, cgal::Cartesian>::
    addPrimitive(const Matrix<Real> & node_coordinates, UInt id,
                 TreeTypeHelper<Triangle<cgal::Cartesian>,
                                cgal::Cartesian>::container_type & list) {

  TreeTypeHelper<Triangle<cgal::Cartesian>, cgal::Cartesian>::point_type a(
      node_coordinates(0, 0), node_coordinates(1, 0), node_coordinates(2, 0)),
      b(node_coordinates(0, 1), node_coordinates(1, 1), node_coordinates(2, 1)),
      c(node_coordinates(0, 2), node_coordinates(1, 2), node_coordinates(2, 2)),
      d(node_coordinates(0, 3), node_coordinates(1, 3), node_coordinates(2, 3));

  Triangle<cgal::Cartesian> t1(a, b, c), t2(b, c, d), t3(c, d, a), t4(d, a, b);

  t1.setId(id);
  t2.setId(id);
  t3.setId(id);
  t4.setId(id);

  list.push_back(t1);
  list.push_back(t2);
  list.push_back(t3);
  list.push_back(t4);
}

} // namespace akantu

#endif // __AKANTU_MESH_GEOM_FACTORY_TMPL_HH__

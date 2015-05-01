/**
 * @file   mesh_geom_factory_inline_impl.cc
 *
 * @author Lucas Frérot <lucas.frerot@epfl.ch>
 *
 * @date creation: Thu Feb 26 2015
 * @date last modification: Fri Mar 6 2015
 *
 * @brief  Class for constructing the CGAL primitives of a mesh
 *
 * @section LICENSE
 *
 * Copyright (©) 2015 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

// 2D and _triangle_3 implementation
template<>
inline void MeshGeomFactory<2, _triangle_3, Triangle<Cartesian>, Cartesian>::addPrimitive(
    const Matrix<Real> & node_coordinates,
    UInt id) {
  TreeTypeHelper<Triangle<Cartesian>, Cartesian>::point_type a(node_coordinates(0, 0), node_coordinates(1, 0), 0.);
  TreeTypeHelper<Triangle<Cartesian>, Cartesian>::point_type b(node_coordinates(0, 1), node_coordinates(1, 1), 0.);
  TreeTypeHelper<Triangle<Cartesian>, Cartesian>::point_type c(node_coordinates(0, 2), node_coordinates(1, 2), 0.);

  Triangle<Cartesian> t(a, b, c);
  t.setId(id);
  primitive_list.push_back(t);
}

// 3D and _tetrahedron_4 with triangles implementation
template<>
inline void MeshGeomFactory<3, _tetrahedron_4, Triangle<Cartesian>, Cartesian>::addPrimitive(
    const Matrix<Real> & node_coordinates,
    UInt id) {
  TreeTypeHelper<Triangle<Cartesian>, Cartesian>::point_type
    a(node_coordinates(0, 0), node_coordinates(1, 0), node_coordinates(2, 0)),
    b(node_coordinates(0, 1), node_coordinates(1, 1), node_coordinates(2, 1)),
    c(node_coordinates(0, 2), node_coordinates(1, 2), node_coordinates(2, 2)),
    d(node_coordinates(0, 3), node_coordinates(1, 3), node_coordinates(2, 3));

  Triangle<Cartesian>
    t1(a, b, c),
    t2(b, c, d),
    t3(c, d, a),
    t4(d, a, b);

  t1.setId(id);
  t2.setId(id);
  t3.setId(id);
  t4.setId(id);

  primitive_list.push_back(t1);
  primitive_list.push_back(t2);
  primitive_list.push_back(t3);
  primitive_list.push_back(t4);
}

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
/* -------------------------------------------------------------------------- */
#include "mesh_geom_common.hh"
#include "mesh_geom_factory.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_MESH_GEOM_FACTORY_TMPL_HH_
#define AKANTU_MESH_GEOM_FACTORY_TMPL_HH_

namespace akantu {

/* -------------------------------------------------------------------------- */
template <UInt dim, ElementType type, class Primitive, class Kernel>
MeshGeomFactory<dim, type, Primitive, Kernel>::MeshGeomFactory(Mesh & mesh)
    : MeshGeomAbstract(mesh) {}

/* -------------------------------------------------------------------------- */
template <UInt dim, ElementType type, class Primitive, class Kernel>
MeshGeomFactory<dim, type, Primitive, Kernel>::~MeshGeomFactory() {
  delete data_tree;
}

/* -------------------------------------------------------------------------- */
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
  auto it = connectivity.begin(nb_nodes_per_element);
  auto end = connectivity.end(nb_nodes_per_element);

  Matrix<Real> node_coordinates(dim, nb_nodes_per_element);

  // This loop builds the list of primitives
  for (; it != end; ++it, ++el_index) {
    const Vector<UInt> & el_connectivity = *it;

    for (UInt i = 0; i < nb_nodes_per_element; i++) {
      for (UInt j = 0; j < dim; j++) {
        node_coordinates(j, i) = nodes(el_connectivity(i), j);
      }
    }

    // the unique elemental id assigned to the primitive is the
    // linearized element index over ghost type
    addPrimitive(node_coordinates, el_index);
  }

  delete data_tree;

  // This condition allows the use of the mesh geom module
  // even if types are not compatible with AABB tree algorithm
  if (TreeTypeHelper_::is_valid) {
    data_tree = new TreeType(primitive_list.begin(), primitive_list.end());
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
namespace {
  namespace details {
    enum class GeometricalType {
      _triangle,
      _tetrahedron,
    };
    template <ElementType element_type> struct GeometricalTypeHelper {};

    template <> struct GeometricalTypeHelper<_triangle_3> {
      static const GeometricalType type{GeometricalType::_triangle};
    };

    template <> struct GeometricalTypeHelper<_triangle_6> {
      static const GeometricalType type{GeometricalType::_triangle};
    };

    template <> struct GeometricalTypeHelper<_tetrahedron_4> {
      static const GeometricalType type{GeometricalType::_triangle};
    };

#if defined(AKANTU_IGFEM)
    template <> struct GeometricalTypeHelper<_igfem_triangle_4> {
      static const GeometricalType type{GeometricalType::_triangle};
    };
    template <> struct GeometricalTypeHelper<_igfem_triangle_5> {
      static const GeometricalType type{GeometricalType::_triangle};
    };
#endif

    template <details::GeometricalType geom_type, class Primitive, class Kernel>
    struct AddPrimitiveHelper {};

    template <class Primitive>
    struct AddPrimitiveHelper<GeometricalType::_triangle, Primitive,
                              cgal::Cartesian> {
      using TreeTypeHelper_ = TreeTypeHelper<Primitive, cgal::Cartesian>;
      using ContainerType = typename TreeTypeHelper_::container_type;
      static void addPrimitive(const Matrix<Real> & node_coordinates, UInt id,
                               ContainerType & list) {
        using Point = typename TreeTypeHelper_::point_type;
        Point a(node_coordinates(0, 0), node_coordinates(1, 0), 0.);
        Point b(node_coordinates(0, 1), node_coordinates(1, 1), 0.);
        Point c(node_coordinates(0, 2), node_coordinates(1, 2), 0.);

        Triangle<cgal::Cartesian> t(a, b, c);
        t.setId(id);
        list.push_back(t);
      }
    };

    template <class Primitive>
    struct AddPrimitiveHelper<GeometricalType::_triangle, Primitive,
                              cgal::Spherical> {
      using TreeTypeHelper_ = TreeTypeHelper<Primitive, cgal::Spherical>;
      using ContainerType = typename TreeTypeHelper_::container_type;
      static void addPrimitive(const Matrix<Real> & node_coordinates, UInt id,
                               ContainerType & list) {
        using Point = typename TreeTypeHelper_::point_type;
        Point a(node_coordinates(0, 0), node_coordinates(1, 0), 0.);
        Point b(node_coordinates(0, 1), node_coordinates(1, 1), 0.);
        Point c(node_coordinates(0, 2), node_coordinates(1, 2), 0.);

        using Line = CGAL::Line_3<cgal::Spherical>;
        Line l1(a, b);
        Line l2(b, c);
        Line l3(c, a);

        using Arc = Line_arc<cgal::Spherical>;
        Arc s1(l1, a, b);
        Arc s2(l2, b, c);
        Arc s3(l3, c, a);

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
    };

    template <class Primitive>
    struct AddPrimitiveHelper<GeometricalType::_tetrahedron, Primitive,
                              cgal::Cartesian> {
      using TreeTypeHelper_ = TreeTypeHelper<Primitive, cgal::Cartesian>;
      using ContainerType = typename TreeTypeHelper_::container_type;
      static void addPrimitive(const Matrix<Real> & node_coordinates, UInt id,
                               ContainerType & list) {
        using Point = typename TreeTypeHelper_::point_type;
        Point a(node_coordinates(0, 0), node_coordinates(1, 0),
                node_coordinates(2, 0));
        Point b(node_coordinates(0, 1), node_coordinates(1, 1),
                node_coordinates(2, 1));
        Point c(node_coordinates(0, 2), node_coordinates(1, 2),
                node_coordinates(2, 2));
        Point d(node_coordinates(0, 3), node_coordinates(1, 3),
                node_coordinates(2, 3));

        Triangle<cgal::Cartesian> t1(a, b, c);
        Triangle<cgal::Cartesian> t2(b, c, d);
        Triangle<cgal::Cartesian> t3(c, d, a);
        Triangle<cgal::Cartesian> t4(d, a, b);

        t1.setId(id);
        t2.setId(id);
        t3.setId(id);
        t4.setId(id);

        list.push_back(t1);
        list.push_back(t2);
        list.push_back(t3);
        list.push_back(t4);
      }
    };
  } // namespace details
} // namespace

/* -------------------------------------------------------------------------- */
template <UInt dim, ElementType type, class Primitive, class Kernel>
void MeshGeomFactory<dim, type, Primitive, Kernel>::addPrimitive(
    const Matrix<Real> & node_coordinates, UInt id, ContainerType & list) {
  details::AddPrimitiveHelper<details::GeometricalTypeHelper<type>::type,
                              Primitive, Kernel>::addPrimitive(node_coordinates,
                                                               id, list);
}

/* -------------------------------------------------------------------------- */
template <UInt dim, ElementType type, class Primitive, class Kernel>
void MeshGeomFactory<dim, type, Primitive, Kernel>::addPrimitive(
    const Matrix<Real> & node_coordinates, UInt id) {
  this->addPrimitive(node_coordinates, id, this->primitive_list);
}

} // namespace akantu

#endif // AKANTU_MESH_GEOM_FACTORY_TMPL_HH_

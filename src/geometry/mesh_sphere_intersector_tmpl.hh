/**
 * @file   mesh_sphere_intersector_tmpl.hh
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 * @author Clement Roux <clement.roux@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Tue Jun 23 2015
 * @date last modification: Tue Feb 20 2018
 *
 * @brief  Computation of mesh intersection with spheres
 *
 * @section LICENSE
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

#ifndef __AKANTU_MESH_SPHERE_INTERSECTOR_TMPL_HH__
#define __AKANTU_MESH_SPHERE_INTERSECTOR_TMPL_HH__

#include "aka_common.hh"
#include "mesh_geom_common.hh"
#include "mesh_sphere_intersector.hh"
#include "tree_type_helper.hh"

namespace akantu {

template <UInt dim, ElementType type>
MeshSphereIntersector<dim, type>::MeshSphereIntersector(Mesh & mesh)
    : parent_type(mesh), tol_intersection_on_node(1e-10) {
#if defined(AKANTU_IGFEM)
  if ((type == _triangle_3) || (type == _igfem_triangle_4) ||
      (type == _igfem_triangle_5)) {
    const_cast<UInt &>(this->nb_seg_by_el) = 3;
  } else {
    AKANTU_ERROR("Not ready for mesh type " << type);
  }
#else
  if ((type != _triangle_3))
    AKANTU_ERROR("Not ready for mesh type " << type);
#endif

  // initialize the intersection pointsss array with the spatial dimension
  this->intersection_points = new Array<Real>(0, dim);
  //  A maximum is set to the number of intersection nodes per element to limit
  //  the size of new_node_per_elem: 2 in 2D and 4 in 3D
  this->new_node_per_elem = new Array<UInt>(0, 1 + 4 * (dim - 1));
}

template <UInt dim, ElementType type>
MeshSphereIntersector<dim, type>::~MeshSphereIntersector() {
  delete this->new_node_per_elem;
  delete this->intersection_points;
}

template <UInt dim, ElementType type>
void MeshSphereIntersector<dim, type>::constructData(GhostType ghost_type) {

  this->new_node_per_elem->resize(this->mesh.getNbElement(type, ghost_type));
  this->new_node_per_elem->clear();

  MeshGeomIntersector<dim, type, Line_arc<SK>, SK::Sphere_3, SK>::constructData(
      ghost_type);
}

template <UInt dim, ElementType type>
void MeshSphereIntersector<dim, type>::computeMeshQueryIntersectionPoint(
    const SK::Sphere_3 & query, UInt nb_old_nodes) {
  /// function to replace computeIntersectionQuery in a more generic geometry
  /// module version
  // The newNodeEvent is not send from this method who only compute the
  // intersection points
  AKANTU_DEBUG_IN();

  Array<Real> & nodes = this->mesh.getNodes();
  UInt nb_node = nodes.size() + this->intersection_points->size();

  // Tolerance for proximity checks should be defined by user
  Real global_tolerance = Math::getTolerance();
  Math::setTolerance(tol_intersection_on_node);
  typedef boost::variant<pair_type> sk_inter_res;

  TreeTypeHelper<Line_arc<cgal::Spherical>, cgal::Spherical>::const_iterator
      it = this->factory.getPrimitiveList().begin(),
      end = this->factory.getPrimitiveList().end();

  for (; it != end; ++it) { // loop on the primitives (segments)
    std::list<sk_inter_res> s_results;
    CGAL::intersection(*it, query, std::back_inserter(s_results));

    if (s_results.size() == 1) { // just one point
      if (pair_type * pair = boost::get<pair_type>(&s_results.front())) {
        if (pair->second == 1) { // not a point tangent to the sphere
          // the intersection point written as a vector
          Vector<Real> new_node(dim, 0.0);
          cgal::Cartesian::Point_3 point(CGAL::to_double(pair->first.x()),
                                         CGAL::to_double(pair->first.y()),
                                         CGAL::to_double(pair->first.z()));
          for (UInt i = 0; i < dim; i++) {
            new_node(i) = point[i];
          }

          /// boolean to decide wheter intersection point is on a standard node
          /// of the mesh or not
          bool is_on_mesh = false;
          /// boolean to decide if this intersection point has been already
          /// computed for a neighbor element
          bool is_new = true;

          /// check if intersection point has already been computed
          UInt n = nb_old_nodes;

          // check if we already compute this intersection and add it as a node
          // for a neighboor element of another type
          auto existing_node = nodes.begin(dim);

          for (; n < nodes.size(); ++n) { // loop on the nodes from nb_old_nodes
            if (Math::are_vector_equal(dim, new_node.storage(),
                                       existing_node[n].storage())) {
              is_new = false;
              break;
            }
          }
          if (is_new) {
            auto intersection_points_it = this->intersection_points->begin(dim);
            auto intersection_points_end = this->intersection_points->end(dim);
            for (; intersection_points_it != intersection_points_end;
                 ++intersection_points_it, ++n) {
              if (Math::are_vector_equal(dim, new_node.storage(),
                                         intersection_points_it->storage())) {
                is_new = false;
                break;
              }
            }
          }

          // get the initial and final points of the primitive (segment) and
          // write them as vectors
          cgal::Cartesian::Point_3 source_cgal(
              CGAL::to_double(it->source().x()),
              CGAL::to_double(it->source().y()),
              CGAL::to_double(it->source().z()));
          cgal::Cartesian::Point_3 target_cgal(
              CGAL::to_double(it->target().x()),
              CGAL::to_double(it->target().y()),
              CGAL::to_double(it->target().z()));
          Vector<Real> source(dim), target(dim);
          for (UInt i = 0; i < dim; i++) {
            source(i) = source_cgal[i];
            target(i) = target_cgal[i];
          }

          // Check if we are close from a node of the primitive (segment)
          if (Math::are_vector_equal(dim, source.storage(),
                                     new_node.storage()) ||
              Math::are_vector_equal(dim, target.storage(),
                                     new_node.storage())) {
            is_on_mesh = true;
            is_new = false;
          }

          if (is_new) { // if the intersection point is a new one add it to the
                        // list
            this->intersection_points->push_back(new_node);
            nb_node++;
          }

          // deduce the element id
          UInt element_id = it->id();

          // fill the new_node_per_elem array
          if (!is_on_mesh) { // if the node is not on a mesh node
            UInt & nb_new_nodes_per_el =
                (*this->new_node_per_elem)(element_id, 0);
            nb_new_nodes_per_el += 1;
            AKANTU_DEBUG_ASSERT(
                2 * nb_new_nodes_per_el <
                    this->new_node_per_elem->getNbComponent(),
                "You might have to interface crossing the same material");
            (*this->new_node_per_elem)(element_id,
                                       (2 * nb_new_nodes_per_el) - 1) = n;
            (*this->new_node_per_elem)(element_id, 2 * nb_new_nodes_per_el) =
                it->segId();
          }
        }
      }
    }
  }

  Math::setTolerance(global_tolerance);

  AKANTU_DEBUG_OUT();
}

} // namespace akantu

#endif // __AKANTU_MESH_SPHERE_INTERSECTOR_TMPL_HH__

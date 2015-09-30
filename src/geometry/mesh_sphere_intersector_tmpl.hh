/**
 * @file mesh_sphere_intersector_tmpl.hh
 *
 * @author Clément Roux-Langlois <clement.roux@epfl.ch>
 *
 * @date creation: Wed june 10 2015
 * @date last modification: Wed June 17 2015
 *
 * @brief Computation of mesh intersection with spheres
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

#ifndef __AKANTU_MESH_SPHERE_INTERSECTOR_TMPL_HH__
#define __AKANTU_MESH_SPHERE_INTERSECTOR_TMPL_HH__

#include "aka_common.hh"
#include "mesh_geom_common.hh"
#include "tree_type_helper.hh"
#include "mesh_sphere_intersector.hh"
#include "static_communicator.hh"

__BEGIN_AKANTU__

template<UInt dim, ElementType type>
MeshSphereIntersector<dim, type>::MeshSphereIntersector(Mesh & mesh,
							const ID & id,
							const MemoryID & memory_id):
  parent_type(mesh, id, memory_id),
  tol_intersection_on_node(1e-10),
  nb_nodes_fem(mesh.getNodes().getSize()),
  nb_prim_by_el(0)
{
  this->intersection_points = new Array<Real>(0,dim);

  for(Mesh::type_iterator it = mesh.firstType(dim); it != mesh.lastType(dim); ++it){
#if defined(AKANTU_IGFEM)
    if(*it == _triangle_3){
      // Addition of the igfem types in the mesh
      this->mesh.addConnectivityType(_igfem_triangle_4, _not_ghost);
      this->mesh.addConnectivityType(_igfem_triangle_4, _ghost);
      this->mesh.addConnectivityType(_igfem_triangle_5, _not_ghost);
      this->mesh.addConnectivityType(_igfem_triangle_5, _ghost);
    }

    else if( (*it != _triangle_3) && (*it != _igfem_triangle_4) && (*it != _igfem_triangle_5) ) {
      AKANTU_DEBUG_ERROR("Not ready for mesh type " << *it);
    }

    if( (type == _triangle_3) || (type == _igfem_triangle_4) || (type == _igfem_triangle_5) ){
      this->nb_prim_by_el = 3;
      const_cast<UInt &>(this->nb_seg_by_el) = 3;
    } else {
      AKANTU_DEBUG_ERROR("Not ready for mesh type " << type);
    }
#else
    if( (*it != _triangle_3) )
      AKANTU_DEBUG_ERROR("Not ready for mesh type " << *it);
    this->nb_prim_by_el = 3;
#endif
  }

  for (ghost_type_t::iterator gt = ghost_type_t::begin();  gt != ghost_type_t::end(); ++gt) {
    GhostType ghost_type = *gt;
    UInt nb_comp = 1 + nb_prim_by_el * 2;
    this->new_node_per_elem.alloc(0, nb_comp, type, ghost_type, 0);
  }

  this->constructData();
}

  template<UInt dim, ElementType type>
MeshSphereIntersector<dim, type>::~MeshSphereIntersector()
{}

template<UInt dim, ElementType type>
void MeshSphereIntersector<dim, type>::constructData() {

  for (ghost_type_t::iterator gt = ghost_type_t::begin();  gt != ghost_type_t::end(); ++gt) {
    GhostType ghost_type = *gt;
    this->new_node_per_elem(type, ghost_type).resize(this->mesh.getNbElement(type, ghost_type));
    this->new_node_per_elem(type, ghost_type).clear();
  }

  MeshGeomIntersector<dim, type, Line_arc<SK>, SK::Sphere_3, SK>::constructData();
}

template<UInt dim, ElementType type>
void MeshSphereIntersector<dim, type>:: computeMeshQueryIntersectionPoint(const SK::Sphere_3 & query) {
  /// function to replace computeIntersectionQuery in a more generic geometry module version
  // The newNodeEvent is not send from this method who only compute the intersection points
  AKANTU_DEBUG_IN();

  Array<Real> & nodes = this->mesh.getNodes();
  UInt nb_node = nodes.getSize() + this->intersection_points->getSize();
  UInt nb_not_ghost_elements = this->mesh.getNbElement(type);

  // Tolerance for proximity checks should be defined by user
  Math::setTolerance(tol_intersection_on_node);
  typedef boost::variant<pair_type> sk_inter_res;

  TreeTypeHelper<Line_arc<Spherical>, Spherical>::const_iterator
    it = this->factory.getPrimitiveList().begin(),
    end= this->factory.getPrimitiveList().end();

  for (; it != end ; ++it) {
    std::list<sk_inter_res> s_results;
    CGAL::intersection(*it, query, std::back_inserter(s_results));

    if (s_results.size() == 1) { // just one point
      if (pair_type * pair = boost::get<pair_type>(&s_results.front())) {
        if (pair->second == 1) { // not a point tangent to the sphere
          // Addition of the new node
          Vector<Real> new_node(dim, 0.0);
          Cartesian::Point_3 point(CGAL::to_double(pair->first.x()),
                                   CGAL::to_double(pair->first.y()),
                                   CGAL::to_double(pair->first.z()));

          for (UInt i = 0 ; i < dim ; i++) {
            new_node(i) = point[i];
          }

          bool is_on_mesh = false, is_new = true;
          // check if we already compute this intersection for a neighboor element
	  UInt n = nb_nodes_fem;
	  Array<Real>::vector_iterator existing_node = nodes.begin(dim);
	  for (; n < nodes.getSize() ; ++n) {
	    if (Math::are_vector_equal(dim, new_node.storage(), existing_node[n].storage())) {
	      is_new = false;
	      break;
	    }
	  }
	  if(is_new){
	    Array<Real>::vector_iterator intersection_points_it = this->intersection_points->begin(dim);
	    Array<Real>::vector_iterator intersection_points_end = this->intersection_points->end(dim);
	    for (; intersection_points_it != intersection_points_end ; ++intersection_points_it, ++n) {
	      if (Math::are_vector_equal(dim, new_node.storage(), intersection_points_it->storage())) {
		is_new = false;
		break;
	      }
	    }
	  }

	  Cartesian::Point_3 source_cgal(CGAL::to_double(it->source().x()),
					 CGAL::to_double(it->source().y()),
					 CGAL::to_double(it->source().z()));
	  Cartesian::Point_3 target_cgal(CGAL::to_double(it->target().x()),
					 CGAL::to_double(it->target().y()),
					 CGAL::to_double(it->target().z()));

	  Vector<Real> source(dim), target(dim);
	  for (UInt i = 0 ; i < dim ; i++) {
	    source(i) = source_cgal[i];
	    target(i) = target_cgal[i];
	  }

	  // Check if we are close from a node of the segment
	  if (Math::are_vector_equal(dim, source.storage(), new_node.storage()) ||
	      Math::are_vector_equal(dim, target.storage(), new_node.storage())) {
	    is_on_mesh = true;
	    is_new = false;
	  }

	  if (is_new) {
	    this->intersection_points->push_back(new_node);
	    nb_node++;
	  }

	  // deduce ghost type and element id
	  UInt element_id = it->id();
	  GhostType ghost_type = _not_ghost;
	  if (element_id >= nb_not_ghost_elements) {
	    element_id -= nb_not_ghost_elements;
	    ghost_type = _ghost;
	  }

	  Array<UInt> & new_node_per_elem_array = this->new_node_per_elem(type, ghost_type);

	  if (!is_on_mesh) {
	    new_node_per_elem_array(element_id, 0) += 1;
	    new_node_per_elem_array(element_id, (2 * new_node_per_elem_array(element_id, 0)) - 1) = n;
	    new_node_per_elem_array(element_id, 2 * new_node_per_elem_array(element_id, 0)) = it->segId();
	  } else {
	    // if intersection is at a node, write node number (in el) in pennultimate position
	    if (Math::are_vector_equal(dim, source.storage(), new_node.storage())) {
	      new_node_per_elem_array(element_id, (new_node_per_elem_array.getNbComponent() - 2)) = it->segId() - 1;
	    } else {
	      new_node_per_elem_array(element_id, (new_node_per_elem_array.getNbComponent() - 2)) =
		it->segId() % this->nb_prim_by_el;
	    }
	  }
	}
      }
    }
  }

  AKANTU_DEBUG_OUT();
}
__END_AKANTU__

#endif // __AKANTU_MESH_SPHERE_INTERSECTOR_TMPL_HH__

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

__BEGIN_AKANTU__

template<UInt dim, ElementType type>
MeshSphereIntersector<dim, type>::MeshSphereIntersector(Mesh & mesh):
  MeshGeomIntersector<dim, type, Line_arc<SK>, SK::Sphere_3, SK>(mesh),
  tol_intersection_on_node(1e-10),
  nb_nodes_fem(mesh.getNodes().getSize()),
  nb_prim_by_el(0)
{
  this->new_node_per_elem = new Array<UInt>(0, 1 + 4*(dim-1));
  this->new_nodes = new Array<Real>(0,dim);
  this->constructData();

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
#endif
  }
}

  template<UInt dim, ElementType type>
MeshSphereIntersector<dim, type>::~MeshSphereIntersector()
{}

template<UInt dim, ElementType type>
void MeshSphereIntersector<dim, type>::constructData() {

  this->new_node_per_elem->resize(this->mesh.getConnectivity(type).getSize());
  this->new_node_per_elem->clear();

  MeshGeomIntersector<dim, type, Line_arc<SK>, SK::Sphere_3, SK>::constructData();
}

template<UInt dim, ElementType type>
void MeshSphereIntersector<dim, type>::computeIntersectionQuery(const SK::Sphere_3 & query) {
  AKANTU_DEBUG_IN();

  Array<Real> & nodes = this->mesh.getNodes();
  UInt nb_node = nodes.getSize() ;
  Array<UInt> & new_node_per_elem = *(this->new_node_per_elem);

  // Tolerance for proximity checks should be defined by user
  Math::setTolerance(tol_intersection_on_node);
  typedef boost::variant<pair_type> sk_inter_res;

  TreeTypeHelper<Line_arc<Spherical>, Spherical>::const_iterator
    it = this->factory.getPrimitiveList().begin(),
    end= this->factory.getPrimitiveList().end();

  NewNodesEvent new_nodes;
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
	  UInt n = nb_nodes_fem-1;
          for (; n < nb_node ; ++n) {
            Array<Real>::vector_iterator existing_node = nodes.begin(dim) + n;
            if (Math::are_vector_equal(dim, new_node.storage(), existing_node->storage())) {
              is_new = false;
              break;
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
              nodes.push_back(new_node);
              new_nodes.getList().push_back(nb_node);
              nb_node++;
            }

            if (!is_on_mesh) {
              new_node_per_elem(it->id(), 0) += 1;
              new_node_per_elem(it->id(), (2 * new_node_per_elem(it->id(), 0)) - 1) = n;
              new_node_per_elem(it->id(), 2 * new_node_per_elem(it->id(), 0)) = it->segId();
            }

            else {
              // if intersection is at a node, write node number (in el) in pennultimate position
              if (Math::are_vector_equal(dim, source.storage(), new_node.storage())) {
                new_node_per_elem(it->id(), (new_node_per_elem.getNbComponent() - 2)) = it->segId() - 1;
              } else {
                new_node_per_elem(it->id(), (new_node_per_elem.getNbComponent() - 2)) =
                  it->segId() % this->nb_prim_by_el;
              }
            }
          }
        }
      }
  }

  if(new_nodes.getList().getSize())
    this->mesh.sendEvent(new_nodes);

  AKANTU_DEBUG_OUT();
}

template<UInt dim, ElementType type>
void MeshSphereIntersector<dim, type>:: computeMeshQueryIntersectionPoint(const SK::Sphere_3 & query) {
  /// function to replace computeIntersectionQuery in a more generic geometry module version
  // The newNodeEvent is not send from this method who only compute the intersection points
  AKANTU_DEBUG_IN();

  Array<Real> & nodes = this->mesh.getNodes();
  UInt nb_node = nodes.getSize() + this->new_nodes->getSize();
  Array<UInt> & new_node_per_elem = *(this->new_node_per_elem);

  // Tolerance for proximity checks should be defined by user
  Math::setTolerance(tol_intersection_on_node);
  typedef boost::variant<pair_type> sk_inter_res;

  TreeTypeHelper<Line_arc<Spherical>, Spherical>::const_iterator
    it = this->factory.getPrimitiveList().begin(),
    end= this->factory.getPrimitiveList().end();

  //NewNodesEvent new_nodes;
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
	  UInt n = nb_nodes_fem-1;
          for (; n < nb_node ; ++n) {
            Array<Real>::vector_iterator existing_node = nodes.begin(dim) + n;
            if (Math::are_vector_equal(dim, new_node.storage(), existing_node->storage())) {
              is_new = false;
              break;
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
              this->new_nodes->push_back(new_node);
              // new_nodes.getList().push_back(nb_node);
              nb_node++;
            }

            if (!is_on_mesh) {
              new_node_per_elem(it->id(), 0) += 1;
              new_node_per_elem(it->id(), (2 * new_node_per_elem(it->id(), 0)) - 1) = n;
              new_node_per_elem(it->id(), 2 * new_node_per_elem(it->id(), 0)) = it->segId();
            }

            else {
              // if intersection is at a node, write node number (in el) in pennultimate position
              if (Math::are_vector_equal(dim, source.storage(), new_node.storage())) {
                new_node_per_elem(it->id(), (new_node_per_elem.getNbComponent() - 2)) = it->segId() - 1;
              } else {
                new_node_per_elem(it->id(), (new_node_per_elem.getNbComponent() - 2)) =
                  it->segId() % this->nb_prim_by_el;
              }
            }
          }
        }
      }
  }

  /*if(new_nodes.getList().getSize())
    this->mesh.sendEvent(new_nodes);*/

  AKANTU_DEBUG_OUT();
}

template<UInt dim, ElementType type>
void MeshSphereIntersector<dim, type>::removeAdditionnalNodes() {
  AKANTU_DEBUG_IN();

  RemovedNodesEvent remove_nodes(this->mesh);
  Array<UInt> & nodes_removed = remove_nodes.getList();
  Array<UInt> & new_numbering = remove_nodes.getNewNumbering();

  for(UInt i = 0 ; i < nb_nodes_fem ; ++i){
    new_numbering(i) = i;
  }

  for(UInt i = nb_nodes_fem ; i < this->mesh.getNodes().getSize(); ++i){
    new_numbering(i) = UInt(-1) ;
    nodes_removed.push_back(i);
  }

  this->mesh.sendEvent(remove_nodes);
  AKANTU_DEBUG_OUT();
}


#if defined(AKANTU_IGFEM)

template<UInt dim, ElementType type>
void MeshSphereIntersector<dim, type>::buildResultFromQueryList(const std::list<SK::Sphere_3> & query_list) {
  AKANTU_DEBUG_IN();

  this->computeIntersectionQueryList(query_list);

  Array<UInt> & connec_type_tmpl = this->mesh.getConnectivity(type);
  const Array<UInt> new_node_per_elem = *(this->new_node_per_elem);

  Array<UInt>
    & connec_igfem_tri4 = this->mesh.getConnectivity(_igfem_triangle_4),
    & connec_igfem_tri5 = this->mesh.getConnectivity(_igfem_triangle_5),
    & connec_tri3 = this->mesh.getConnectivity(_triangle_3);

  Element element_tri3(_triangle_3, 0, _not_ghost, _ek_regular);
  Element element_tri4(_igfem_triangle_4, 0, _not_ghost, _ek_igfem);
  Element element_tri5(_igfem_triangle_5, 0, _not_ghost, _ek_igfem);

  NewIGFEMElementsEvent new_elements;
  UInt n_new_el = 0;
  Array<Element> & new_elements_list = new_elements.getList();
  Array<Element> & old_elements_list = new_elements.getOldElementsList();
 
  RemovedElementsEvent remove_elem(this->mesh);
  remove_elem.getNewNumbering().alloc(connec_type_tmpl.getSize(), 1, type, _not_ghost);
  Array<UInt> & new_numbering = remove_elem.getNewNumbering(type, _not_ghost);
  UInt new_nb_type_tmpl = 0; // Meter of the number of element (type) will we loop on elements

  for (UInt nel = 0 ; nel != new_node_per_elem.getSize(); ++nel) {
    if( (type != _triangle_3)
        && ( (new_node_per_elem(nel,0)==0)
             || ( (new_node_per_elem(nel,0) == 1) 
                  && ( ( (new_node_per_elem(nel,2)+1) % this->nb_prim_by_el ) 
                       != new_node_per_elem(nel, new_node_per_elem.getNbComponent() - 2) ) ) ) ){
      Element element_type_tmpl(type, 0, _not_ghost, Mesh::getKind(type));  
      new_elements_list.resize(n_new_el+1);
      old_elements_list.resize(n_new_el+1);
      Vector<UInt> ctri3(3);
      ctri3(0) =  connec_type_tmpl(nel,0);
      ctri3(1) =  connec_type_tmpl(nel,1);
      ctri3(2) =  connec_type_tmpl(nel,2);
      /// add the new element in the mesh
      UInt nb_tri3 = connec_tri3.getSize();
      connec_tri3.push_back(ctri3);
      element_tri3.element = nb_tri3;
      new_elements_list(n_new_el) = element_tri3;
      /// the number of the old element in the mesh
      element_type_tmpl.element = nel;
      old_elements_list(n_new_el) = element_type_tmpl;
      new_numbering(nel) =  UInt(-1);
      ++n_new_el;

      remove_elem.getList().push_back(element_type_tmpl);
  
    }
    else if( (new_node_per_elem(nel,0)!=0)
             && !( (new_node_per_elem(nel,0) == 1) 
                   && ( ( (new_node_per_elem(nel,2)+1) % this->nb_prim_by_el ) 
                        != new_node_per_elem(nel, new_node_per_elem.getNbComponent() - 2) ) ) ){
      Element element_type_tmpl(type, 0, _not_ghost);
      element_type_tmpl.kind = Mesh::getKind(type);
      new_elements_list.resize(n_new_el+1);
      old_elements_list.resize(n_new_el+1);
      switch(new_node_per_elem(nel,0)){
      case 1 :{
	Vector<UInt> ctri4(4);
	switch(new_node_per_elem(nel,2)){
	case 1 :
	  ctri4(0) = connec_type_tmpl(nel,2);
	  ctri4(1) = connec_type_tmpl(nel,0);
	  ctri4(2) = connec_type_tmpl(nel,1);
	  break;
	case 2 :
	  ctri4(0) = connec_type_tmpl(nel,0);
	  ctri4(1) = connec_type_tmpl(nel,1);
	  ctri4(2) = connec_type_tmpl(nel,2);
	  break;
	case 3 :
	  ctri4(0) = connec_type_tmpl(nel,1);
	  ctri4(1) = connec_type_tmpl(nel,2);
	  ctri4(2) = connec_type_tmpl(nel,0);
	  break;
	default :
	  AKANTU_DEBUG_ERROR("A triangle have only 3 segments not "<< new_node_per_elem(nel,2));
	  break;
	}
	ctri4(3) = new_node_per_elem(nel,1);
	UInt nb_tri4 = connec_igfem_tri4.getSize();
	connec_igfem_tri4.push_back(ctri4);
	element_tri4.element = nb_tri4;
	//new_elements.getList().push_back(element_tri4);
	new_elements_list(n_new_el) = element_tri4;
	if(type == _igfem_triangle_4)
	   new_numbering.push_back(connec_igfem_tri4.getSize()-2);
	break;
      }
      case 2 :{
	Vector<UInt> ctri5(5);
	if( (new_node_per_elem(nel,2)==1) && (new_node_per_elem(nel,4)==2) ){
	  ctri5(0) = connec_type_tmpl(nel,1);
	  ctri5(1) = connec_type_tmpl(nel,2);
	  ctri5(2) = connec_type_tmpl(nel,0);
	  ctri5(3) = new_node_per_elem(nel,3);
	  ctri5(4) = new_node_per_elem(nel,1);
	}
	else if((new_node_per_elem(nel,2)==1) && (new_node_per_elem(nel,4)==3)){
	  ctri5(0) = connec_type_tmpl(nel,0);
	  ctri5(1) = connec_type_tmpl(nel,1);
	  ctri5(2) = connec_type_tmpl(nel,2);
	  ctri5(3) = new_node_per_elem(nel,1);
	  ctri5(4) = new_node_per_elem(nel,3);
	}
	else if((new_node_per_elem(nel,2)==2) && (new_node_per_elem(nel,4)==3)){
	  ctri5(0) = connec_type_tmpl(nel,2);
	  ctri5(1) = connec_type_tmpl(nel,0);
	  ctri5(2) = connec_type_tmpl(nel,1);
	  ctri5(3) = new_node_per_elem(nel,3);
	  ctri5(4) = new_node_per_elem(nel,1);
	}
	else if((new_node_per_elem(nel,2)==2) && (new_node_per_elem(nel,4)==1)){
	  ctri5(0) = connec_type_tmpl(nel,1);
	  ctri5(1) = connec_type_tmpl(nel,2);
	  ctri5(2) = connec_type_tmpl(nel,0);
	  ctri5(3) = new_node_per_elem(nel,1);
	  ctri5(4) = new_node_per_elem(nel,3);
	}
	else if((new_node_per_elem(nel,2)==3) && (new_node_per_elem(nel,4)==1)){
	  ctri5(0) = connec_type_tmpl(nel,0);
	  ctri5(1) = connec_type_tmpl(nel,1);
	  ctri5(2) = connec_type_tmpl(nel,2);
	  ctri5(3) = new_node_per_elem(nel,3);
	  ctri5(4) = new_node_per_elem(nel,1);
	}
	else if((new_node_per_elem(nel,2)==3) && (new_node_per_elem(nel,4)==2)){
	  ctri5(0) = connec_type_tmpl(nel,2);
	  ctri5(1) = connec_type_tmpl(nel,0);
	  ctri5(2) = connec_type_tmpl(nel,1);
	  ctri5(3) = new_node_per_elem(nel,1);
	  ctri5(4) = new_node_per_elem(nel,3);
	}
	else{
	  AKANTU_DEBUG_ERROR("A triangle have only 3 segments not "<< new_node_per_elem(nel,2) << "and" << new_node_per_elem(nel,4));
	}	      
	UInt nb_tri5 = connec_igfem_tri5.getSize();
	connec_igfem_tri5.push_back(ctri5);
	element_tri5.element = nb_tri5;
	//new_elements.getList().push_back(element_tri5);
	new_elements_list(n_new_el) = element_tri5;
	if(type == _igfem_triangle_5){
	  new_numbering.push_back(new_nb_type_tmpl);
	  ++new_nb_type_tmpl;
	}
	break;
      }
      default:
	AKANTU_DEBUG_ERROR("Igfem cannot add "<< new_node_per_elem(nel,0) << " nodes to triangles");
	break;
      }
      element_type_tmpl.element = nel;
      old_elements_list(n_new_el) = element_type_tmpl;
      remove_elem.getList().push_back(element_type_tmpl);
      new_numbering(nel) =  UInt(-1);
      ++n_new_el;
    }
    else{
      new_numbering(nel) = new_nb_type_tmpl;
      ++new_nb_type_tmpl;
    }
  }

  // for(UInt nel=new_node_per_elem.getSize(); nel < this->mesh.getNbElement(type); ++nel) {
  //     new_numbering(nel) = new_nb_type_tmpl;
  //     ++new_nb_type_tmpl;
  // }

  UInt el_index = 0;
  for (UInt e = 0; e < this->mesh.getNbElement(type); ++e) {
    if (new_numbering(e) != UInt(-1)) {
      new_numbering(e) = el_index;
      ++el_index;
    }
  }
    
  if(n_new_el > 0){
    this->mesh.sendEvent(new_elements);
    this->mesh.sendEvent(remove_elem);
  }

  AKANTU_DEBUG_OUT();
}

#else

template<UInt dim, ElementType type>
void MeshSphereIntersector<dim, type>::buildResultFromQueryList(const std::list<SK::Sphere_3> & query_list) {
  AKANTU_DEBUG_TO_IMPLEMENT();
}


#endif

__END_AKANTU__

#endif // __AKANTU_MESH_SPHERE_INTERSECTOR_TMPL_HH__


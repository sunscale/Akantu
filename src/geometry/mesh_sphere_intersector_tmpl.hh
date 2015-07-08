/**
 * @file mesh_sphere_intersector_tmpl.hh
 *
 * @author Clement Roux-Langlois <clement.roux@epfl.ch>
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

__BEGIN_AKANTU__

template<UInt dim, ElementType type>
MeshSphereIntersector<dim, type>::MeshSphereIntersector(const Mesh & mesh):
  parent_type(mesh),
  new_node_per_elem(mesh.getConnectivity(_triangle_3).getSize(), 1 + 4*(dim-1)),
  nb_nodes_init(mesh.getNodes().getSize())
{
  this->constructData();
  this->new_node_per_elem.clear();

  if( ++mesh.firstType() != mesh.lastType() )
    AKANTU_DEBUG_ERROR("first_type : " << *(mesh.firstType()) <<
		       ", last_type : " << *(mesh.lastType()) <<
		       "Not ready for mesh type different than triangle_3");
}

template<UInt dim, ElementType type>
MeshSphereIntersector<dim, type>::~MeshSphereIntersector()
{}

template<UInt dim, ElementType type>
void MeshSphereIntersector<dim, type>::computeIntersectionQuery(const SK::Sphere_3 & query) {
  AKANTU_DEBUG_IN();

  Array<Real> & nodes = const_cast<Array<Real> &>(this->mesh.getNodes());
  UInt nb_node = nodes.getSize() ;
  Real tol = 1e-10;
  ///Array<UInt> & connectivity_tri3 = this->mesh.getConnectivity(_triangle_3);
  ///Array<UInt> & connectivity_IGFEMtri4 = this->mesh.getConnectivity(_IGFEM_triangle_3);
  typedef boost::variant<pair_type> sk_inter_res;

  TreeTypeHelper<Line_arc<Spherical>, Spherical>::const_iterator
    it = this->factory.getPrimitiveList().begin(),
    end= this->factory.getPrimitiveList().end();

  NewNodesEvent new_nodes;
  for(; it!=end; ++it){
    std::list<sk_inter_res> s_results;
    CGAL::intersection(*it, query, std::back_inserter(s_results));
    if(s_results.size()==1){ // just one point
      if (pair_type * pair = boost::get<pair_type>(&s_results.front())){
	if(pair->second==1){ // not a point tangent to the sphere
	  // Addition of the new node
	  Vector<Real> new_node(dim, 0.0);
	  SK::Circular_arc_point_3 arc_point = pair->first;
	  new_node(0) = arc_point.x();
	  new_node(1) = arc_point.y();
	  bool is_on_mesh = false, is_new = true;
	  UInt n = nb_nodes_init-1;
	  for(; n != nb_node; ++n){
	    if( ( pow((new_node(0)-nodes(n,0)),2) +
		  pow((new_node(1)-nodes(n,1)),2) ) < pow(tol,2) ){
	      is_new = false;
	      break;
	    }
	  }
	  Real x1 = it->source().x(), y1 = it->source().y(), x2 = it->target().x(),
	    y2 = it->target().y(), l = pow( pow((new_node(0)-x1),2) + pow((new_node(1)-y1),2), 1/2);
	  if( ( ( pow((new_node(0)-x1),2) + pow((new_node(1)-y1),2) ) < pow(l*tol,2) )
	      || ( ( pow((new_node(0)-x2),2) + pow((new_node(1)-y2),2) ) < pow(l*tol,2) ) ){
	    is_on_mesh = true;
	    is_new = false;
	  }
	  if(is_new){
	    nodes.push_back(new_node);
	    new_nodes.getList().push_back(nb_node);
	    nb_node += 1 ;
	  }
	  if(!is_on_mesh){
	    new_node_per_elem(it->id(), 0) += 1;
	    new_node_per_elem(it->id(), ( 2*new_node_per_elem(it->id(),0)) - 1 ) = n ;
	    new_node_per_elem(it->id(), 2*new_node_per_elem(it->id(),0) ) = it->segId() ;
	  }
	}
      }
    }
  }
  const_cast<Mesh &>(this->mesh).sendEvent(new_nodes);

  AKANTU_DEBUG_OUT();
}

template<UInt dim, ElementType type>
void MeshSphereIntersector<dim, type>::computeIntersectionQueryList(const std::list<SK::Sphere_3> & query_list) {
  AKANTU_DEBUG_IN();
  
  parent_type::computeIntersectionQueryList(query_list);
  
  AKANTU_DEBUG_OUT();
}

#if defined(AKANTU_IGFEM)

template<UInt dim, ElementType type>
void MeshSphereIntersector<dim, type>::buildIgfemMesh(const std::list<SK::Sphere_3> & query_list) {
  AKANTU_DEBUG_IN();
  
  computeIntersectionQueryList(query_list);
  // Addition of the segment type in the mesh
  //TODO Mesh & mesh_no_const = const_cast<Mesh &>(this->mesh); or remove init const...
  const_cast<Mesh &>(this->mesh).addConnectivityType(_igfem_triangle_4, _not_ghost);
  const_cast<Mesh &>(this->mesh).addConnectivityType(_igfem_triangle_4, _ghost);
  const_cast<Mesh &>(this->mesh).addConnectivityType(_igfem_triangle_5, _not_ghost);
  const_cast<Mesh &>(this->mesh).addConnectivityType(_igfem_triangle_5, _ghost);

  Array<UInt> &
    connec_igfem_tri4 = const_cast<Array<UInt> &>(this->mesh.getConnectivity(_igfem_triangle_4)),
    & connec_igfem_tri5 = const_cast<Array<UInt> &>(this->mesh.getConnectivity(_igfem_triangle_5));
  Element element_tri4(_igfem_triangle_4, 0, _not_ghost, _ek_igfem),
    element_tri5(_igfem_triangle_5, 0, _not_ghost, _ek_igfem);
  NewElementsEvent new_elements;
  new_elements.getList().extendComponentsInterlaced(2, 1);

  Array<UInt> ctri3 = this->mesh.getConnectivity(_triangle_3);
  RemovedElementsEvent remove_elem(this->mesh);
  remove_elem.getNewNumbering().alloc(ctri3.getSize(), 1, _triangle_3, _not_ghost);
  Array<UInt> & new_numbering = remove_elem.getNewNumbering(_triangle_3, _not_ghost);
  UInt new_nb_tri3 = 0;

  for(UInt nel=0; nel != ctri3.getSize(); ++nel){
    if(new_node_per_elem(nel,0)!=0){
      Element element_tri3(_triangle_3, 0, _not_ghost);
      Array<Element> & new_elements_list = new_elements.getList() ;
      UInt n_new_el = new_elements_list.getSize();
      new_elements_list.resize(n_new_el+1);
      switch(new_node_per_elem(nel,0)){
      case 1 :{
	Vector<UInt> ctri4(4);
	switch(new_node_per_elem(nel,2)){
	case 1 :
	  ctri4(0) = ctri3(nel,2);
	  ctri4(1) = ctri3(nel,0);
	  ctri4(2) = ctri3(nel,1);
	  break;
	case 2 :
	  ctri4(0) = ctri3(nel,0);
	  ctri4(1) = ctri3(nel,1);
	  ctri4(2) = ctri3(nel,2);
	  break;
	case 3 :
	  ctri4(0) = ctri3(nel,1);
	  ctri4(1) = ctri3(nel,2);
	  ctri4(2) = ctri3(nel,0);
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
	new_elements_list(n_new_el,0) = element_tri4;
	break;
      }
      case 2 :{
	Vector<UInt> ctri5(5);
	if( (new_node_per_elem(nel,2)==1) && (new_node_per_elem(nel,4)==2) ){
	  ctri5(0) = ctri3(nel,1);
	  ctri5(1) = ctri3(nel,2);
	  ctri5(2) = ctri3(nel,0);
	  ctri5(3) = new_node_per_elem(nel,3);
	  ctri5(4) = new_node_per_elem(nel,1);
	}
	else if((new_node_per_elem(nel,2)==1) && (new_node_per_elem(nel,4)==3)){
	  ctri5(0) = ctri3(nel,0);
	  ctri5(1) = ctri3(nel,1);
	  ctri5(2) = ctri3(nel,2);
	  ctri5(3) = new_node_per_elem(nel,1);
	  ctri5(4) = new_node_per_elem(nel,3);
	}
	else if((new_node_per_elem(nel,2)==2) && (new_node_per_elem(nel,4)==3)){
	  ctri5(0) = ctri3(nel,2);
	  ctri5(1) = ctri3(nel,0);
	  ctri5(2) = ctri3(nel,1);
	  ctri5(3) = new_node_per_elem(nel,3);
	  ctri5(4) = new_node_per_elem(nel,1);
	}
	else if((new_node_per_elem(nel,2)==2) && (new_node_per_elem(nel,4)==1)){
	  ctri5(0) = ctri3(nel,1);
	  ctri5(1) = ctri3(nel,2);
	  ctri5(2) = ctri3(nel,0);
	  ctri5(3) = new_node_per_elem(nel,1);
	  ctri5(4) = new_node_per_elem(nel,3);
	}
	else if((new_node_per_elem(nel,2)==3) && (new_node_per_elem(nel,4)==1)){
	  ctri5(0) = ctri3(nel,0);
	  ctri5(1) = ctri3(nel,1);
	  ctri5(2) = ctri3(nel,2);
	  ctri5(3) = new_node_per_elem(nel,3);
	  ctri5(4) = new_node_per_elem(nel,1);
	}
	else if((new_node_per_elem(nel,2)==3) && (new_node_per_elem(nel,4)==2)){
	  ctri5(0) = ctri3(nel,2);
	  ctri5(1) = ctri3(nel,0);
	  ctri5(2) = ctri3(nel,1);
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
	new_elements_list(n_new_el,0) = element_tri5;
	break;
      }
      default:
	AKANTU_DEBUG_ERROR("Igfem cannot add "<< new_node_per_elem(nel,0) << " nodes to triangles");
	break;
      }
      element_tri3.element = nel;
      new_elements_list(n_new_el,1) = element_tri3;
      remove_elem.getList().push_back(element_tri3);
      new_numbering(nel) =  UInt(-1);
    }
    else{
      new_numbering(nel) = new_nb_tri3;
      ++new_nb_tri3 ;
    }
  }
  const_cast<Mesh &>(this->mesh).sendEvent(new_elements);
  const_cast<Mesh &>(this->mesh).sendEvent(remove_elem);

  AKANTU_DEBUG_OUT();
}

#endif

__END_AKANTU__

#endif // __AKANTU_MESH_SPHERE_INTERSECTOR_TMPL_HH__


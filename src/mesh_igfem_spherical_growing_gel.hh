/**
 * @file mesh_igfem_spherical_growing_gel.hh
 *
 * @author Clement Roux-Langlois <clement.roux@epfl.ch>
 *
 * @date creation: Mon Jul 13 2015
 *
 * @brief Computation of mesh intersection with sphere(s) and growing of these
 *        spheres. This class handle the intersectors templated for every element
 *        types.
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

//#if 0
#ifndef __AKANTU_MESH_IGFEM_SPHERICAL_GROWING_GEL_HH__
#define __AKANTU_MESH_IGFEM_SPHERICAL_GROWING_GEL_HH__

#include "aka_common.hh"
#include "mesh_sphere_intersector.hh"
#include "mesh_utils.hh"
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__


template<UInt dim>
class MeshIgfemSphericalGrowingGel {

  // definition of the element list
#define ELEMENT_LIST \
  (_triangle_3)		      \
  (_igfem_triangle_4)	      \
  (_igfem_triangle_5)

  // Solution 1
  /* #define INTERSECTOR_DEFINITION(type)			\
     MeshSphereIntersector<2, type> intersector##type(mesh);*/

  // Solution 2
#define INSTANTIATOR(_type)						\
  									\
  std::map<ElementType, MeshAbstractIntersector<SK::Sphere_3> *>::iterator iit = intersectors.find(_type); \
  if(iit != intersectors.end()) {					\
    inter = iit->second;						\
  } else {								\
    inter = new MeshSphereIntersector<dim, _type>(this->mesh);		\
    intersectors[_type] = inter;					\
  }

public:
  /// Construct from mesh
  MeshIgfemSphericalGrowingGel(Mesh & mesh):
    mesh(mesh),
    nb_nodes_fem(mesh.getNodes().getSize()),
    nb_enriched_nodes(0)
  {
    // Solution 1
    //   /// the mesh sphere intersector for the supported element type
    //   AKANTU_BOOST_APPLY_ON_LIST(INTERSECTOR_DEFINITION, ELEMENT_LIST)
    ElementTypeMapArray<UInt>::type_iterator tit = mesh.firstType(dim);
    ElementTypeMapArray<UInt>::type_iterator tend = mesh.lastType(dim);
    for(;tit != tend; ++tit) { // loop to add corresponding IGFEM element types
      if(*tit == _triangle_3){
        this->mesh.addConnectivityType(_igfem_triangle_4, _not_ghost);
	this->mesh.addConnectivityType(_igfem_triangle_4, _ghost);
        this->mesh.addConnectivityType(_igfem_triangle_5, _not_ghost);
	this->mesh.addConnectivityType(_igfem_triangle_5, _ghost);
      }	
      else
	AKANTU_DEBUG_ERROR("Not ready for mesh type " << *tit);
    }
    ElementTypeMapArray<UInt>::type_iterator tcit = mesh.getConnectivities().firstType(dim);
    ElementTypeMapArray<UInt>::type_iterator tcend = mesh.getConnectivities().lastType(dim);
    for(;tcit != tcend; ++tcit) {
      MeshAbstractIntersector<SK::Sphere_3> * inter;
      AKANTU_BOOST_LIST_SWITCH(INSTANTIATOR, ELEMENT_LIST, *tcit);
    }
  }

  /// Destructor
  ~MeshIgfemSphericalGrowingGel() {
    // Solution 2
    std::map<ElementType, MeshAbstractIntersector<SK::Sphere_3> *>::iterator iit
      = intersectors.begin();
    std::map<ElementType, MeshAbstractIntersector<SK::Sphere_3> *>::iterator iend
      = intersectors.end();
    for(;iit != iend; ++iit) {
      delete iit->second;
    }
  }

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// get the new_node_per_elem array
  // AKANTU_GET_MACRO(NewNodePerElem, new_node_per_elem, const Array<UInt>)

public:
  /// Construct the primitive tree object
  void constructData(){
    // Solution 1 TODO

    // Solution 2
    std::map<ElementType, MeshAbstractIntersector<SK::Sphere_3> *>::iterator iit
      = intersectors.begin();
    std::map<ElementType, MeshAbstractIntersector<SK::Sphere_3> *>::iterator iend
      = intersectors.end();
    for(;iit != iend; ++iit) {
      MeshAbstractIntersector<SK::Sphere_3> & intersector = *(iit->second);
      intersector.constructData();
    }
  }

  /// Remove the additionnal nodes
  void removeAdditionnalNodes(){
    AKANTU_DEBUG_IN();

    RemovedNodesEvent remove_nodes(this->mesh);
    Array<UInt> & nodes_removed = remove_nodes.getList();
    Array<UInt> & new_numbering = remove_nodes.getNewNumbering();
    UInt total_nodes = this->mesh.getNbNodes();
    UInt nb_new_enriched_nodes  = total_nodes - this->nb_enriched_nodes - this->nb_nodes_fem;
    UInt old_total_nodes = this->nb_nodes_fem + this->nb_enriched_nodes; 

    for(UInt nnod = 0; nnod < this->nb_nodes_fem ; ++nnod){
      new_numbering(nnod) = nnod ;
    }

    for(UInt nnod = nb_nodes_fem ; nnod < old_total_nodes; ++nnod){
      new_numbering(nnod) = UInt(-1) ;
      nodes_removed.push_back(nnod);
    }

    for(UInt nnod = 0; nnod < nb_new_enriched_nodes ; ++nnod){
      new_numbering(nnod + old_total_nodes) = this->nb_nodes_fem + nnod ;
    }
    for (UInt gt = _not_ghost; gt <= _ghost; ++gt) {
      GhostType ghost_type = (GhostType) gt;

      Mesh::type_iterator it  = mesh.firstType(_all_dimensions, ghost_type, _ek_not_defined);
      Mesh::type_iterator end = mesh.lastType(_all_dimensions, ghost_type, _ek_not_defined);
      for(; it != end; ++it) {

	ElementType type(*it);
	UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);

	const Array<UInt> & connectivity_vect = mesh.getConnectivity(type, ghost_type);
	UInt nb_element(connectivity_vect.getSize());
	UInt * connectivity = connectivity_vect.storage();

	UInt nb_nodes = nb_element*nb_nodes_per_element;

	for (UInt n = 0; n < nb_nodes; ++n, ++connectivity) {
	  UInt & node = *connectivity;
	  UInt old_node = node;
	  node = new_numbering(old_node);
	}	
      }
    }
    this->mesh.sendEvent(remove_nodes);
    AKANTU_DEBUG_OUT();
  }

  //Solution 1
  /* #define INTERSECTORS(type)					\
     intersector##type.buildResultFromQueryList(query_list)*/

  /// Build the IGFEM mesh
  void buildResultFromQueryList(const std::list<SK::Sphere_3> & query_list) {
    /// store number of currently enriched nodes
    this->nb_enriched_nodes = mesh.getNbNodes() - nb_nodes_fem;
    constructData();  
    //Solution 1
    //     it = mesh.firstType();
    //     end = mesh.lastType();
    //     for(;it != end;++it) {
    //       AKANTU_BOOST_LIST_SWITCH(INTERSECTORS, ELEMENT_LIST, *it)
    //     } }

    //Solution 2
    std::map<ElementType, MeshAbstractIntersector<SK::Sphere_3> *>::iterator iit
      = intersectors.begin();
    std::map<ElementType, MeshAbstractIntersector<SK::Sphere_3> *>::iterator iend
      = intersectors.end();
    for(;iit != iend; ++iit) {
      MeshAbstractIntersector<SK::Sphere_3> & intersector = *(iit->second);
      intersector.buildResultFromQueryList(query_list);
    }
    removeAdditionnalNodes();
    ///MeshUtils::purifyMesh(mesh);
  }

  /// increase sphere radius and build the IGFEM mesh
  void buildResultFromQueryList(const std::list<SK::Sphere_3> & query_list, Real inf_fact) {
    std::list<SK::Sphere_3>::const_iterator query_it = query_list.begin(),
      query_end = query_list.end();
    std::list<SK::Sphere_3> sphere_list;
    for (; query_it != query_end ; ++query_it) {
      SK::Sphere_3 sphere(query_it->center(),
			  query_it->squared_radius() * inf_fact * inf_fact );
      sphere_list.push_back(sphere);
    }
    buildResultFromQueryList(sphere_list);
  }

  /// Compute the intersection points between the mesh and the query list for all element types and send the NewNodeEvent
  void computeMeshQueryListIntersectionPoint(const std::list<SK::Sphere_3> & query_list) {
    /// store number of currently enriched nodes
    this->nb_enriched_nodes = mesh.getNbNodes() - nb_nodes_fem;
    Array<Real> new_nodes(0,dim); // new nodes from intersection of all types
    constructData();  
    
    std::map<ElementType, MeshAbstractIntersector<SK::Sphere_3> *>::iterator iit
      = intersectors.begin();
    std::map<ElementType, MeshAbstractIntersector<SK::Sphere_3> *>::iterator iend
      = intersectors.end();
    for(;iit != iend; ++iit) {
      MeshAbstractIntersector<SK::Sphere_3> & intersector = *(iit->second);
      intersector.computeMeshQueryListIntersectionPoint(query_list);
      const Array<Real> new_node_current_type = *(intersector.getNewNodes());
      const ElementTypeMapUInt & new_node_per_elem = intersector.getNewNodePerElem();
      
      /// Send the new node event
      if(new_node_current_type.getSize()) {
	Array<Real> & nodes = this->mesh.getNodes();
	UInt nb_node = nodes.getSize() ;
	NewIGFEMNodesEvent new_nodes_event;
	for(UInt in = 0; in < new_node_current_type.getSize(); ++in ) {
	  Vector<Real> new_node(dim, 0.0);
	  for(UInt id = 0; id < dim; ++id)
	    new_node(id) = new_node_current_type(in,id);
	  nodes.push_back(new_node);
	  new_nodes_event.getList().push_back(nb_node);
	  ++nb_node;
	}
	new_nodes_event.setNewNodePerElem(new_node_per_elem);
	new_nodes_event.setType(iit->first);
	this->mesh.sendEvent(new_nodes_event);
      }
    }
    
    //removeAdditionnalNodes();
    ///MeshUtils::purifyMesh(mesh);
  }

/// increase sphere radius and build the new intersection points between the mesh and the query list for all element types and send the NewNodeEvent
  void computeMeshQueryListIntersectionPoint(const std::list<SK::Sphere_3> & query_list, Real inf_fact) {
    std::list<SK::Sphere_3>::const_iterator query_it = query_list.begin(),
      query_end = query_list.end();
    std::list<SK::Sphere_3> sphere_list;
    for (; query_it != query_end ; ++query_it) {
      SK::Sphere_3 sphere(query_it->center(),
			  query_it->squared_radius() * inf_fact * inf_fact );
      sphere_list.push_back(sphere);
    }
    computeMeshQueryListIntersectionPoint(sphere_list);
  }

  /// Build the IGFEM mesh from intersection points
  void buildIgfemMesh(){
    AKANTU_DEBUG_IN();

    std::map<ElementType, MeshAbstractIntersector<SK::Sphere_3> *>::iterator iit
      = intersectors.begin();
    std::map<ElementType, MeshAbstractIntersector<SK::Sphere_3> *>::iterator iend
      = intersectors.end();
    for(;iit != iend; ++iit) {
      MeshAbstractIntersector<SK::Sphere_3> & intersector = *(iit->second);
      const ElementTypeMapUInt new_node_per_elem = intersector.getNewNodePerElem();
      const UInt nb_prim_by_el = intersector.getNbSegByEl();

      ElementType type = iit->first;
      if( (type !=_triangle_3) && (type !=_igfem_triangle_4) && (type !=_igfem_triangle_5) ) {
	AKANTU_DEBUG_ERROR("Not ready for mesh type " << type);
      }

      NewIGFEMElementsEvent new_elements;
      UInt total_new_elements = 0;
      Array<Element> & new_elements_list = new_elements.getList();
      Array<Element> & old_elements_list = new_elements.getOldElementsList();
      RemovedElementsEvent remove_elem(this->mesh);
      Int nb_proc = StaticCommunicator::getStaticCommunicator().getNbProc();

      for (ghost_type_t::iterator gt = ghost_type_t::begin();  gt != ghost_type_t::end(); ++gt) {
	GhostType ghost_type = *gt;

	if (nb_proc == 1 && ghost_type == _ghost) continue;
	if (!this->mesh.getConnectivities().exists(type, ghost_type)) continue;
	UInt n_new_el = 0;
	Array<UInt> & connec_type_tmpl = this->mesh.getConnectivity(type, ghost_type);

	Array<UInt>
	  & connec_igfem_tri4 = this->mesh.getConnectivity(_igfem_triangle_4, ghost_type),
	  & connec_igfem_tri5 = this->mesh.getConnectivity(_igfem_triangle_5, ghost_type),
	  & connec_tri3 = this->mesh.getConnectivity(_triangle_3, ghost_type);

	Element element_tri3(_triangle_3, 0, ghost_type, _ek_regular);
	Element element_tri4(_igfem_triangle_4, 0, ghost_type, _ek_igfem);
	Element element_tri5(_igfem_triangle_5, 0, ghost_type, _ek_igfem);


	remove_elem.getNewNumbering().alloc(connec_type_tmpl.getSize(), 1, type, ghost_type);
	Array<UInt> & new_numbering = remove_elem.getNewNumbering(type, ghost_type);
	UInt new_nb_type_tmpl = 0; // Meter of the number of element (type) will we loop on elements
	const Array<UInt> & new_node_per_elem_array = new_node_per_elem(type, ghost_type);

	for (UInt nel = 0 ; nel != new_node_per_elem_array.getSize(); ++nel) {
	  if( (type != _triangle_3)
	      && ( (new_node_per_elem_array(nel,0)==0)
		   || ( (new_node_per_elem_array(nel,0) == 1)
			&& ( ( (new_node_per_elem_array(nel,2)+1) % nb_prim_by_el )
			     != new_node_per_elem_array(nel, new_node_per_elem_array.getNbComponent() - 2) ) ) ) ){
	    Element element_type_tmpl(type, 0, ghost_type, Mesh::getKind(type));
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
	  else if( (new_node_per_elem_array(nel,0)!=0)
		   && !( (new_node_per_elem_array(nel,0) == 1)
			 && ( ( (new_node_per_elem_array(nel,2)+1) % nb_prim_by_el )
			      != new_node_per_elem_array(nel, new_node_per_elem_array.getNbComponent() - 2) ) ) ){
	    Element element_type_tmpl(type, 0, ghost_type);
	    element_type_tmpl.kind = Mesh::getKind(type);
	    new_elements_list.resize(n_new_el+1);
	    old_elements_list.resize(n_new_el+1);
	    switch(new_node_per_elem_array(nel,0)){
	    case 1 :{
	      Vector<UInt> ctri4(4);
	      switch(new_node_per_elem_array(nel,2)){
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
		AKANTU_DEBUG_ERROR("A triangle have only 3 segments not "<< new_node_per_elem_array(nel,2));
		break;
	      }
	      ctri4(3) = new_node_per_elem_array(nel,1);
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
	      if( (new_node_per_elem_array(nel,2)==1) && (new_node_per_elem_array(nel,4)==2) ){
		ctri5(0) = connec_type_tmpl(nel,1);
		ctri5(1) = connec_type_tmpl(nel,2);
		ctri5(2) = connec_type_tmpl(nel,0);
		ctri5(3) = new_node_per_elem_array(nel,3);
		ctri5(4) = new_node_per_elem_array(nel,1);
	      }
	      else if((new_node_per_elem_array(nel,2)==1) && (new_node_per_elem_array(nel,4)==3)){
		ctri5(0) = connec_type_tmpl(nel,0);
		ctri5(1) = connec_type_tmpl(nel,1);
		ctri5(2) = connec_type_tmpl(nel,2);
		ctri5(3) = new_node_per_elem_array(nel,1);
		ctri5(4) = new_node_per_elem_array(nel,3);
	      }
	      else if((new_node_per_elem_array(nel,2)==2) && (new_node_per_elem_array(nel,4)==3)){
		ctri5(0) = connec_type_tmpl(nel,2);
		ctri5(1) = connec_type_tmpl(nel,0);
		ctri5(2) = connec_type_tmpl(nel,1);
		ctri5(3) = new_node_per_elem_array(nel,3);
		ctri5(4) = new_node_per_elem_array(nel,1);
	      }
	      else if((new_node_per_elem_array(nel,2)==2) && (new_node_per_elem_array(nel,4)==1)){
		ctri5(0) = connec_type_tmpl(nel,1);
		ctri5(1) = connec_type_tmpl(nel,2);
		ctri5(2) = connec_type_tmpl(nel,0);
		ctri5(3) = new_node_per_elem_array(nel,1);
		ctri5(4) = new_node_per_elem_array(nel,3);
	      }
	      else if((new_node_per_elem_array(nel,2)==3) && (new_node_per_elem_array(nel,4)==1)){
		ctri5(0) = connec_type_tmpl(nel,0);
		ctri5(1) = connec_type_tmpl(nel,1);
		ctri5(2) = connec_type_tmpl(nel,2);
		ctri5(3) = new_node_per_elem_array(nel,3);
		ctri5(4) = new_node_per_elem_array(nel,1);
	      }
	      else if((new_node_per_elem_array(nel,2)==3) && (new_node_per_elem_array(nel,4)==2)){
		ctri5(0) = connec_type_tmpl(nel,2);
		ctri5(1) = connec_type_tmpl(nel,0);
		ctri5(2) = connec_type_tmpl(nel,1);
		ctri5(3) = new_node_per_elem_array(nel,1);
		ctri5(4) = new_node_per_elem_array(nel,3);
	      }
	      else{
		AKANTU_DEBUG_ERROR("A triangle have only 3 segments not "<< new_node_per_elem_array(nel,2) << "and" << new_node_per_elem_array(nel,4));
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
	      AKANTU_DEBUG_ERROR("Igfem cannot add "<< new_node_per_elem_array(nel,0) << " nodes to triangles");
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

	// for(UInt nel=new_node_per_elem_array.getSize(); nel < this->mesh.getNbElement(type); ++nel) {
	//     new_numbering(nel) = new_nb_type_tmpl;
	//     ++new_nb_type_tmpl;
	// }

	UInt el_index = 0;
	for (UInt e = 0; e < this->mesh.getNbElement(type, ghost_type); ++e) {
	  if (new_numbering(e) != UInt(-1)) {
	    new_numbering(e) = el_index;
	    ++el_index;
	  }
	}

	total_new_elements += n_new_el;
      }

      if(total_new_elements > 0){
	this->mesh.sendEvent(new_elements);
	this->mesh.sendEvent(remove_elem);
      }
    }

    removeAdditionnalNodes();
    AKANTU_DEBUG_OUT();
    }

  /// Build the IGFEM mesh from spheres
  void buildIgfemMeshFromSpheres(const std::list<SK::Sphere_3> & query_list){
    AKANTU_DEBUG_IN();
    
    computeMeshQueryListIntersectionPoint(query_list);
    buildIgfemMesh();

    AKANTU_DEBUG_OUT();
    }

  /// Build the IGFEM mesh from spheres with increase factor
  void buildIgfemMeshFromSpheres(const std::list<SK::Sphere_3> & query_list, Real inf_fact){
    AKANTU_DEBUG_IN();
    
    computeMeshQueryListIntersectionPoint(query_list, inf_fact);
    buildIgfemMesh();

    AKANTU_DEBUG_OUT();
    }

protected:
  /// Mesh used to construct the primitives
  Mesh & mesh;

  /// number of fem nodes in the initial mesh
  const UInt nb_nodes_fem;

  /// number of enriched nodes before intersection
  UInt nb_enriched_nodes;

  //Solution 2
  /// map of the elements types in the mesh and the corresponding intersectors
  std::map<ElementType, MeshAbstractIntersector<SK::Sphere_3> *> intersectors;
};
 
__END_AKANTU__

#endif // __AKANTU_MESH_IGFEM_SPHERICAL_GROWING_GEL_HH__

//#endif //

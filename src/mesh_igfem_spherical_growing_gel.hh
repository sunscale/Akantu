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
    nb_nodes_fem(mesh.getNodes().getSize())
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
    UInt nnod = 0;
    for(; nnod != nb_nodes_fem ; ++nnod){
      new_numbering(nnod) = nnod ;
    }
    for(nnod = nb_nodes_fem ; nnod != this->mesh.getNodes().getSize(); ++nnod){
      new_numbering(nnod) = UInt(-1) ;
      nodes_removed.push_back(nnod);
    }
    this->mesh.sendEvent(remove_nodes);
    AKANTU_DEBUG_OUT();
  }

  //Solution 1
  /* #define INTERSECTORS(type)					\
     intersector##type.buildResultFromQueryList(query_list)*/

  /// Build the IGFEM mesh
  void buildResultFromQueryList(const std::list<SK::Sphere_3> & query_list) {
    constructData();
    removeAdditionnalNodes();

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

protected:
  /// Mesh used to construct the primitives
  Mesh & mesh;

  /// number of fem nodes in the initial mesh
  const UInt nb_nodes_fem;

  //Solution 2
  /// map of the elements types in the mesh and the corresponding intersectors
  std::map<ElementType, MeshAbstractIntersector<SK::Sphere_3> *> intersectors;
};
 
__END_AKANTU__

#endif // __AKANTU_MESH_IGFEM_SPHERICAL_GROWING_GEL_HH__

//#endif //

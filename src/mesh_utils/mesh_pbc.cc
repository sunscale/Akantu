/**
 * @file   mesh_pbc.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author David Simon Kammer <david.kammer@epfl.ch>
 *
 * @date creation: Wed Feb 09 2011
 * @date last modification: Fri Aug 01 2014
 *
 * @brief  periodic boundary condition connectivity tweak
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include <map>
/* -------------------------------------------------------------------------- */
#include "mesh_utils.hh"
#include "element_group.hh"
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
/// class that sorts a set of nodes of same coordinates in 'dir' direction
class CoordinatesComparison {
public:
  CoordinatesComparison (const UInt dimension,
			 const UInt dirx, const UInt diry,
			 Real normalization,
			 Real tolerance,
			 Real * coords):
    dim(dimension),dir_x(dirx),dir_y(diry),normalization(normalization),
    tolerance(tolerance),coordinates(coords){}
  // answers the question whether n2 is larger or equal to n1
  bool operator() (UInt n1, UInt n2){
    Real p1_x = coordinates[dim*n1+dir_x];
    Real p2_x = coordinates[dim*n2+dir_x];
    Real diff_x = p1_x - p2_x;
    if (dim == 2 || std::abs(diff_x)/normalization > tolerance)
      return diff_x > 0.0 ? false : true;
    else if (dim > 2){
      Real p1_y = coordinates[dim*n1+dir_y];
      Real p2_y = coordinates[dim*n2+dir_y];
      Real diff_y = p1_y - p2_y;
      return diff_y > 0 ? false : true;
    }
    return true;
  }
private:
  UInt dim;
  UInt dir_x;
  UInt dir_y;
  Real normalization;
  Real tolerance;
  Real * coordinates;
};

/* -------------------------------------------------------------------------- */
// void MeshUtils::tweakConnectivityForPBC(Mesh & mesh,
//					bool flag_x,
//					bool flag_y,
//					bool flag_z){
//   std::map<UInt,UInt> pbc_pair;
//   mesh.computeBoundingBox();
//   mesh.pbc_directions[0] = flag_x;
//   mesh.pbc_directions[1] = flag_y;
//   mesh.pbc_directions[2] = flag_z;

//   if (flag_x) computePBCMap(mesh,0,pbc_pair);
//   if (flag_y) computePBCMap(mesh,1,pbc_pair);
//   if (flag_z) computePBCMap(mesh,2,pbc_pair);

//   {
//     std::map<UInt,UInt>::iterator it = pbc_pair.begin();
//     std::map<UInt,UInt>::iterator end = pbc_pair.end();

//     Real * coords = mesh.nodes->values;
//     UInt dim = mesh.getSpatialDimension();
//     while(it != end){
//       UInt i1 = (*it).first;
//       UInt i2 = (*it).second;

//       AKANTU_DEBUG_INFO("pairing " << i1 << "("
//			<< coords[dim*i1] << "," << coords[dim*i1+1] << ","
//			<< coords[dim*i1+2]
//			<< ") with"
//			<< i2 << "("
//			<< coords[dim*i2] << "," << coords[dim*i2+1] << ","
//			<< coords[dim*i2+2]
//			<< ")");
//       ++it;
//     }
//   }

//   //allocate and initialize list of reversed elements
//   mesh.initElementTypeMapArray<UInt>Array(mesh.reversed_elements_pbc,1,0,mesh.id,"reversed");
//   // now loop over the elements to change the connectivity of some elements
//   const Mesh::ConnectivityTypeList & type_list = mesh.getConnectivityTypeList();
//   Mesh::ConnectivityTypeList::const_iterator it;
//   for(it = type_list.begin(); it != type_list.end(); ++it) {
//     ElementType type = *it;
//     UInt nb_elem = mesh.getNbElement(type);
//     UInt nb_nodes_per_elem = mesh.getNbNodesPerElement(type);
//     UInt * conn = mesh.getConnectivityPointer(type)->values;
//     UInt index = 0;
//     Array<UInt> & list = *(mesh.reversed_elements_pbc[type]);
//     for (UInt el = 0; el < nb_elem; el++) {
//       UInt flag_should_register_elem = false;
//       for (UInt k = 0; k < nb_nodes_per_elem; ++k,++index){
//	if (pbc_pair.count(conn[index])){
//	  flag_should_register_elem = true;
//	  AKANTU_DEBUG_INFO("elem list size " << list.getSize());
//	  AKANTU_DEBUG_INFO("node " << conn[index] +1
//			    << " switch to "
//			    << pbc_pair[conn[index]]+1);
//	  // for (UInt toto = 0; toto < 3; ++toto) {
//	  //   AKANTU_DEBUG_INFO("dir " << toto << " coords "
//	  //		      << mesh.nodes->values[conn[index]*3+toto]
//	  //		      << " switch to "
//	  //		      << mesh.nodes->values[pbc_pair[conn[index]]*3+toto]);
//	  // }
//	  std::stringstream str_temp;
//	  str_temp << "initial elem(" << el << ") is ";
//	  for (UInt l = 0 ; l < nb_nodes_per_elem ; ++ l){
//	    str_temp << conn[el*nb_nodes_per_elem+l]+1 << " ";
//	  }
//	  AKANTU_DEBUG_INFO(str_temp.str());
//	  conn[index] = pbc_pair[conn[index]];
//	}
//       }
//       if (flag_should_register_elem) list.push_back(el);
//     }
//   }
// }

/* -------------------------------------------------------------------------- */
void MeshUtils::computePBCMap(const Mesh & mymesh,
			      const UInt dir,
			      std::map<UInt,UInt> & pbc_pair){
  Array<UInt> selected_left;
  Array<UInt> selected_right;

  Real * coords = mymesh.nodes->storage();
  const UInt nb_nodes = mymesh.nodes->getSize();
  const UInt dim = mymesh.getSpatialDimension();

  if (dim <= dir) return;

  AKANTU_DEBUG_INFO("min " << mymesh.lower_bounds[dir]);
  AKANTU_DEBUG_INFO("max " << mymesh.upper_bounds[dir]);

  for (UInt i = 0; i < nb_nodes; ++i) {
    AKANTU_DEBUG_TRACE("treating " << coords[dim*i+dir]);
    if(Math::are_float_equal(coords[dim*i+dir], mymesh.lower_bounds[dir])){
      AKANTU_DEBUG_TRACE("pushing node " << i << " on the left side");
      selected_left.push_back(i);
    }
    else if(Math::are_float_equal(coords[dim*i+dir], mymesh.upper_bounds[dir])){
      selected_right.push_back(i);
      AKANTU_DEBUG_TRACE("pushing node " << i << " on the right side");
    }
  }

  AKANTU_DEBUG_INFO("found " << selected_left.getSize() << " and " << selected_right.getSize()
	       << " nodes at each boundary for direction " << dir);

  // match pairs
  MeshUtils::matchPBCPairs(mymesh, dir, selected_left, selected_right, pbc_pair);

}

/* -------------------------------------------------------------------------- */
void MeshUtils::computePBCMap(const Mesh & mymesh,
			      const SurfacePair & surface_pair,
			      std::map<UInt,UInt> & pbc_pair) {

  Array<UInt> selected_first;
  Array<UInt> selected_second;

  // find nodes on surfaces
  const ElementGroup & first_surf = mymesh.getElementGroup(surface_pair.first);
  const ElementGroup & second_surf = mymesh.getElementGroup(surface_pair.second);

  // if this surface pair is not on this proc
  if (first_surf.getNbNodes() == 0 || second_surf.getNbNodes() == 0) {
    AKANTU_DEBUG_WARNING("computePBCMap has at least one surface without any nodes. I will ignore it.");
    return;
  }

  // copy nodes from element group
  selected_first.copy(first_surf.getNodeGroup().getNodes());
  selected_second.copy(second_surf.getNodeGroup().getNodes());

  // coordinates
  const Array<Real> & coords = mymesh.getNodes();
  const UInt dim = mymesh.getSpatialDimension();

  // variables to find min and max of surfaces
  Real first_max[3], first_min[3];
  Real second_max[3], second_min[3];
  for (UInt i=0; i<dim; ++i) {
    first_min[i] = std::numeric_limits<Real>::max();
    second_min[i] = std::numeric_limits<Real>::max();
    first_max[i] = -std::numeric_limits<Real>::max();
    second_max[i] = -std::numeric_limits<Real>::max();
  }

  // find min and max of surface nodes
  for (Array<UInt>::scalar_iterator it = selected_first.begin();
       it != selected_first.end();
       ++it) {
    for (UInt i=0; i<dim; ++i) {
      if (first_min[i] > coords(*it,i)) first_min[i] = coords(*it,i);
      if (first_max[i] < coords(*it,i)) first_max[i] = coords(*it,i);
    }
  }
  for (Array<UInt>::scalar_iterator it = selected_second.begin();
       it != selected_second.end();
       ++it) {
    for (UInt i=0; i<dim; ++i) {
      if (second_min[i] > coords(*it,i)) second_min[i] = coords(*it,i);
      if (second_max[i] < coords(*it,i)) second_max[i] = coords(*it,i);
    }
  }

  // find direction of pbc
  Int first_dir = -1;
#ifndef AKANTU_NDEBUG
  Int second_dir = -2;
#endif
  for (UInt i=0; i<dim; ++i) {
    if (Math::are_float_equal(first_min[i], first_max[i])) {
      first_dir = i;
    }
#ifndef AKANTU_NDEBUG
    if (Math::are_float_equal(second_min[i], second_max[i])) {
      second_dir = i;
    }
#endif
  }

  AKANTU_DEBUG_ASSERT(first_dir == second_dir, "Surface pair has not same direction. Surface "
		      << surface_pair.first << " dir=" << first_dir << " ; Surface "
		      << surface_pair.second << " dir=" << second_dir);
  UInt dir = first_dir;

  // match pairs
  if (first_min[dir] < second_min[dir])
    MeshUtils::matchPBCPairs(mymesh, dir, selected_first, selected_second, pbc_pair);
  else
    MeshUtils::matchPBCPairs(mymesh, dir, selected_second, selected_first, pbc_pair);
}


/* -------------------------------------------------------------------------- */
void MeshUtils::matchPBCPairs(const Mesh & mymesh,
			      const UInt dir,
			      Array<UInt> & selected_left,
			      Array<UInt> & selected_right,
			      std::map<UInt,UInt> & pbc_pair) {


  // tolerance is that large because most meshers generate points coordinates 
  // with single precision only (it is the case of GMSH for instance)
  Real tolerance = 1e-7;
  Real * coords = mymesh.nodes->storage();
  const UInt dim = mymesh.getSpatialDimension();
  Real normalization = mymesh.upper_bounds[dir]-mymesh.lower_bounds[dir];

  AKANTU_DEBUG_ASSERT(std::abs(normalization) > Math::getTolerance(), 
		      "In matchPBCPairs: The normalization is zero. "
		      << "Did you compute the bounding box of the mesh?");

  UInt dir_x = UInt(-1) ,dir_y = UInt(-1);

  if (dim == 3){
    if (dir == 0){
      dir_x = 1;dir_y = 2;
    }
    else if (dir == 1){
      dir_x = 0;dir_y = 2;
    }
    else if (dir == 2){
      dir_x = 0;dir_y = 1;
    }
  }
  else if (dim == 2){
    if (dir == 0){
      dir_x = 1;
    }
    else if (dir == 1){
      dir_x = 0;
    }
  }

  CoordinatesComparison compare_nodes(dim,dir_x,dir_y,normalization,tolerance,coords);

  std::sort(selected_left.begin(),selected_left.end(),compare_nodes);
  std::sort(selected_right.begin(),selected_right.end(),compare_nodes);


  Array<UInt>::scalar_iterator it_left = selected_left.begin();
  Array<UInt>::scalar_iterator end_left = selected_left.end();

  Array<UInt>::scalar_iterator it_right = selected_right.begin();
  Array<UInt>::scalar_iterator end_right = selected_right.end();

  while ((it_left != end_left) && (it_right != end_right) ){
    UInt i1 = *it_left;
    UInt i2 = *it_right;

    AKANTU_DEBUG_TRACE("do I pair? " << i1 << "("
		      << coords[dim*i1] << "," << coords[dim*i1+1] << ","
		       << coords[dim*i1+2]
		      << ") with"
		      << i2 << "("
		      << coords[dim*i2] << "," << coords[dim*i2+1] << ","
		       << coords[dim*i2+2]
		      << ") in direction " << dir);


    Real dx = 0.0;
    Real dy = 0.0;
    if (dim >= 2) dx = coords[dim*i1 + dir_x] - coords[dim*i2 + dir_x];
    if (dim == 3) dy = coords[dim*i1 + dir_y] - coords[dim*i2 + dir_y];

    if (fabs(dx*dx+dy*dy)/normalization < tolerance) {
      //then i match these pairs
      if (pbc_pair.count(i2)){
	i2 = pbc_pair[i2];
      }
      pbc_pair[i1] = i2;

      AKANTU_DEBUG_TRACE("pairing " << i1 << "("
			 << coords[dim*i1] << "," << coords[dim*i1+1] << ","
			 << coords[dim*i1+2]
			 << ") with"
			 << i2 << "("
			 << coords[dim*i2] << "," << coords[dim*i2+1] << ","
			 << coords[dim*i2+2]
			 << ") in direction " << dir);
      ++it_left;
      ++it_right;
    } else if (compare_nodes(i1, i2)) {
      ++it_left;
    } else {
      ++it_right;
    }

  }
  AKANTU_DEBUG_INFO("found " <<  pbc_pair.size() << " pairs for direction " << dir);

}

__END_AKANTU__

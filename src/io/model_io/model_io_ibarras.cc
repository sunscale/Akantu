/**
 * @file   model_io_ibarras.cc
 *
 * @author Fabian Barras <fabian.barras@epfl.ch>
 *
 * @date creation: Wed Jan 16 2013
 * @date last modification: Thu Jun 05 2014
 *
 * @brief  Mesh Reader specially created for Wood's Tower analysis performed by the Institute I-Barras
 *
 * @section LICENSE
 *
 * Copyright (©) 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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


/* -------------------------------------------------------------------------- */


/* -------------------------------------------------------------------------- */
#include <fstream>
#include <iostream>

/* -------------------------------------------------------------------------- */
#include "model_io_ibarras.hh"
#include "static_communicator.hh"
#include "aka_math.hh"
#include "static_communicator.hh"
#include "sparse_matrix.hh"
#include "solver.hh"
#include "integration_scheme_2nd_order.hh"
/* -------------------------------------------------------------------------- */
#include <string.h>
/* -------------------------------------------------------------------------- */
#include <stdio.h>

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
/*   Methods Implentations                                                    */
/* -------------------------------------------------------------------------- */

void ModelIOIBarras::read(const std::string & filename, Model & mod) {

  if (! dynamic_cast<StructuralMechanicsModel*>(&mod)){
    AKANTU_DEBUG_ERROR("");
  }

  StructuralMechanicsModel & model = static_cast<StructuralMechanicsModel&>(mod);
  Mesh * mesh = new Mesh(3);

  std::ifstream infile;
  infile.open(filename.c_str());

  std::string line;
  UInt current_line = 0;


  if(!infile.good()) {
    AKANTU_DEBUG_ERROR("Cannot open file " << filename);
  }

  // Get Nodes Position
  std::getline(infile, line);
  current_line++;
  std::stringstream sstr(line);
  UInt nb_nodes;
  sstr >> nb_nodes;

  UInt spatial_dimension = 3;
  Real coord[spatial_dimension];
  Real * temp_nodes = new Real[nb_nodes*spatial_dimension];
  UInt * connect_to_akantu = new UInt[nb_nodes];
  std::fill_n(connect_to_akantu ,0, nb_nodes);
  std::fill_n(temp_nodes ,0., nb_nodes * spatial_dimension);
  for (UInt i = 0; i < nb_nodes; ++i) {
    UInt offset = i * spatial_dimension;

    std::getline(infile, line);
    std::stringstream sstr_node(line);
    sstr_node >> coord[0] >> coord[1] >> coord[2];
    current_line++;

    /// read the coordinates of structural nodes and help nodes
    for(UInt j = 0; j < spatial_dimension; ++j)
      temp_nodes[offset + j] = coord[j];
      }

  // Get Connectivities
  std::getline(infile, line);
  current_line++;
  std::stringstream sstr_elem(line);
  UInt nb_elements;
  sstr_elem >> nb_elements;

  mesh->addConnectivityType(_bernoulli_beam_3);
  Array<UInt> & connectivity = const_cast<Array<UInt> &>(mesh->getConnectivity(_bernoulli_beam_3));

  connectivity.resize(nb_elements);

  UInt nonodes[2];
  UInt nb_struct_nodes=1;
  for (UInt i = 0; i < nb_elements; ++i){

    std::getline(infile, line);
    std::stringstream sstr_element(line);
    sstr_element >> nonodes[0] >> nonodes[1];
    current_line++;


    /// read the connectivities
    for(UInt j = 0; j < 2; ++j){

      if (connect_to_akantu[nonodes[j]-1]==0){
	  connect_to_akantu[nonodes[j]-1]=nb_struct_nodes;
	  ++nb_struct_nodes;
      }
      connectivity(i,j)=connect_to_akantu[nonodes[j]-1]-1;
    }
  }
  nb_struct_nodes-=1;

  /// read the coordinates of structural nodes
  Array<Real> & nodes = const_cast<Array<Real> &>(mesh->getNodes());
  nodes.resize(nb_struct_nodes);

  for(UInt k = 0; k < nb_nodes; ++k){
    if (connect_to_akantu[k]!=0){
      for(UInt j = 0; j < spatial_dimension; ++j)
	nodes(connect_to_akantu[k]-1,j) = temp_nodes[k*spatial_dimension+j];
    }
  }

  //MeshPartitionScotch partition(*mesh, spatial_dimension);
  //partition.reorder();


  ///Apply Boundaries
  model.registerFEEngineObject<StructuralMechanicsModel::MyFEEngineType>("StructuralMechanicsModel", *mesh, spatial_dimension);

  model.initModel();

  model.initArrays();

  Array<bool> & blocked_dofs = model.getBlockedDOFs();

  std::getline(infile, line);
  std::stringstream sstr_nb_boundaries(line);
  UInt nb_boundaries;
  sstr_nb_boundaries >> nb_boundaries;
  current_line++;

  for (UInt i = 0; i < nb_boundaries; ++i){
    std::getline(infile, line);
    std::stringstream sstr_boundary(line);
    UInt boundnary_node;
    sstr_boundary >> boundnary_node;
    current_line++;

    for (UInt j = 0; j < spatial_dimension; ++j)
    blocked_dofs(connect_to_akantu[boundnary_node-1]-1 ,j)=true;
  }

  ///Define Materials

  std::getline(infile, line);
  std::stringstream sstr_nb_materials(line);
  UInt nb_materials;
  sstr_nb_materials >> nb_materials;
  current_line++;

  for (UInt i = 0; i < nb_materials; ++i){

    std::getline(infile, line);
    std::stringstream sstr_material(line);
    Real material[6];
    sstr_material >> material[0] >> material[1] >> material[2] >> material[3] >> material[4] >> material[5];
    current_line++;

    StructuralMaterial mat;
    mat.E = material[0];
    mat.GJ = material[1] * material[2];
    mat.Iy = material[3];
    mat.Iz = material[4];
    mat.A = material[5];

    model.addMaterial(mat);
  }

  /// Apply normals and Material TO IMPLEMENT

  UInt property[2];

  Array<UInt> & element_material = model.getElementMaterial(_bernoulli_beam_3);

  mesh->initNormals();
  Array<Real> & normals = const_cast<Array<Real> &>(mesh->getNormals(_bernoulli_beam_3));
  normals.resize(nb_elements);

  for (UInt i = 0; i < nb_elements; ++i){

    std::getline(infile, line);
    std::stringstream sstr_properties(line);
    sstr_properties >> property[0] >> property[1];
    current_line++;


    /// Assign material
    element_material(i)=property[0]-1;

    /// Compute normals

    Real x[3];
    Real v[3];
    Real w[3];
    Real n[3];

    if (property[1]==0){
      for (UInt j = 0; j < spatial_dimension; ++j){
	x[j] = nodes(connectivity(i,1),j) - nodes(connectivity(i,0),j);
      }
      n[0] = x[1];
      n[1] = -x[0];
      n[2] = 0.;

     }

    else{
      for (UInt j = 0; j < spatial_dimension; ++j){
	x[j] = nodes(connectivity(i,1),j) - nodes(connectivity(i,0),j);
	v[j] = nodes(connectivity(i,1),j) - temp_nodes[(property[1]-1) * spatial_dimension + j];
	Math::vectorProduct3(x, v, w);
	Math::vectorProduct3(x, w, n);

      }
    }

    Math::normalize3(n);
    for (UInt j = 0; j < spatial_dimension; ++j){
      normals(i,j) = n[j];
      }
  }

  model.computeRotationMatrix(_bernoulli_beam_3);
  infile.close();

}
/* -------------------------------------------------------------------------- */
void ModelIOIBarras::assign_sets(const std::string & filename, StructuralMechanicsModel & model){

  std::ifstream infile;
  infile.open(filename.c_str());

  std::string line;
  UInt current_line = 0;


  if(!infile.good()) {
    AKANTU_DEBUG_ERROR("Cannot open file " << filename);
  }

  // Define Sets of Beams

  Array<UInt> & set_ID = model.getSet_ID(_bernoulli_beam_3);
  set_ID.clear();

  std::getline(infile, line);
  std::stringstream sstr_nb_sets(line);
  UInt nb_sets;
  sstr_nb_sets >> nb_sets;
  current_line++;

  UInt no_element[2];

  for (UInt i = 0; i < nb_sets; ++i){
    std::getline(infile, line);
    std::stringstream sstr_set(line);
    sstr_set >> no_element[0];
    no_element[1]=no_element[0];
    sstr_set >> no_element[1];

    while (no_element[0]!=0) {

      for (UInt j = no_element[0]-1 ; j < no_element[1] ; ++j){
	set_ID(j) = i+1;
      }
      std::getline(infile, line);
      std::stringstream sstr_sets(line);
      sstr_sets >> no_element[0];
      no_element[1]=no_element[0];
      sstr_sets >> no_element[1];
    }
  }
}

__END_AKANTU__

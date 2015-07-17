/**
 * @file   test_solid_mechanics_model_igfem.cc
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 *
 *
 * @brief  test the solidmechancis model for IGFEM analysis
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
// #include "mesh_segment_intersector.hh"
// #include "mesh_sphere_intersector.hh"
// #include "geom_helper_functions.hh"
// #include "mesh_geom_common.hh"
#include "solid_mechanics_model.hh"
/* -------------------------------------------------------------------------- */
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <cmath>
#include <math.h> 
#include "dumper_paraview.hh"
/* -------------------------------------------------------------------------- */
using namespace akantu;
int main(int argc, char *argv[]) {

  initialize("material.dat", argc, argv);

  /// problem dimension
  UInt spatial_dimension = 2;  

  /// mesh creation
  Mesh mesh(spatial_dimension);
  mesh.read(std::string(argv[1]));

  /// geometry of inclusion
  Real radius_inclusion = 0.36;
  /// model creation
  SolidMechanicsModel model(mesh);
  MeshDataMaterialSelector<UInt> * mat_selector;
  mat_selector = new MeshDataMaterialSelector<UInt>("tag_1", model, 5);
 
  model.setMaterialSelector(*mat_selector);
  model.initFull(SolidMechanicsModelOptions(_static));

  /// boundary conditions
  mesh.computeBoundingBox();
  const Vector<Real> & lowerBounds = mesh.getLowerBounds();
  const Vector<Real> & upperBounds = mesh.getUpperBounds();
  Real bottom  = lowerBounds(1);
  Real top = upperBounds(1);
  Real left = lowerBounds(0);
  Real right = upperBounds(0);

  Real eps = std::abs((top - bottom) * 1e-12);
  const Array<Real> & pos = mesh.getNodes();
  Array<Real> & disp = model.getDisplacement();
  Array<bool> & boun = model.getBlockedDOFs();
  Real radius = 0;
  Real phi = 0;
  Real alpha = 0.929131442935370;
  
  /// absolute confinement
  for (UInt i = 0; i < mesh.getNbNodes(); ++i) {
 
    if(std::abs(pos(i,0) - left) < eps) {
      radius = std::sqrt(pos(i,0)*pos(i,0) + pos(i,1)*pos(i,1));
      phi = std::atan2(pos(i,1), pos(i,0));
      boun(i,0) = true;
      disp(i,0) = cos(phi) * ( (radius - 4./radius) * alpha + 4./radius );
      boun(i,1) = true;
      disp(i,1) = sin(phi) * ( (radius - 4./radius) * alpha + 4./radius );
    }

    if(std::abs(pos(i,0) - right) < eps) {
      radius = std::sqrt(pos(i,0)*pos(i,0) + pos(i,1)*pos(i,1));
      phi = std::atan2(pos(i,1), pos(i,0));
      boun(i,0) = true;
      disp(i,0) = cos(phi) * ( (radius - 4./radius) * alpha + 4./radius );
      boun(i,1) = true;
      disp(i,1) = sin(phi) * ( (radius - 4./radius) * alpha + 4./radius );
    }

    if(std::abs(pos(i,1) - top) < eps) {
      radius = std::sqrt(pos(i,0)*pos(i,0) + pos(i,1)*pos(i,1));
      phi = std::atan2(pos(i,1), pos(i,0));
      boun(i,0) = true;
      disp(i,0) = cos(phi) * ( (radius - 4./radius) * alpha + 4./radius );
      boun(i,1) = true;
      disp(i,1) = sin(phi) * ( (radius - 4./radius) * alpha + 4./radius );
    }

    if(std::abs(pos(i,1) - bottom) < eps) {
      radius = std::sqrt(pos(i,0)*pos(i,0) + pos(i,1)*pos(i,1));
      phi = std::atan2(pos(i,1), pos(i,0));
      boun(i,0) = true;
      disp(i,0) = cos(phi) * ( (radius - 4./radius) * alpha + 4./radius );
      boun(i,1) = true;
      disp(i,1) = sin(phi) * ( (radius - 4./radius) * alpha + 4./radius );
    }

  }


  model.setBaseName("regular_mesh_test");
  model.addDumpField("material_index");
  model.addDumpFieldVector("displacement");
  model.addDumpField("blocked_dofs");
  model.addDumpField("stress");
  model.addDumpField("damage");
  model.addDumpField("Sc");

  Math::setTolerance(1e-14);
  //  dumper_igfem.dump();

  model.dump();


  bool factorize = false;
  bool converged = false;
  Real error; 
  converged = model.solveStep<_scm_newton_raphson_tangent, _scc_increment>(1e-4, error, 2, factorize);
  model.dump();


  finalize();
  return EXIT_SUCCESS;
}

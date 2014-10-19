/**
 * @file   test_material_stiffness_proportional_damping.cc
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 *
 * @date   Mon Nov 14 09:15:53 2011
 *
 * @brief  test of the material elastic caughey - physical aspecte
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include <iostream>

/* -------------------------------------------------------------------------- */
#include <limits>
#include <fstream>

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "mesh.hh"
#include "mesh_io.hh"
#include "mesh_io_msh.hh"
#include "solid_mechanics_model.hh"
#include "material.hh"

/* -------------------------------------------------------------------------- */
#ifdef AKANTU_USE_IOHELPER
#  include "io_helper_tools.hh"
#endif //AKANTU_USE_IOHELPER

using namespace akantu;

static bool testFloat(Real a, Real b, Real adm_error);

int main(int argc, char *argv[])
{
  akantu::initialize(argc, argv);
  akantu::debug::setDebugLevel(akantu::dblWarning);
  
  const ElementType element_type = TYPE;
  const UInt dim = ElementClass<TYPE>::getSpatialDimension();
  
  /// load mesh
  Mesh mesh(dim);
  MeshIOMSH mesh_io;
  std::stringstream meshname_sstr; 
  meshname_sstr << "single_" << element_type << ".msh";
  mesh_io.read(meshname_sstr.str().c_str(), mesh);
  
  UInt max_steps = 1000;
  Real time_factor = 0.1;
  UInt nb_nodes = mesh.getNbNodes();
    
  SolidMechanicsModel model(mesh);
  
  /* ------------------------------------------------------------------------ */
  /* Initialization                                                           */
  /* ------------------------------------------------------------------------ */
  model.initArrays();
  model.getForce().clear();
  model.getVelocity().clear();
  model.getAcceleration().clear();
  model.getDisplacement().clear();
  model.updateResidual();
  
  model.initExplicit();
  model.initModel();
  model.readMaterials("material_elastic_caughey_damping.dat");
  Material & my_mat = model.getMaterial(0);
  Real a_value = 1e-6;//3.836e-06;
  my_mat.setProperty("alpha", a_value);
  model.initMaterials();
  
  std::cout << model.getMaterial(0) << std::endl;
  
  model.assembleMassLumped();
  
  /* ------------------------------------------------------------------------ */
  /* Boundary + initial conditions                                            */
  /* ------------------------------------------------------------------------ */
  Array<UInt> the_nodes(0,1);
  Real imposed_disp = 0.1;
  for (UInt i = 0; i < nb_nodes; ++i) {
    // block lower nodes
    if(mesh.getNodes().storage()[i*dim+1] < 0.5) {
      for (UInt j=0; j<dim; ++j) 
	model.getBlockedDOFs().storage()[dim*i + j] = true;
    }
    // impose displacement
    else {
      model.getBlockedDOFs().storage()[dim*i + 0] = true;
      model.getDisplacement().storage()[dim*i + 1] = imposed_disp;
      the_nodes.push_back(i);
    }
  }
  model.updateResidual();

#ifdef AKANTU_USE_IOHELPER
  iohelper::DumperParaview dumper;
  paraviewInit(dumper, model, element_type, "test_mat_el_cau_dump");
#endif //AKANTU_USE_IOHELPER

  /// Setting time step
  Real time_step = model.getStableTimeStep() * time_factor;
  std::cout << "Time Step = " << time_step << "s" << std::endl;
  model.setTimeStep(time_step);

  /* ------------------------------------------------------------------------ */
  /* Main loop                                                                */
  /* ------------------------------------------------------------------------ */
  for(UInt s = 1; s <= max_steps; ++s) {

    if (s%100 == 0)
      std::cout << "passing step " << s << "/" << max_steps << std::endl;

    model.explicitPred();
    model.updateResidual();
    model.updateAcceleration();
    model.explicitCorr();

#ifdef AKANTU_USE_IOHELPER
    if(s % 100 == 0) dumper.Dump();
#endif //AKANTU_USE_IOHELPER
  }

  /*
  for (UInt i=0; i<the_nodes.getSize(); ++i) {
    std::cout << "disp " <<  model.getDisplacement().storage()[the_nodes(i)*dim+1] << "; vel " << model.getVelocity().storage()[the_nodes(i)*dim+1] << std::endl;
  }
  */

  /* ------------------------------------------------------------------------ */
  /* Test solution                                                            */
  /* ------------------------------------------------------------------------ */
  Real disp_tol = 1e-07;
  Real velo_tol = 1e-03;
  // solution triangle_3
  Array<Real> disp_triangle_3(0,1); disp_triangle_3.push_back(-0.0344941);
  Array<Real> velo_triangle_3(0,1); velo_triangle_3.push_back(-433.9);
  // solution quadrangle_4
  Array<Real> disp_quadrangle_4(0,1); 
  disp_quadrangle_4.push_back(0.0338388);
  disp_quadrangle_4.push_back(0.0338388);
  Array<Real> velo_quadrangle_4(0,1); 
  velo_quadrangle_4.push_back(-307.221);
  velo_quadrangle_4.push_back(-307.221);

  // pointer to solution
  Array<Real> * disp = NULL;
  Array<Real> * velo = NULL;
  if (element_type == _triangle_3) {
    disp = &disp_triangle_3;
    velo = &velo_triangle_3;
  }
  else if (element_type == _quadrangle_4) {
    disp = &disp_quadrangle_4;
    velo = &velo_quadrangle_4;
  }

  for (UInt i=0; i<the_nodes.getSize(); ++i) {
    UInt node = the_nodes.storage()[i];
    if (!testFloat(model.getDisplacement().storage()[node*dim+1], disp->storage()[i], disp_tol)) {
      std::cout << "Node " << node << " has wrong disp. Computed = " << model.getDisplacement().storage()[node*dim+1] << " Solution = " << disp->storage()[i] << std::endl;
      return EXIT_FAILURE;
    }
    if (!testFloat(model.getVelocity().storage()[node*dim+1], velo->storage()[i], velo_tol)) {
      std::cout << "Node " << node << " has wrong velo. Computed = " << model.getVelocity().storage()[node*dim+1] << " Solution = " << velo->storage()[i] << std::endl;
      return EXIT_FAILURE;
    }
  }

  finalize();
  std::cout << "Patch test successful!" << std::endl;
  return EXIT_SUCCESS;
}

/* -------------------------------------------------------------------------- */
bool testFloat(Real a, Real b, Real adm_error) {                                      
  if (fabs(a-b) < adm_error)
    return true;
  else
    return false;
}

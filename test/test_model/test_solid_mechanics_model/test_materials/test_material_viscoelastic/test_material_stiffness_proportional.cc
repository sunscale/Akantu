/**
 * @file   test_material_stiffness_proportional.cc
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Sun Sep 12 23:37:43 2010
 *
 * @brief  test of the material elastic caughey
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
#include <limits>
#include <fstream>
#include <iostream>

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

int main(int argc, char *argv[])
{
  akantu::initialize(argc, argv);
  akantu::debug::setDebugLevel(akantu::dblWarning);

  const ElementType element_type = TYPE;
  UInt dim = Mesh::getSpatialDimension(element_type);

  /// load mesh
  Mesh mesh(dim);
  MeshIOMSH mesh_io;
  std::stringstream meshname_sstr; 
  meshname_sstr << "test_material_elastic_caughey_" << element_type << ".msh";
  mesh_io.read(meshname_sstr.str().c_str(), mesh);

  UInt max_steps = 2000;
  Real time_factor = 0.8;
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
  model.readMaterials("material_elastic_caughey.dat");
  model.initMaterials();

  std::cout << model.getMaterial(0) << std::endl;

  model.assembleMassLumped();

  /* ------------------------------------------------------------------------ */
  /* Boundary + initial conditions                                            */
  /* ------------------------------------------------------------------------ */
  Real eps = 1e-16;
  for (UInt i = 0; i < nb_nodes; ++i) {
    if(mesh.getNodes().storage()[dim*i] >= 9)
      model.getDisplacement().storage()[dim*i]
     	= (mesh.getNodes().storage()[dim*i] - 9) / 100.;

    if(mesh.getNodes().storage()[dim*i] <= eps)
      model.getBlockedDOFs().storage()[dim*i] = true;

    if(mesh.getNodes().storage()[dim*i + 1] <= eps ||
       mesh.getNodes().storage()[dim*i + 1] >= 1 - eps ) {
      model.getBlockedDOFs().storage()[dim*i + 1] = true;
    }
  }

  /// dump facet and surface information to paraview
#ifdef AKANTU_USE_IOHELPER
  iohelper::DumperParaview dumper;
  paraviewInit(dumper, model, element_type, "test_mat_el_caughey");
#endif //AKANTU_USE_IOHELPER

  std::stringstream filename_sstr;
  filename_sstr << "test_material_elastic_caughey_" << element_type << ".out";
  std::ofstream energy;
  energy.open(filename_sstr.str().c_str());
  energy << "id epot ekin tot" << std::endl;

  /// Setting time step
  Real time_step = model.getStableTimeStep() * time_factor;
  std::cout << "Time Step = " << time_step << "s" << std::endl;
  model.setTimeStep(time_step);

  /* ------------------------------------------------------------------------ */
  /* Main loop                                                                */
  /* ------------------------------------------------------------------------ */
  for(UInt s = 1; s <= max_steps; ++s) {

    model.explicitPred();
    model.updateResidual();
    model.updateAcceleration();
    model.explicitCorr();

    Real epot = model.getPotentialEnergy();
    Real ekin = model.getKineticEnergy();

    if(s % 10 == 0) {
       std::cerr << "passing step " << s << "/" << max_steps << std::endl;
       energy << s << " " << epot << " " << ekin << " " << epot + ekin
	      << std::endl;
    }

#ifdef AKANTU_USE_IOHELPER
    if(s % 10 == 0) paraviewDump(dumper);
#endif //AKANTU_USE_IOHELPER
  }

  finalize();
  return EXIT_SUCCESS;
}

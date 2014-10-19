/**
 * @file   test_solid_mechanics_model_bar_traction2d.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Sep 03 2010
 * @date last modification: Tue Sep 02 2014
 *
 * @brief  test of the class SolidMechanicsModel
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
#  include "io_helper.hh"


iohelper::ElemType paraview_type = iohelper::TRIANGLE2;
#endif //AKANTU_USE_IOHELPER

//#define CHECK_STRESS
akantu::ElementType type = akantu::_triangle_6;

akantu::SolidMechanicsModel * model;
akantu::UInt spatial_dimension = 2;
akantu::UInt nb_nodes;
akantu::UInt nb_element;
akantu::UInt nb_quadrature_points;

int main(int argc, char *argv[])
{
  akantu::initialize("material.dat", argc, argv);
  akantu::UInt max_steps = 5000;
  akantu::Real time_factor = 0.8;

  akantu::Mesh mesh(spatial_dimension);
  mesh.read("bar2.msh");

  akantu::SolidMechanicsModel model(mesh);

  nb_nodes = mesh.getNbNodes();
  nb_element = mesh.getNbElement(type);

  /// model initialization
  model.initFull();

  std::cout << model.getMaterial(0) << std::endl;

  model.initMaterials();
  model.assembleMassLumped();
  mesh.createGroupsFromMeshData<std::string>("physical_names");

  /// boundary conditions
  akantu::Real eps = 1e-16;
  const akantu::Array<akantu::Real> & pos = mesh.getNodes();
  akantu::Array<akantu::Real> & disp = model.getDisplacement();
  akantu::Array<bool> & boun = model.getBlockedDOFs();

  for (akantu::UInt i = 0; i < nb_nodes; ++i) {
    if(pos(i, 0) >= 9.) disp(i, 0) = (pos(i, 0) - 9) / 100.;
    if(pos(i) <= eps)   boun(i, 0) = true;
    if(pos(i, 1) <= eps || pos(i, 1) >= 1 - eps ) boun(i, 1) = true;
  }

  /// set the time step
  akantu::Real time_step = model.getStableTimeStep() * time_factor;
  std::cout << "Time Step = " << time_step << "s" << std::endl;
  model.setTimeStep(time_step);

  /// initialize the paraview output
  model.updateResidual();
  mesh.setBaseName("bar_traction_2d");
  model.addDumpField("displacement");
  model.addDumpField("mass"        );
  model.addDumpField("velocity"    );
  model.addDumpField("acceleration");
  model.addDumpField("force"       );
  model.addDumpField("residual"    );
  model.addDumpFieldTensor("stress");
  model.addDumpField("grad_u"       );

  model.addDumpGroupField("displacement","Top");
  model.dumpGroup("Top");
  model.dump();

#ifdef CHECK_STRESS
  std::ofstream outfile;
  outfile.open("stress");
#endif // CHECK_STRESS

  std::ofstream energy;
  energy.open("energy_bar_2d.csv");
  energy << "id,rtime,epot,ekin,tot" << std::endl;

  for(akantu::UInt s = 1; s <= max_steps; ++s) {
    model.explicitPred();

    model.updateResidual();
    model.updateAcceleration();
    model.explicitCorr();

    akantu::Real epot = model.getPotentialEnergy();
    akantu::Real ekin = model.getKineticEnergy();

    energy << s << "," << (s-1)*time_step << "," << epot << "," << ekin << "," << epot + ekin
	   << std::endl;

#ifdef CHECK_STRESS
    /// search the position of the maximum of stress to determine the wave speed
    akantu::Real max_stress = std::numeric_limits<akantu::Real>::min();
    akantu::Real * stress = model.getMaterial(0).getStress(type).storage();
    for (akantu::UInt i = 0; i < nb_element; ++i) {
      if(max_stress < stress[i*spatial_dimension*spatial_dimension]) {
	max_stress = stress[i*spatial_dimension*spatial_dimension];
      }
    }

    akantu::Real * coord    = model.getFEEngine().getMesh().getNodes().storage();
    akantu::Real * disp_val = model.getDisplacement().storage();
    akantu::UInt * conn     = model.getFEEngine().getMesh().getConnectivity(type).storage();
    akantu::UInt nb_nodes_per_element = model.getFEEngine().getMesh().getNbNodesPerElement(type);
    akantu::Real * coords = new akantu::Real[spatial_dimension];
    akantu::Real min_x = std::numeric_limits<akantu::Real>::max();
    akantu::Real max_x = std::numeric_limits<akantu::Real>::min();

    akantu::Real stress_range = 5e7;
    for (akantu::UInt el = 0; el < nb_element; ++el) {
      if(stress[el*spatial_dimension*spatial_dimension] > max_stress - stress_range) {
	akantu::UInt el_offset  = el * nb_nodes_per_element;
	memset(coords, 0, spatial_dimension*sizeof(akantu::Real));
	for (akantu::UInt n = 0; n < nb_nodes_per_element; ++n) {
	  for (akantu::UInt i = 0; i < spatial_dimension; ++i) {
	    akantu::UInt node = conn[el_offset + n] * spatial_dimension;
	    coords[i] += (coord[node + i] + disp_val[node + i])
	      / ((akantu::Real) nb_nodes_per_element);
	  }
	}
	min_x = min_x < coords[0] ? min_x : coords[0];
	max_x = max_x > coords[0] ? max_x : coords[0];
      }
    }

    outfile << s << " " << .5 * (min_x + max_x) << " " << min_x << " " << max_x << " " << max_x - min_x << " " << max_stress << std::endl;

    delete [] coords;
#endif // CHECK_STRESS


#ifdef AKANTU_USE_IOHELPER
    if(s % 100 == 0) {
      model.dump();
      model.dumpGroup();
    }
#endif //AKANTU_USE_IOHELPER
    if(s % 100 == 0) std::cout << "passing step " << s << "/" << max_steps << std::endl;
  }

  energy.close();

#ifdef CHECK_STRESS
  outfile.close();
#endif // CHECK_STRESS

  akantu::finalize();

  return EXIT_SUCCESS;
}

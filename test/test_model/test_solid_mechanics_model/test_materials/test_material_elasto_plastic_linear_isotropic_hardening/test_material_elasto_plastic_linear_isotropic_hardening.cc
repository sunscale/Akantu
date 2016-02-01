/**
 * @file   test_material_elasto_plastic_linear_isotropic_hardening.cc
 *
 * @author Jaehyun Cho <jaehyun.cho@epfl.ch>
 *
 * @date creation: Thu Dec 03 2015
 *
 * @brief  test for material type elasto plastic linear isotropic hardening using
 * #         tension-compression test
 *
 * @section LICENSE
 *
 * Copyright (©) 2015 EPFL (Ecole Polytechnique Fédérale de Lausanne) Laboratory
 * (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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
#include "solid_mechanics_model.hh"
#include <iostream>
using namespace akantu;

// /* -------------------------------------------------------------------------- */
const UInt spatial_dimension = 2;
const Real time_step = 1e-4;
const Real max_time  = 0.15;
/* -------------------------------------------------------------------------- */

int main(int argc, char *argv[]) {
  initialize("test_material_elasto_plastic_linear_isotropic_hardening.dat", argc, argv);
  Mesh mesh(spatial_dimension);
  StaticCommunicator & comm = StaticCommunicator::getStaticCommunicator();
  Int psize = comm.getNbProc();
  Int prank = comm.whoAmI();
  MeshPartition * partition = NULL;
  if(prank == 0) {
    mesh.read("test_material_elasto_plastic_linear_isotropic_hardening.msh");
    partition = new MeshPartitionScotch(mesh, spatial_dimension);
    partition->partitionate(psize);
  }

  SolidMechanicsModel model(mesh);
  model.initParallel(partition);
  mesh.createGroupsFromMeshData<std::string>("physical_names");

  /// model initialization
  model.initFull(SolidMechanicsModelOptions(_implicit_dynamic));
  Material &mat = model.getMaterial(0);
  Real E = mat.getParam<Real>("E");
  Real rho = mat.getParam<Real>("rho");
  Real h = mat.getParam<Real>("h");
  Real sigma_y = mat.getParam<Real>("sigma_y");
  
  mat.setParam("E",   E);
  mat.setParam("rho", rho);
  mat.setParam("h", h);
  mat.setParam("sigma_y", sigma_y);

  model.applyBC(BC::Dirichlet::FixedValue(0.0, _x), "left");
  model.applyBC(BC::Dirichlet::FixedValue(0.0, _y), "bottom");
  
  Mesh::type_iterator it = mesh.firstType(spatial_dimension);
  Mesh::type_iterator lastType = mesh.lastType(spatial_dimension);

  model.setBaseName("dynamic");
  model.addDumpFieldVector("displacement");
  model.addDumpField("velocity"    );
  model.addDumpField("acceleration");
  model.addDumpField("stress"       );
  model.addDumpField("strain"    );
  model.dump();


  std::ofstream output2;
  output2.open("strain-stress.txt");
  if(!output2.good()) AKANTU_DEBUG_ERROR("Cannot open file \"strain-stress.txt\"");
  output2 << "timestep, strain(0), stress(0)" << std::endl;
  output2 << "0 0.0 0.0 " << std::endl;
  model.setTimeStep(time_step);
  double dz = 1e-4;
  Real time = 0.;
  
  ElementTypeMapArray<Real> & byel_stress = model.flattenInternal("stress", _ek_regular);

  for (UInt s = 1; time < max_time; ++s, time += time_step) {
    if(prank == 0)
      std::cout << "Traction by " << dz*s << " \r" << std::flush;
    model.solveStep<_scm_newton_raphson_tangent_modified, _scc_increment>(1e-12, 100);
    model.applyBC(BC::Dirichlet::IncrementValue(dz, _x), "right");

    if(s % 10 == 0){
      model.dump();

      Real strainxx = 0.0;
      Real stressxx = 0.0;
      
      for(it = mesh.firstType(spatial_dimension); it != lastType; ++it){
  	Array<Real> & stress = byel_stress(*it);
 	for(UInt quad = 0; quad < stress.getSize(); ++quad){
	  stressxx += stress(quad,0);
  	}
	strainxx = dz*s/10.0;
	stressxx = stressxx/stress.getSize();
	output2 << s << " " << strainxx << " " << stressxx << std::endl;
      }
    }
  }
  output2.close();
  
  finalize();


  return EXIT_SUCCESS;
}

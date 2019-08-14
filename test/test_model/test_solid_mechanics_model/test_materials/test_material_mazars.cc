/**
 * @file   test_material_mazars.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Clement Roux <clement.roux@epfl.ch>
 *
 * @date creation: Thu Oct 08 2015
 * @date last modification: Sun Jul 09 2017
 *
 * @brief  test for material mazars, dissymmetric
 *
 * @section LICENSE
 *
 * Copyright (©) 2015-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#include "mesh_accessor.hh"
#include "solid_mechanics_model.hh"
/* -------------------------------------------------------------------------- */
#include <fstream>
/* -------------------------------------------------------------------------- */

using namespace akantu;

/* -------------------------------------------------------------------------- */
int main(int argc, char * argv[]) {
  debug::setDebugLevel(akantu::dblWarning);

  akantu::initialize("material_mazars.dat", argc, argv);
  const UInt spatial_dimension = 3;

  //  ElementType type = _quadrangle_4;
  ElementType type = _hexahedron_8;

  Mesh mesh(spatial_dimension);

  MeshAccessor mesh_accessor(mesh);

  Array<Real> & nodes = mesh_accessor.getNodes();
  Array<UInt> & connectivity = mesh_accessor.getConnectivity(type);

  const Real width = 1;
  UInt nb_dof = 0;

  connectivity.resize(1);

  if (type == _hexahedron_8) {
    nb_dof = 8;
    nodes.resize(nb_dof);

    nodes(0, 0) = 0.;
    nodes(0, 1) = 0.;
    nodes(0, 2) = 0.;

    nodes(1, 0) = width;
    nodes(1, 1) = 0.;
    nodes(1, 2) = 0.;

    nodes(2, 0) = width;
    nodes(2, 1) = width;
    nodes(2, 2) = 0;

    nodes(3, 0) = 0;
    nodes(3, 1) = width;
    nodes(3, 2) = 0;

    nodes(4, 0) = 0.;
    nodes(4, 1) = 0.;
    nodes(4, 2) = width;

    nodes(5, 0) = width;
    nodes(5, 1) = 0.;
    nodes(5, 2) = width;

    nodes(6, 0) = width;
    nodes(6, 1) = width;
    nodes(6, 2) = width;

    nodes(7, 0) = 0;
    nodes(7, 1) = width;
    nodes(7, 2) = width;

    connectivity(0, 0) = 0;
    connectivity(0, 1) = 1;
    connectivity(0, 2) = 2;
    connectivity(0, 3) = 3;
    connectivity(0, 4) = 4;
    connectivity(0, 5) = 5;
    connectivity(0, 6) = 6;
    connectivity(0, 7) = 7;
  } else if (type == _quadrangle_4) {
    nb_dof = 4;
    nodes.resize(nb_dof);

    nodes(0, 0) = 0.;
    nodes(0, 1) = 0.;

    nodes(1, 0) = width;
    nodes(1, 1) = 0;

    nodes(2, 0) = width;
    nodes(2, 1) = width;

    nodes(3, 0) = 0.;
    nodes(3, 1) = width;

    connectivity(0, 0) = 0;
    connectivity(0, 1) = 1;
    connectivity(0, 2) = 2;
    connectivity(0, 3) = 3;
  }

  mesh_accessor.makeReady();

  SolidMechanicsModel model(mesh);
  model.initFull();
  Material & mat = model.getMaterial(0);
  std::cout << mat << std::endl;

  /// boundary conditions
  Array<Real> & disp = model.getDisplacement();
  Array<Real> & velo = model.getVelocity();
  Array<bool> & boun = model.getBlockedDOFs();

  for (UInt i = 0; i < nb_dof; ++i) {
    boun(i, 0) = true;
  }

  Real time_step = model.getStableTimeStep() * .5;
  // Real time_step = 1e-5;

  std::cout << "Time Step = " << time_step
            << "s - nb elements : " << mesh.getNbElement(type) << std::endl;
  model.setTimeStep(time_step);

  std::ofstream energy;
  energy.open("energies_and_scalar_mazars.csv");
  energy << "id,rtime,epot,ekin,u,wext,kin+pot,D,strain_xx,strain_yy,stress_xx,"
            "stress_yy,edis,tot"
         << std::endl;

  /// Set dumper
  model.setBaseName("uniaxial_comp-trac_mazars");
  model.addDumpFieldVector("displacement");
  model.addDumpField("velocity");
  model.addDumpField("acceleration");
  model.addDumpField("damage");
  model.addDumpField("strain");
  model.addDumpField("stress");
  model.addDumpField("partitions");
  model.dump();

  const Array<Real> & strain = mat.getGradU(type);
  const Array<Real> & stress = mat.getStress(type);
  const Array<Real> & damage = mat.getArray<Real>("damage", type);

  /* ------------------------------------------------------------------------ */
  /* Main loop                                                                */
  /* ------------------------------------------------------------------------ */
  Real wext = 0.;
  Real sigma_max = 0, sigma_min = 0;

  Real max_disp;
  Real stress_xx_compression_1;
  UInt nb_steps = 7e5 / 150;

  Real adisp = 0;
  for (UInt s = 0; s < nb_steps; ++s) {
    if (s == 0) {
      max_disp = 0.003;
      adisp = -(max_disp * 8. / nb_steps) / 2.;
      std::cout << "Step " << s << " compression: " << max_disp << std::endl;
    }

    if (s == (2 * nb_steps / 8)) {
      stress_xx_compression_1 = stress(0, 0);
      max_disp = 0 + max_disp;
      adisp = max_disp * 8. / nb_steps;
      std::cout << "Step " << s << " discharge" << std::endl;
    }

    if (s == (3 * nb_steps / 8)) {
      max_disp = 0.004;
      adisp = -max_disp * 8. / nb_steps;
      std::cout << "Step " << s << " compression: " << max_disp << std::endl;
    }

    if (s == (4 * nb_steps / 8)) {
      if (stress(0, 0) < stress_xx_compression_1) {
        std::cout << "after this second compression step softening should have "
                     "started"
                  << std::endl;
        return EXIT_FAILURE;
      }
      max_disp = 0.0015 + max_disp;
      adisp = max_disp * 8. / nb_steps;
      std::cout << "Step " << s << " discharge tension: " << max_disp
                << std::endl;
    }

    if (s == (5 * nb_steps / 8)) {
      max_disp = 0. + 0.0005;
      adisp = -max_disp * 8. / nb_steps;
      std::cout << "Step " << s << " discharge" << std::endl;
    }

    if (s == (6 * nb_steps / 8)) {
      max_disp = 0.0035 - 0.001;
      adisp = max_disp * 8. / nb_steps;
      std::cout << "Step " << s << " tension: " << max_disp << std::endl;
    }

    if (s == (7 * nb_steps / 8)) {
      // max_disp = max_disp;
      adisp = -max_disp * 8. / nb_steps;
      std::cout << "Step " << s << " discharge" << std::endl;
    }

    for (UInt i = 0; i < nb_dof; ++i) {
      if (std::abs(nodes(i, 0) - width) <
          std::numeric_limits<Real>::epsilon()) {
        disp(i, 0) += adisp;
        velo(i, 0) = adisp / time_step;
      }
    }

    std::cout << "S: " << s << "/" << nb_steps << " inc disp: " << adisp
              << " disp: " << std::setw(5) << disp(2, 0) << "\r" << std::flush;

    model.solveStep();

    Real epot = model.getEnergy("potential");
    Real ekin = model.getEnergy("kinetic");
    Real edis = model.getEnergy("dissipated");
    wext += model.getEnergy("external work");

    sigma_max = std::max(sigma_max, stress(0, 0));
    sigma_min = std::min(sigma_min, stress(0, 0));
    if (s % 10 == 0)
      energy << s << ","             // 1
             << s * time_step << "," // 2
             << epot << ","          // 3
             << ekin << ","          // 4
             << disp(2, 0) << ","    // 5
             << wext << ","          // 6
             << epot + ekin << ","   // 7
             << damage(0) << ","     // 8
             << strain(0, 0) << ","  // 9
             << strain(0, 3) << ","  // 11
             << stress(0, 0) << ","  // 10
             << stress(0, 3) << ","  // 10
             << edis << ","          // 12
             << epot + ekin + edis   // 13
             << std::endl;

    if (s % 100 == 0)
      model.dump();
  }

  std::cout << std::endl
            << "sigma_max = " << sigma_max << ", sigma_min = " << sigma_min
            << std::endl;
  /// Verif the maximal/minimal stress values
  if ((std::abs(sigma_max) > std::abs(sigma_min)) or
      (std::abs(sigma_max - 6.24e6) > 1e5) or
      (std::abs(sigma_min + 2.943e7) > 1e6))
    return EXIT_FAILURE;
  energy.close();

  akantu::finalize();

  return EXIT_SUCCESS;
}

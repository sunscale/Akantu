/**
 * @file   patch_test_explicit.cc
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Cyprien Wolff <cyprien.wolff@epfl.ch>
 *
 * @date creation: Sat Apr 16 2011
 * @date last modification: Thu Oct 15 2015
 *
 * @brief  patch test for elastic material in solid mechanics model
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2012, 2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de
 * Lausanne)  Laboratory (LSMS  -  Laboratoire de  Simulation  en Mécanique  des
 * Solides)
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
#include "mesh_utils.hh"
/* -------------------------------------------------------------------------- */
#include <fstream>
/* -------------------------------------------------------------------------- */

using namespace akantu;

Real alpha[3][4] = {{0.01, 0.02, 0.03, 0.04},
                    {0.05, 0.06, 0.07, 0.08},
                    {0.09, 0.10, 0.11, 0.12}};

/* -------------------------------------------------------------------------- */
template <ElementType type, bool plane_strain>
static Matrix<Real> prescribed_strain() {
  UInt spatial_dimension = ElementClass<type>::getSpatialDimension();
  Matrix<Real> strain(spatial_dimension, spatial_dimension);

  for (UInt i = 0; i < spatial_dimension; ++i) {
    for (UInt j = 0; j < spatial_dimension; ++j) {
      strain(i, j) = alpha[i][j + 1];
    }
  }
  return strain;
}

template <ElementType type, bool is_plane_strain>
static Matrix<Real> prescribed_stress() {
  UInt spatial_dimension = ElementClass<type>::getSpatialDimension();
  Matrix<Real> stress(spatial_dimension, spatial_dimension);

  // plane strain in 2d
  Matrix<Real> strain(spatial_dimension, spatial_dimension);
  Matrix<Real> pstrain;
  pstrain = prescribed_strain<type, is_plane_strain>();
  Real nu = 0.3;
  Real E = 2.1e11;
  Real trace = 0;

  /// symetric part of the strain tensor
  for (UInt i = 0; i < spatial_dimension; ++i)
    for (UInt j = 0; j < spatial_dimension; ++j)
      strain(i, j) = 0.5 * (pstrain(i, j) + pstrain(j, i));

  for (UInt i = 0; i < spatial_dimension; ++i)
    trace += strain(i, i);

  if (spatial_dimension == 1) {
    stress(0, 0) = E * strain(0, 0);
  } else {
    if (is_plane_strain) {
      Real Ep = E / (1 + nu);
      for (UInt i = 0; i < spatial_dimension; ++i)
        for (UInt j = 0; j < spatial_dimension; ++j) {
          stress(i, j) = Ep * strain(i, j);
          if (i == j)
            stress(i, j) += Ep * (nu / (1 - 2 * nu)) * trace;
        }
    } else {
      Real Ep = E / (1 + nu);
      for (UInt i = 0; i < spatial_dimension; ++i)
        for (UInt j = 0; j < spatial_dimension; ++j) {
          stress(i, j) = Ep * strain(i, j);
          if (i == j)
            stress(i, j) += (nu * E) / (1 - (nu * nu)) * trace;
        }
    }
  }

  return stress;
}

/* -------------------------------------------------------------------------- */
int main(int argc, char * argv[]) {
  std::string input_file;
  if (PLANE_STRAIN)
    input_file = "material_check_stress_plane_strain.dat";
  else
    input_file = "material_check_stress_plane_stress.dat";

  initialize(input_file, argc, argv);

  debug::setDebugLevel(dblWarning);
  UInt dim = ElementClass<TYPE>::getSpatialDimension();
  const ElementType element_type = TYPE;

  UInt damping_steps = 600000;
  UInt damping_interval = 50;
  Real damping_ratio = 0.99;

  UInt additional_steps = 20000;
  UInt max_steps = damping_steps + additional_steps;

  /// load mesh
  Mesh my_mesh(dim);

  std::stringstream filename;
  filename << TYPE << ".msh";
  my_mesh.read(filename.str());

  UInt nb_nodes = my_mesh.getNbNodes();

  /// declaration of model
  SolidMechanicsModel my_model(my_mesh);
  /// model initialization
  my_model.initFull();

  std::cout << my_model.getMaterial(0) << std::endl;
  Real time_step = my_model.getStableTimeStep() / 100.;
  my_model.setTimeStep(time_step);
  my_model.assembleMassLumped();

  std::cout << "The number of time steps is: " << max_steps << " (" << time_step
            << "s)" << std::endl;

  // boundary conditions
  const Array<Real> & coordinates = my_mesh.getNodes();
  Array<Real> & displacement = my_model.getDisplacement();
  Array<bool> & boundary = my_model.getBlockedDOFs();

  MeshUtils::buildFacets(my_mesh);

  my_mesh.createBoundaryGroupFromGeometry();

  // Loop over (Sub)Boundar(ies) to block the nodes
  for (GroupManager::const_element_group_iterator it(
           my_mesh.element_group_begin());
       it != my_mesh.element_group_end(); ++it)
    for (const auto & node : it->second->getNodeGroup())
      for (UInt i = 0; i < dim; ++i)
        boundary(node, i) = true;

  // set the position of all nodes to the static solution
  for (UInt n = 0; n < nb_nodes; ++n) {
    for (UInt i = 0; i < dim; ++i) {
      displacement(n, i) = alpha[i][0];
      for (UInt j = 0; j < dim; ++j) {
        displacement(n, i) += alpha[i][j + 1] * coordinates(n, j);
      }
    }
  }

  Array<Real> & velocity = my_model.getVelocity();

  std::ofstream energy;
  std::stringstream energy_filename;
  energy_filename << "energy_" << TYPE << ".csv";
  energy.open(energy_filename.str().c_str());
  energy << "id,time,ekin" << std::endl;
  Real ekin_mean = 0.;

//#define DEBUG_TEST
#ifdef DEBUG_TEST
  my_model.solveStep();
  my_model.addDumpField("strain");
  my_model.addDumpField("stress");
  my_model.addDumpField("external_force");
  my_model.addDumpField("internal_force");
  my_model.addDumpField("velocity");
  my_model.addDumpField("acceleration");
  my_model.addDumpField("displacement");
  my_model.dump();
#endif
  /* ------------------------------------------------------------------------ */
  /* Main loop                                                                */
  /* ------------------------------------------------------------------------ */
  UInt s;
  for (s = 1; s <= max_steps; ++s) {
    if (s % 10000 == 0)
      std::cout << "passing step " << s << "/" << max_steps << " ("
                << s * time_step << "s)" << std::endl;

    // damp velocity in order to find equilibrium
    if ((s < damping_steps) && (s % damping_interval == 0)) {
      velocity *= damping_ratio;
    }

    if (s % 1000 == 0) {
      ekin_mean = ekin_mean / 1000.;
      std::cout << "Ekin mean = " << ekin_mean << std::endl;
      if (ekin_mean < 1e-10)
        break;
      ekin_mean = 0.;
    }

    my_model.solveStep();

    akantu::Real ekin = my_model.getEnergy("kinetic");
    ekin_mean += ekin;

    if (s % 1000 == 0)
      energy << s << "," << s * time_step << "," << ekin << std::endl;
  }

  energy.close();

  UInt nb_quadrature_points =
      my_model.getFEEngine().getNbIntegrationPoints(TYPE);
  Array<Real> & stress_vect = const_cast<Array<Real> &>(
      my_model.getMaterial(0).getStress(element_type));
  Array<Real> & strain_vect =
      const_cast<Array<Real> &>(my_model.getMaterial(0).getGradU(element_type));

  Array<Real>::matrix_iterator stress_it = stress_vect.begin(dim, dim);
  Array<Real>::matrix_iterator strain_it = strain_vect.begin(dim, dim);

  Matrix<Real> presc_stress;
  presc_stress = prescribed_stress<TYPE, PLANE_STRAIN>();
  Matrix<Real> presc_strain;
  presc_strain = prescribed_strain<TYPE, PLANE_STRAIN>();

  UInt nb_element = my_mesh.getNbElement(TYPE);

  Real strain_tolerance = 1e-13;
  Real stress_tolerance = 1e-13;

  for (UInt el = 0; el < nb_element; ++el) {
    for (UInt q = 0; q < nb_quadrature_points; ++q) {
      Matrix<Real> & stress = *stress_it;
      Matrix<Real> & strain = *strain_it;

      Matrix<Real> diff(dim, dim);

      diff = strain;
      diff -= presc_strain;
      Real strain_error = diff.norm<L_inf>() / strain.norm<L_inf>();

      if (strain_error > strain_tolerance) {
        std::cerr << "strain error: " << strain_error << " > "
                  << strain_tolerance << std::endl;
        std::cerr << "strain: " << strain << std::endl
                  << "prescribed strain: " << presc_strain << std::endl;
        return EXIT_FAILURE;
      } else {
        std::cerr << "strain error: " << strain_error << " < "
                  << strain_tolerance << std::endl;
      }

      diff = stress;
      diff -= presc_stress;
      Real stress_error = diff.norm<L_inf>() / stress.norm<L_inf>();

      if (stress_error > stress_tolerance) {
        std::cerr << "stress error: " << stress_error << " > "
                  << stress_tolerance << std::endl;
        std::cerr << "stress: " << stress << std::endl
                  << "prescribed stress: " << presc_stress << std::endl;
        return EXIT_FAILURE;
      } else {
        std::cerr << "stress error: " << stress_error << " < "
                  << stress_tolerance << std::endl;
      }

      ++stress_it;
      ++strain_it;
    }
  }

  for (UInt n = 0; n < nb_nodes; ++n) {
    for (UInt i = 0; i < dim; ++i) {
      Real disp = alpha[i][0];
      for (UInt j = 0; j < dim; ++j) {
        disp += alpha[i][j + 1] * coordinates(n, j);
      }

      if (!(std::abs(displacement(n, i) - disp) < 1e-7)) {
        std::cerr << "displacement(" << n << ", " << i
                  << ")=" << displacement(n, i) << " should be equal to "
                  << disp << "(" << displacement(n, i) - disp << ")"
                  << std::endl;
        return EXIT_FAILURE;
      }
    }
  }

  finalize();

  return EXIT_SUCCESS;
}

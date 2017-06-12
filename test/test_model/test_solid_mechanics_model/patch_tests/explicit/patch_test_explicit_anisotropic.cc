/**
 * @file   patch_test_explicit_anisotropic.cc
 *
 * @author Till Junge <till.junge@epfl.ch>
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
#include <iostream>

#include "solid_mechanics_model.hh"

using namespace akantu;

Real alpha[3][4] = {{0.01, 0.02, 0.03, 0.04},
                    {0.05, 0.06, 0.07, 0.08},
                    {0.09, 0.10, 0.11, 0.12}};

// Stiffness tensor, rotated by hand
Real C[3][3][3][3] = {
    {{{112.93753505, 1.85842452538e-10, -4.47654358027e-10},
      {1.85847317471e-10, 54.2334345331, -3.69840984824},
      {-4.4764768395e-10, -3.69840984824, 56.848605217}},
     {{1.85847781609e-10, 25.429294233, -3.69840984816},
      {25.429294233, 3.31613847493e-10, -8.38797920011e-11},
      {-3.69840984816, -8.38804581349e-11, -1.97875715813e-10}},
     {{-4.47654358027e-10, -3.69840984816, 28.044464917},
      {-3.69840984816, 2.09374961813e-10, 9.4857455224e-12},
      {28.044464917, 9.48308098714e-12, -2.1367885239e-10}}},
    {{{1.85847781609e-10, 25.429294233, -3.69840984816},
      {25.429294233, 3.31613847493e-10, -8.38793479119e-11},
      {-3.69840984816, -8.38795699565e-11, -1.97876381947e-10}},
     {{54.2334345331, 3.31617400207e-10, 2.09372075233e-10},
      {3.3161562385e-10, 115.552705733, -3.15093728886e-10},
      {2.09372075233e-10, -3.15090176173e-10, 54.2334345333}},
     {{-3.69840984824, -8.38795699565e-11, 9.48219280872e-12},
      {-8.38795699565e-11, -3.1509195253e-10, 25.4292942335},
      {9.48441325477e-12, 25.4292942335, 3.69840984851}}},
    {{{-4.47653469848e-10, -3.69840984816, 28.044464917},
      {-3.69840984816, 2.09374073634e-10, 9.48752187924e-12},
      {28.044464917, 9.48552347779e-12, -2.1367885239e-10}},
     {{-3.69840984824, -8.3884899027e-11, 9.48219280872e-12},
      {-8.3884899027e-11, -3.150972816e-10, 25.4292942335},
      {9.48041645188e-12, 25.4292942335, 3.69840984851}},
     {{56.848605217, -1.97875493768e-10, -2.13681516925e-10},
      {-1.97877270125e-10, 54.2334345333, 3.69840984851},
      {-2.13683293282e-10, 3.69840984851, 112.93753505}}}};

/* -------------------------------------------------------------------------- */
template <ElementType type> static Matrix<Real> prescribed_grad_u() {
  UInt spatial_dimension = ElementClass<type>::getSpatialDimension();
  Matrix<Real> grad_u(spatial_dimension, spatial_dimension);

  for (UInt i = 0; i < spatial_dimension; ++i) {
    for (UInt j = 0; j < spatial_dimension; ++j) {
      grad_u(i, j) = alpha[i][j + 1];
    }
  }
  return grad_u;
}

template <ElementType type> static Matrix<Real> prescribed_stress() {
  UInt spatial_dimension = ElementClass<type>::getSpatialDimension();
  Matrix<Real> stress(spatial_dimension, spatial_dimension);

  // plane strain in 2d
  Matrix<Real> strain(spatial_dimension, spatial_dimension);
  Matrix<Real> pstrain;
  pstrain = prescribed_grad_u<type>();
  Real trace = 0;

  /// symetric part of the strain tensor
  for (UInt i = 0; i < spatial_dimension; ++i)
    for (UInt j = 0; j < spatial_dimension; ++j)
      strain(i, j) = 0.5 * (pstrain(i, j) + pstrain(j, i));

  for (UInt i = 0; i < spatial_dimension; ++i)
    trace += strain(i, i);

  for (UInt i = 0; i < spatial_dimension; ++i) {
    for (UInt j = 0; j < spatial_dimension; ++j) {
      stress(i, j) = 0;
      for (UInt k = 0; k < spatial_dimension; ++k) {
        for (UInt l = 0; l < spatial_dimension; ++l) {
          stress(i, j) += C[i][j][k][l] * strain(k, l);
        }
      }
    }
  }
  return stress;
}

#define TYPE _tetrahedron_4
/* -------------------------------------------------------------------------- */
int main(int argc, char * argv[]) {
  initialize("material_anisotropic.dat", argc, argv);

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
  SolidMechanicsModel model(my_mesh);
  /// model initialization
  model.initFull(SolidMechanicsModelOptions(_explicit_lumped_mass));

  std::cout << model.getMaterial(0) << std::endl;
  Real time_step = model.getStableTimeStep() / 5.;
  model.setTimeStep(time_step);

  std::cout << "The number of time steps is: " << max_steps << " (" << time_step
            << "s)" << std::endl;

  // boundary conditions
  const Array<Real> & coordinates = my_mesh.getNodes();
  Array<Real> & displacement = model.getDisplacement();
  Array<bool> & boundary = model.getBlockedDOFs();

  MeshUtils::buildFacets(my_mesh);

  my_mesh.createBoundaryGroupFromGeometry();

  // Loop over (Sub)Boundar(ies)
  // Loop over (Sub)Boundar(ies)
  for (GroupManager::const_element_group_iterator it(
           my_mesh.element_group_begin());
       it != my_mesh.element_group_end(); ++it) {
    for (ElementGroup::const_node_iterator nodes_it(it->second->node_begin());
         nodes_it != it->second->node_end(); ++nodes_it) {
      UInt n(*nodes_it);
      std::cout << "Node " << *nodes_it << std::endl;
      for (UInt i = 0; i < dim; ++i) {
        displacement(n, i) = alpha[i][0];
        for (UInt j = 0; j < dim; ++j) {
          displacement(n, i) += alpha[i][j + 1] * coordinates(n, j);
        }
        boundary(n, i) = true;
      }
    }
  }

  // Actually, loop over all nodes, since I wanna test a static solution
  for (UInt n = 0; n < nb_nodes; ++n) {
    for (UInt i = 0; i < dim; ++i) {
      displacement(n, i) = alpha[i][0];
      for (UInt j = 0; j < dim; ++j) {
        displacement(n, i) += alpha[i][j + 1] * coordinates(n, j);
      }
    }
  }

  Array<Real> & velocity = model.getVelocity();

  std::ofstream energy;
  std::stringstream energy_filename;
  energy_filename << "energy_" << TYPE << ".csv";
  energy.open(energy_filename.str().c_str());
  energy << "id,time,ekin" << std::endl;
  Real ekin_mean = 0.;

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

    model.solveStep();

    Real ekin = model.getEnergy("kinetic");
    ekin_mean += ekin;

    if (s % 1000 == 0)
      energy << s << "," << s * time_step << "," << ekin << std::endl;
  }

  energy.close();

  UInt nb_quadrature_points = model.getFEEngine().getNbIntegrationPoints(TYPE);
  Array<Real> & stress_vect =
      const_cast<Array<Real> &>(model.getMaterial(0).getStress(element_type));
  Array<Real> & gradu_vect =
      const_cast<Array<Real> &>(model.getMaterial(0).getGradU(element_type));

  Array<Real>::iterator<Matrix<Real>> stress_it = stress_vect.begin(dim, dim);
  Array<Real>::iterator<Matrix<Real>> gradu_it = gradu_vect.begin(dim, dim);

  Matrix<Real> presc_stress;
  presc_stress = prescribed_stress<TYPE>();
  Matrix<Real> presc_gradu;
  presc_gradu = prescribed_grad_u<TYPE>();

  UInt nb_element = my_mesh.getNbElement(TYPE);

  Real gradu_tolerance = 1e-9;
  Real stress_tolerance = 1e-2;
  if (s > max_steps) {
    stress_tolerance = 1e-4;
    gradu_tolerance = 1e-7;
  }

  for (UInt el = 0; el < nb_element; ++el) {
    for (UInt q = 0; q < nb_quadrature_points; ++q) {
      Matrix<Real> & stress = *stress_it;
      Matrix<Real> & gradu = *gradu_it;

      for (UInt i = 0; i < dim; ++i) {
        for (UInt j = 0; j < dim; ++j) {
          if (!(std::abs(gradu(i, j) - presc_gradu(i, j)) < gradu_tolerance)) {
            std::cerr << "gradu[" << i << "," << j << "] = " << gradu(i, j)
                      << " but should be = " << presc_gradu(i, j) << " (-"
                      << std::abs(gradu(i, j) - presc_gradu(i, j))
                      << ") [el : " << el << " - q : " << q << "]" << std::endl;
            std::cerr << gradu << presc_gradu << std::endl;
            return EXIT_FAILURE;
          }

          if (!(std::abs(stress(i, j) - presc_stress(i, j)) <
                stress_tolerance)) {
            std::cerr << "stress[" << i << "," << j << "] = " << stress(i, j)
                      << " but should be = " << presc_stress(i, j) << " (-"
                      << std::abs(stress(i, j) - presc_stress(i, j))
                      << ") [el : " << el << " - q : " << q << "]" << std::endl;
            std::cerr << stress << presc_stress << std::endl;
            return EXIT_FAILURE;
          }
        }
      }

      ++stress_it;
      ++gradu_it;
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

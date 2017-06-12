/**
 * @file   patch_test_implicit.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
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
#include "non_linear_solver.hh"
#include "solid_mechanics_model.hh"
#include "sparse_matrix.hh"
/* -------------------------------------------------------------------------- */

using namespace akantu;

Real alpha[3][4] = {{0.01, 0.02, 0.03, 0.04},
                    {0.05, 0.06, 0.07, 0.08},
                    {0.09, 0.10, 0.11, 0.12}};

/* -------------------------------------------------------------------------- */
template <ElementType type, bool is_plane_strain>
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

  Real lambda = nu * E / ((1 + nu) * (1 - 2 * nu));
  Real mu = E / (2 * (1 + nu));

  if (!is_plane_strain) {
    std::cout << "toto" << std::endl;
    lambda = nu * E / (1 - nu * nu);
  }

  if (spatial_dimension == 1) {
    stress(0, 0) = E * strain(0, 0);
  } else {
    for (UInt i = 0; i < spatial_dimension; ++i)
      for (UInt j = 0; j < spatial_dimension; ++j) {
        stress(i, j) = (i == j) * lambda * trace + 2 * mu * strain(i, j);
      }
  }

  return stress;
}

/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
int main(int argc, char * argv[]) {
  std::string input_file;
  if (PLANE_STRAIN)
    input_file = "material_check_stress_plane_strain.dat";
  else
    input_file = "material_check_stress_plane_stress.dat";

  initialize(input_file, argc, argv);

  UInt dim = ElementClass<TYPE>::getSpatialDimension();
  const ElementType element_type = TYPE;

  /// load mesh
  Mesh mesh(dim);

  std::stringstream filename;
  filename << TYPE << ".msh";
  mesh.read(filename.str());

  UInt nb_nodes = mesh.getNbNodes();

  /// declaration of model
  SolidMechanicsModel model(mesh);
  /// model initialization
  model.initFull(SolidMechanicsModelOptions(_static));

  const Array<Real> & coordinates = mesh.getNodes();
  Array<Real> & displacement = model.getDisplacement();
  Array<bool> & boundary = model.getBlockedDOFs();
  MeshUtils::buildFacets(mesh);

  mesh.createBoundaryGroupFromGeometry();

  // Loop over (Sub)Boundar(ies)
  for (GroupManager::const_element_group_iterator it(
           mesh.element_group_begin());
       it != mesh.element_group_end(); ++it) {
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

  /* ------------------------------------------------------------------------ */
  /* Static solve                                                             */
  /* ------------------------------------------------------------------------ */
  auto & solver = model.getNonLinearSolver();
  solver.set("max_iterations", 2);
  solver.set("threshold", 2e-4);
  solver.set("convergence_type", _scc_residual);

  model.solveStep();
  model.getDOFManager().getMatrix("K").saveMatrix("clown_matrix.mtx");

  /* ------------------------------------------------------------------------ */
  /* Checks                                                                   */
  /* ------------------------------------------------------------------------ */
  UInt nb_quadrature_points =
      model.getFEEngine().getNbIntegrationPoints(element_type);

  Array<Real> & stress_vect =
      const_cast<Array<Real> &>(model.getMaterial(0).getStress(element_type));
  Array<Real> & strain_vect =
      const_cast<Array<Real> &>(model.getMaterial(0).getGradU(element_type));

  Array<Real>::matrix_iterator stress_it = stress_vect.begin(dim, dim);
  Array<Real>::matrix_iterator strain_it = strain_vect.begin(dim, dim);

  Matrix<Real> presc_stress;
  presc_stress = prescribed_stress<TYPE, PLANE_STRAIN>();
  Matrix<Real> presc_strain;
  presc_strain = prescribed_strain<TYPE, PLANE_STRAIN>();

  UInt nb_element = mesh.getNbElement(TYPE);

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

      if (!(std::abs(displacement(n, i) - disp) < 2e-15)) {
        std::cerr << "displacement(" << n << ", " << i
                  << ")=" << displacement(n, i) << " should be equal to "
                  << disp << std::endl;
        return EXIT_FAILURE;
      }
    }
  }

  finalize();

  return EXIT_SUCCESS;
}

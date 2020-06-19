/**
 * @file   test_solid_mechanics_model_material_eigenstrain.cc
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Sat Apr 16 2011
 * @date last modification: Thu Feb 01 2018
 *
 * @brief  test the internal field prestrain
 *
 *
 * Copyright (©)  2010-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "mesh_utils.hh"
#include "non_linear_solver.hh"
#include "solid_mechanics_model.hh"
/* -------------------------------------------------------------------------- */

using namespace akantu;

Real alpha[3][4] = {{0.01, 0.02, 0.03, 0.04},
                    {0.05, 0.06, 0.07, 0.08},
                    {0.09, 0.10, 0.11, 0.12}};

/* -------------------------------------------------------------------------- */
template <ElementType type> static Matrix<Real> prescribed_strain() {
  UInt spatial_dimension = ElementClass<type>::getSpatialDimension();
  Matrix<Real> strain(spatial_dimension, spatial_dimension);

  for (UInt i = 0; i < spatial_dimension; ++i) {
    for (UInt j = 0; j < spatial_dimension; ++j) {
      strain(i, j) = alpha[i][j + 1];
    }
  }
  return strain;
}

template <ElementType type>
static Matrix<Real> prescribed_stress(Matrix<Real> prescribed_eigengradu) {
  UInt spatial_dimension = ElementClass<type>::getSpatialDimension();
  Matrix<Real> stress(spatial_dimension, spatial_dimension);

  // plane strain in 2d
  Matrix<Real> strain(spatial_dimension, spatial_dimension);
  Matrix<Real> pstrain;
  pstrain = prescribed_strain<type>();
  Real nu = 0.3;
  Real E = 2.1e11;
  Real trace = 0;

  /// symetric part of the strain tensor
  for (UInt i = 0; i < spatial_dimension; ++i)
    for (UInt j = 0; j < spatial_dimension; ++j)
      strain(i, j) = 0.5 * (pstrain(i, j) + pstrain(j, i));

  // elastic strain is equal to elastic strain minus the eigenstrain
  strain -= prescribed_eigengradu;
  for (UInt i = 0; i < spatial_dimension; ++i)
    trace += strain(i, i);

  Real lambda = nu * E / ((1 + nu) * (1 - 2 * nu));
  Real mu = E / (2 * (1 + nu));

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
  initialize("material_elastic_plane_strain.dat", argc, argv);

  UInt dim = 3;
  const ElementType element_type = _tetrahedron_4;

  Matrix<Real> prescribed_eigengradu(dim, dim);
  prescribed_eigengradu.set(0.1);

  /// load mesh
  Mesh mesh(dim);
  mesh.read("cube_3d_tet_4.msh");

  /// declaration of model
  SolidMechanicsModel model(mesh);
  /// model initialization
  model.initFull(_analysis_method = _static);

  // model.getNewSolver("static", TimeStepSolverType::_static,
  // NonLinearSolverType::_newton_raphson_modified);
  auto & solver = model.getNonLinearSolver("static");
  solver.set("threshold", 2e-4);
  solver.set("max_iterations", 2);
  solver.set("convergence_type", SolveConvergenceCriteria::_residual);

  const Array<Real> & coordinates = mesh.getNodes();
  Array<Real> & displacement = model.getDisplacement();
  Array<bool> & boundary = model.getBlockedDOFs();
  MeshUtils::buildFacets(mesh);

  mesh.createBoundaryGroupFromGeometry();

  // Loop over (Sub)Boundar(ies)
  for (auto & group : mesh.iterateElementGroups()) {
    for (const auto & n : group.getNodeGroup()) {
      std::cout << "Node " << n << std::endl;
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
  /* Apply eigenstrain in each element */
  /* ------------------------------------------------------------------------ */
  Array<Real> & eigengradu_vect =
      model.getMaterial(0).getInternal<Real>("eigen_grad_u")(element_type);
  auto eigengradu_it = eigengradu_vect.begin(dim, dim);
  auto eigengradu_end = eigengradu_vect.end(dim, dim);

  for (; eigengradu_it != eigengradu_end; ++eigengradu_it) {
    *eigengradu_it = prescribed_eigengradu;
  }

  /* ------------------------------------------------------------------------ */
  /* Static solve                                                             */
  /* ------------------------------------------------------------------------ */
  model.solveStep();

  std::cout << "Converged in " << Int(solver.get("nb_iterations")) << " ("
            << Real(solver.get("error")) << ")" << std::endl;

  /* ------------------------------------------------------------------------ */
  /* Checks                                                                   */
  /* ------------------------------------------------------------------------ */
  const Array<Real> & stress_vect =
      model.getMaterial(0).getStress(element_type);

  auto stress_it = stress_vect.begin(dim, dim);
  auto stress_end = stress_vect.end(dim, dim);

  Matrix<Real> presc_stress;
  presc_stress = prescribed_stress<element_type>(prescribed_eigengradu);

  Real stress_tolerance = 1e-13;

  for (; stress_it != stress_end; ++stress_it) {
    const auto & stress = *stress_it;
    Matrix<Real> diff(dim, dim);

    diff = stress;
    diff -= presc_stress;
    Real stress_error = diff.norm<L_inf>() / stress.norm<L_inf>();

    if (stress_error > stress_tolerance) {
      std::cerr << "stress error: " << stress_error << " > " << stress_tolerance
                << std::endl;
      std::cerr << "stress: " << stress << std::endl
                << "prescribed stress: " << presc_stress << std::endl;
      return EXIT_FAILURE;
    } // else {
    //   std::cerr << "stress error: " << stress_error
    //             << " < " << stress_tolerance << std::endl;
    // }
  }

  finalize();

  return EXIT_SUCCESS;
}

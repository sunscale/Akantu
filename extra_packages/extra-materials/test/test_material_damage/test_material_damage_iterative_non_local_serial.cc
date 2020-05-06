/**
 * @file   test_material_damage_iterative_non_local_serial.cc
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @date   Thu Nov 26 12:20:15 2015
 *
 * @brief  test the material damage iterative non local in serial
 *
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
#include "material_damage_iterative_non_local.hh"
#include "solid_mechanics_model.hh"
/* -------------------------------------------------------------------------- */
using namespace akantu;

/* -------------------------------------------------------------------------- */
/* Main                                                                       */
/* -------------------------------------------------------------------------- */
int main(int argc, char * argv[]) {
  Math::setTolerance(1e-13);
  debug::setDebugLevel(dblWarning);

  initialize("material_non_local.dat", argc, argv);

  const UInt spatial_dimension = 2;
  ElementType element_type = _triangle_3;
  /// read the mesh and partion it
  Mesh mesh(spatial_dimension);
  mesh.read("plate.msh");

  /// model creation
  SolidMechanicsModel model(mesh);

  /// initialization of the model
  model.initFull(SolidMechanicsModelOptions(_static));

  /// boundary conditions
  /// Dirichlet BC
  mesh.createGroupsFromMeshData<std::string>(
      "physical_names"); // creates groups from mesh names
  model.applyBC(BC::Dirichlet::FixedValue(0, _x), "left");
  model.applyBC(BC::Dirichlet::FixedValue(0, _y), "bottom");
  model.applyBC(BC::Dirichlet::FixedValue(2., _y), "top");

  /// add fields that should be dumped
  model.setBaseName("material_damage_iterative_test");
  model.addDumpFieldVector("displacement");
  ;
  model.addDumpField("stress");
  model.addDumpField("blocked_dofs");
  model.addDumpField("residual");
  model.addDumpField("grad_u");
  model.addDumpField("grad_u non local");
  model.addDumpField("damage");
  model.addDumpField("partitions");
  model.addDumpField("material_index");
  model.addDumpField("Sc");
  model.addDumpField("force");
  model.addDumpField("equivalent_stress");

  model.dump();

  MaterialDamageIterativeNonLocal<spatial_dimension> & material =
      dynamic_cast<MaterialDamageIterativeNonLocal<spatial_dimension> &>(
          model.getMaterial(0));

  Real error;
  bool converged = false;
  Real max_eq_stress = 0;

  /// solve the system
  converged =
      model.solveStep<_scm_newton_raphson_tangent_modified,
                      SolveConvergenceCriteria::_increment>(1e-4, error, 2);

  if (converged == false) {
    std::cout << "The error is: " << error << std::endl;
    AKANTU_DEBUG_ASSERT(converged, "Did not converge");
  }

  model.dump();

  /// check the non-local grad_u: since grad_u is constant everywhere
  /// also the grad_u non-local has to be constant
  Array<Real> & grad_u_nl =
      material.getInternal<Real>("grad_u non local")(element_type, _not_ghost);
  Array<Real>::const_matrix_iterator grad_u_nl_it =
      grad_u_nl.begin(spatial_dimension, spatial_dimension);
  Array<Real>::const_matrix_iterator grad_u_nl_end =
      grad_u_nl.end(spatial_dimension, spatial_dimension);
  Real diff = 0.;
  Matrix<Real> diff_matrix(spatial_dimension, spatial_dimension);
  Matrix<Real> const_grad_u(spatial_dimension, spatial_dimension, 0.);
  const_grad_u(1, 1) = 1.;

  for (; grad_u_nl_it != grad_u_nl_end; ++grad_u_nl_it) {
    diff_matrix = (*grad_u_nl_it) - const_grad_u;
    diff += diff_matrix.norm<L_2>();
  }

  if (diff > 10.e-13) {
    std::cout << "Error in the non-local grad_u computation" << std::endl;
    return EXIT_FAILURE;
  }

  /// change the displacement in one node to modify grad_u
  Array<Real> & displ = model.getDisplacement();
  displ(0, 1) = 2.6;

  /// compute stresses: this will average grad_u and compute the max. eq. stress
  model.updateResidual();
  model.dump();
  /// due to the change in the displacement element 33 and 37 will
  /// have a grad_u different then one
  const Array<Real> & grad_u =
      material.getInternal<Real>("grad_u")(element_type, _not_ghost);
  Array<Real>::const_matrix_iterator grad_u_it =
      grad_u.begin(spatial_dimension, spatial_dimension);
  Array<Real>::const_matrix_iterator grad_u_end =
      grad_u.end(spatial_dimension, spatial_dimension);
  diff = 0.;
  diff_matrix.clear();

  UInt counter = 0;
  for (; grad_u_it != grad_u_end; ++grad_u_it) {
    diff_matrix = (*grad_u_it) - const_grad_u;
    if (counter == 34 || counter == 38) {
      if ((diff_matrix.norm<L_2>()) < 0.1) {
        std::cout << "Error in the grad_u computation" << std::endl;
        return EXIT_FAILURE;
      }
    } else
      diff += diff_matrix.norm<L_2>();
    ++counter;
  }

  if (diff > 10.e-13) {
    std::cout << "Error in the grad_u computation" << std::endl;
    return EXIT_FAILURE;
  }

  /// check that the non-local grad_u
  diff = 0.;
  diff_matrix.clear();
  Real nl_radius = 1.0; /// same values as in material file
  grad_u_nl_it = grad_u_nl.begin(spatial_dimension, spatial_dimension);
  ElementTypeMapReal quad_coords("quad_coords");
  mesh.initElementTypeMapArray(quad_coords, spatial_dimension,
                               spatial_dimension, false, _ek_regular, true);
  model.getFEEngine().computeIntegrationPointsCoordinates(quad_coords);
  UInt nb_elements = mesh.getNbElement(element_type, _not_ghost);
  UInt nb_quads = model.getFEEngine().getNbIntegrationPoints(element_type);
  Array<Real> & coords = quad_coords(element_type, _not_ghost);
  auto coord_it = coords.begin(spatial_dimension);
  Vector<Real> q1(spatial_dimension);
  Vector<Real> q2(spatial_dimension);
  q1 = coord_it[34];
  q2 = coord_it[38];
  for (UInt e = 0; e < nb_elements; ++e) {
    for (UInt q = 0; q < nb_quads; ++q, ++coord_it, ++grad_u_nl_it) {
      diff_matrix = (*grad_u_nl_it) - const_grad_u;
      if ((q1.distance(*coord_it) <= (nl_radius + Math::getTolerance())) ||
          (q2.distance(*coord_it) <= (nl_radius + Math::getTolerance()))) {
        if ((diff_matrix.norm<L_2>()) < 1.e-6) {
          std::cout << (diff_matrix.norm<L_2>()) << std::endl;
          std::cout << "Error in the non-local grad_u computation" << std::endl;
          return EXIT_FAILURE;
        }
      } else
        diff += diff_matrix.norm<L_2>();
    }
  }

  if (diff > 10.e-13) {
    std::cout << "Error in the non-local grad_u computation" << std::endl;
    return EXIT_FAILURE;
  }

  /// make sure that the normalized equivalent stress is based on the
  /// non-local grad_u for this test check the elements that have the
  /// constant stress of 1 but different non-local gradu because they
  /// are in the neighborhood of the modified elements
  coord_it = coords.begin(spatial_dimension);
  const Array<Real> & eq_stress =
      material.getInternal<Real>("equivalent_stress")(element_type, _not_ghost);
  Array<Real>::const_scalar_iterator eq_stress_it = eq_stress.begin();
  counter = 0;
  for (UInt e = 0; e < nb_elements; ++e) {
    for (UInt q = 0; q < nb_quads;
         ++q, ++coord_it, ++grad_u_nl_it, ++eq_stress_it) {
      if (counter == 34 || counter == 38)
        continue;
      if (((q1.distance(*coord_it) <= (nl_radius + Math::getTolerance())) ||
           (q2.distance(*coord_it) <= (nl_radius + Math::getTolerance()))) &&
          Math::are_float_equal(*eq_stress_it, 0.1)) {
        std::cout << "the normalized equivalent stress is most likely based on "
                     "the local, not the non-local grad_u!!!!"
                  << std::endl;
        finalize();
        return EXIT_FAILURE;
      }
      ++counter;
    }
  }

  max_eq_stress = material.getNormMaxEquivalentStress();

  if (!Math::are_float_equal(max_eq_stress, 0.1311267235941873)) {
    std::cout << "the maximum equivalent stress is wrong" << std::endl;
    finalize();
    return EXIT_FAILURE;
  }

  model.dump();
  finalize();

  return EXIT_SUCCESS;
}

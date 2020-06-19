/**
 * @file   test_material_damage_iterative_non_local_parallel.cc
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @date   Thu Nov 26 12:20:15 2015
 *
 * @brief  test the material damage iterative non local in parallel
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
#include "material_damage_iterative.hh"
#include "solid_mechanics_model.hh"
/* -------------------------------------------------------------------------- */
using namespace akantu;

bool checkDisplacement(SolidMechanicsModel & model, ElementType type,
                       std::ofstream & error_output, UInt step,
                       bool barycenters);

/* -------------------------------------------------------------------------- */
/* Main                                                                       */
/* -------------------------------------------------------------------------- */
int main(int argc, char * argv[]) {

  debug::setDebugLevel(dblWarning);
  ElementType element_type = _triangle_3;

  initialize("two_materials.dat", argc, argv);

  const UInt spatial_dimension = 2;
  StaticCommunicator & comm =
      akantu::StaticCommunicator::getStaticCommunicator();
  Int psize = comm.getNbProc();
  Int prank = comm.whoAmI();

  /// read the mesh and partion it
  Mesh mesh(spatial_dimension);
  akantu::MeshPartition * partition = NULL;

  if (prank == 0) {

    mesh.read("one_circular_inclusion.msh");

    /// partition the mesh
    partition = new MeshPartitionScotch(mesh, spatial_dimension);

    partition->partitionate(psize);
  }

  /// model creation
  SolidMechanicsModel model(mesh);
  model.initParallel(partition);
  delete partition;

  /// assign the material
  MeshDataMaterialSelector<std::string> * mat_selector;
  mat_selector =
      new MeshDataMaterialSelector<std::string>("physical_names", model);
  model.setMaterialSelector(*mat_selector);
  mesh.createGroupsFromMeshData<std::string>(
      "physical_names"); // creates groups from mesh names
  /// initialization of the model
  model.initFull(SolidMechanicsModelOptions(_static));

  /// boundary conditions
  /// Dirichlet BC
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
  model.addDumpField("damage");
  model.addDumpField("partitions");
  model.addDumpField("material_index");
  model.addDumpField("Sc");
  model.addDumpField("force");
  model.addDumpField("equivalent_stress");

  model.dump();

  std::stringstream error_stream;
  error_stream << "error"
               << ".csv";
  std::ofstream error_output;
  error_output.open(error_stream.str().c_str());
  error_output << "# Step, Average, Max, Min" << std::endl;

  checkDisplacement(model, element_type, error_output, 0, true);

  MaterialDamageIterative<spatial_dimension> & aggregate =
      dynamic_cast<MaterialDamageIterative<spatial_dimension> &>(
          model.getMaterial(0));
  MaterialDamageIterative<spatial_dimension> & paste =
      dynamic_cast<MaterialDamageIterative<spatial_dimension> &>(
          model.getMaterial(1));

  Real error;
  bool converged = false;
  UInt nb_damaged_elements = 0;
  Real max_eq_stress_agg = 0;
  Real max_eq_stress_paste = 0;

  /// solve the system
  converged =
      model.solveStep<_scm_newton_raphson_tangent_modified,
                      SolveConvergenceCriteria::_increment>(1e-12, error, 2);

  if (converged == false) {
    std::cout << "The error is: " << error << std::endl;
    AKANTU_DEBUG_ASSERT(converged, "Did not converge");
  }

  if (!checkDisplacement(model, element_type, error_output, 1, false)) {
    finalize();
    return EXIT_FAILURE;
  }

  model.dump();

  /// get the maximum equivalent stress in both materials
  max_eq_stress_agg = aggregate.getNormMaxEquivalentStress();
  max_eq_stress_paste = paste.getNormMaxEquivalentStress();

  nb_damaged_elements = 0;
  if (max_eq_stress_agg > max_eq_stress_paste)
    nb_damaged_elements = aggregate.updateDamage();
  else
    nb_damaged_elements = paste.updateDamage();

  if (prank == 0 && nb_damaged_elements)
    std::cout << nb_damaged_elements << " elements damaged" << std::endl;

  /// resolve the system
  converged =
      model.solveStep<_scm_newton_raphson_tangent_modified,
                      SolveConvergenceCriteria::_increment>(1e-12, error, 2);

  if (converged == false) {
    std::cout << "The error is: " << error << std::endl;
    AKANTU_DEBUG_ASSERT(converged, "Did not converge");
  }

  if (!checkDisplacement(model, element_type, error_output, 2, false)) {
    finalize();
    return EXIT_FAILURE;
  }

  model.dump();

  finalize();

  return EXIT_SUCCESS;
}

/* -------------------------------------------------------------------------- */
bool checkDisplacement(SolidMechanicsModel & model, ElementType type,
                       std::ofstream & error_output, UInt step,
                       bool barycenters) {

  Mesh & mesh = model.getMesh();
  UInt spatial_dimension = mesh.getSpatialDimension();
  const Array<UInt> & connectivity = mesh.getConnectivity(type);
  const Array<Real> & displacement = model.getDisplacement();
  UInt nb_element = mesh.getNbElement(type);
  UInt nb_nodes_per_elem = Mesh::getNbNodesPerElement(type);

  StaticCommunicator & comm = StaticCommunicator::getStaticCommunicator();
  Int psize = comm.getNbProc();
  Int prank = comm.whoAmI();

  if (psize == 1) {
    std::stringstream displacement_file;
    displacement_file << "displacement/displacement_" << std::setfill('0')
                      << std::setw(6) << step;
    std::ofstream displacement_output;
    displacement_output.open(displacement_file.str().c_str());

    for (UInt el = 0; el < nb_element; ++el) {
      for (UInt n = 0; n < nb_nodes_per_elem; ++n) {
        UInt node = connectivity(el, n);

        for (UInt dim = 0; dim < spatial_dimension; ++dim) {
          displacement_output << std::setprecision(15)
                              << displacement(node, dim) << " ";
        }
        displacement_output << std::endl;
      }
    }

    displacement_output.close();

    if (barycenters) {
      std::stringstream barycenter_file;
      barycenter_file << "displacement/barycenters";
      std::ofstream barycenter_output;
      barycenter_output.open(barycenter_file.str().c_str());

      Element element(type, 0);
      Vector<Real> bary(spatial_dimension);

      for (UInt el = 0; el < nb_element; ++el) {
        element.element = el;
        mesh.getBarycenter(element, bary);

        for (UInt dim = 0; dim < spatial_dimension; ++dim) {
          barycenter_output << std::setprecision(15) << bary(dim) << " ";
        }
        barycenter_output << std::endl;
      }

      barycenter_output.close();
    }
  } else {

    if (barycenters)
      return true;

    /// read data
    std::stringstream displacement_file;
    displacement_file << "displacement/displacement_" << std::setfill('0')
                      << std::setw(6) << step;
    std::ifstream displacement_input;
    displacement_input.open(displacement_file.str().c_str());

    Array<Real> displacement_serial(0, spatial_dimension);
    Vector<Real> disp_tmp(spatial_dimension);

    while (displacement_input.good()) {
      for (UInt i = 0; i < spatial_dimension; ++i)
        displacement_input >> disp_tmp(i);

      displacement_serial.push_back(disp_tmp);
    }

    std::stringstream barycenter_file;
    barycenter_file << "displacement/barycenters";
    std::ifstream barycenter_input;
    barycenter_input.open(barycenter_file.str().c_str());

    Array<Real> barycenter_serial(0, spatial_dimension);

    while (barycenter_input.good()) {
      for (UInt dim = 0; dim < spatial_dimension; ++dim)
        barycenter_input >> disp_tmp(dim);

      barycenter_serial.push_back(disp_tmp);
    }

    Element element(type, 0);
    Vector<Real> bary(spatial_dimension);

    Array<Real>::iterator<Vector<Real>> it;
    Array<Real>::iterator<Vector<Real>> begin =
        barycenter_serial.begin(spatial_dimension);
    Array<Real>::iterator<Vector<Real>> end =
        barycenter_serial.end(spatial_dimension);

    Array<Real>::const_iterator<Vector<Real>> disp_it;
    Array<Real>::iterator<Vector<Real>> disp_serial_it;

    Vector<Real> difference(spatial_dimension);
    Array<Real> error;

    /// compute error
    for (UInt el = 0; el < nb_element; ++el) {
      element.element = el;
      mesh.getBarycenter(element, bary);

      /// find element
      for (it = begin; it != end; ++it) {
        UInt matched_dim = 0;

        while (matched_dim < spatial_dimension &&
               Math::are_float_equal(bary(matched_dim), (*it)(matched_dim)))
          ++matched_dim;

        if (matched_dim == spatial_dimension)
          break;
      }

      if (it == end) {
        std::cout << "Element barycenter not found!" << std::endl;
        return false;
      }

      UInt matched_el = it - begin;

      disp_serial_it = displacement_serial.begin(spatial_dimension) +
                       matched_el * nb_nodes_per_elem;

      for (UInt n = 0; n < nb_nodes_per_elem; ++n, ++disp_serial_it) {
        UInt node = connectivity(el, n);
        if (!mesh.isLocalOrMasterNode(node))
          continue;

        disp_it = displacement.begin(spatial_dimension) + node;

        difference = *disp_it;
        difference -= *disp_serial_it;

        error.push_back(difference.norm());
      }
    }

    /// compute average error
    Real average_error = std::accumulate(error.begin(), error.end(), 0.);
    comm.allReduce(&average_error, 1, _so_sum);

    UInt error_size = error.getSize();
    comm.allReduce(&error_size, 1, _so_sum);

    average_error /= error_size;

    /// compute maximum and minimum
    Real max_error = *std::max_element(error.begin(), error.end());
    comm.allReduce(&max_error, 1, _so_max);

    Real min_error = *std::min_element(error.begin(), error.end());
    comm.allReduce(&min_error, 1, _so_min);

    /// output data
    if (prank == 0) {
      error_output << step << ", " << average_error << ", " << max_error << ", "
                   << min_error << std::endl;
    }

    if (max_error > 1.e-9) {
      std::cout << "Displacement error of " << max_error << " is too big!"
                << std::endl;
      return false;
    }
  }

  return true;
}

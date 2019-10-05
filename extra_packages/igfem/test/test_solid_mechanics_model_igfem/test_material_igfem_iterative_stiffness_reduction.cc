/**
 * @file   test_material_igfem_iterative_strength_reduction.cc
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @date   Thu Nov 26 12:20:15 2015
 *
 * @brief  test the material iterative stiffness reduction
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
#include "material_damage_iterative.hh"
#include "material_igfem_saw_tooth_damage.hh"
#include "solid_mechanics_model_igfem.hh"
/* -------------------------------------------------------------------------- */
using namespace akantu;

/// function declaration
bool checkDamageState(UInt step, const SolidMechanicsModelIGFEM & model,
                      bool igfem_analysis);

class TestMaterialSelector : public MaterialSelector {
public:
  TestMaterialSelector(SolidMechanicsModelIGFEM & model)
      : MaterialSelector(), model(model),
        spatial_dimension(model.getSpatialDimension()) {}

  UInt operator()(const Element & element) {
    if (Mesh::getKind(element.type) == _ek_igfem)
      return 2;
    else {
      /// regular elements
      const Mesh & mesh = model.getMesh();
      Vector<Real> barycenter(this->spatial_dimension);
      mesh.getBarycenter(element, barycenter);
      /// check if element belongs to ASR gel
      if (model.isInside(barycenter))
        return 1;
    }
    return 0;
  }

protected:
  SolidMechanicsModelIGFEM & model;
  UInt spatial_dimension;
};

/* -------------------------------------------------------------------------- */
/* Main                                                                       */
/* -------------------------------------------------------------------------- */
int main(int argc, char * argv[]) {

  Math::setTolerance(1e-13);
  debug::setDebugLevel(dblWarning);

  initialize("material_stiffness_reduction.dat", argc, argv);

  bool igfem_analysis;
  std::string action(argv[1]);
  if (action == "igfem") {
    igfem_analysis = true;
  } else if (action == "standard_fem") {
    igfem_analysis = false;
  } else {
    std::cerr << "invalid option" << std::endl;
  }

  const UInt spatial_dimension = 2;
  StaticCommunicator & comm =
      akantu::StaticCommunicator::getStaticCommunicator();
  Int psize = comm.getNbProc();
  Int prank = comm.whoAmI();

  /// read the mesh and partion it
  Mesh mesh(spatial_dimension);
  akantu::MeshPartition * partition = NULL;

  if (prank == 0) {

    if (igfem_analysis)
      mesh.read("igfem_mesh.msh");
    else
      mesh.read("regular_mesh.msh");

    /// partition the mesh
    partition = new MeshPartitionScotch(mesh, spatial_dimension);

    partition->partitionate(psize);
  }

  /// model creation
  SolidMechanicsModelIGFEM model(mesh);
  model.initParallel(partition);
  delete partition;

  Math::setTolerance(1.e-14);
  /// intialize the geometry and set the material selector
  std::list<SK::Sphere_3> inclusions_list;
  model.registerGeometryObject(inclusions_list, "inclusion");
  Real val = 1000000000;
  Real radius_squared = val * val;
  Vector<Real> center(spatial_dimension);
  center(0) = 0;
  center(1) = val;
  SK::Sphere_3 sphere(SK::Point_3(center(0), center(1), 0.), radius_squared);
  inclusions_list.push_back(sphere);
  TestMaterialSelector * mat_selector = new TestMaterialSelector(model);
  model.setMaterialSelector(*mat_selector);

  /// initialization of the model
  model.initFull();

  /// create the interface
  if (igfem_analysis)
    model.update("inclusion");

  /// boundary conditions
  mesh.computeBoundingBox();
  const Vector<Real> & lowerBounds = mesh.getLowerBounds();
  const Vector<Real> & upperBounds = mesh.getUpperBounds();
  Real bottom = lowerBounds(1);
  Real top = upperBounds(1);
  Real left = lowerBounds(0);
  Real eps = std::abs((top - bottom) * 1e-6);
  const Array<Real> & pos = mesh.getNodes();
  Array<bool> & boun = model.getBlockedDOFs();
  Array<Real> & disp = model.getDisplacement();
  for (UInt n = 0; n < mesh.getNbNodes(); ++n) {
    if (std::abs(pos(n, 1) - bottom) < eps) {
      boun(n, 1) = true;
      disp(n, 1) = 0.;
    }
    if (std::abs(pos(n, 1) - top) < eps) {
      boun(n, 1) = true;
      disp(n, 1) = 1.e-3;
    }
    if (std::abs(pos(n, 0) - left) < eps) {
      boun(n, 0) = true;
      disp(n, 0) = 0.;
    }
  }

  /// add fields that should be dumped
  model.setBaseName("regular");
  model.addDumpField("material_index");
  model.addDumpFieldVector("displacement");
  ;
  model.addDumpField("stress");
  model.addDumpField("blocked_dofs");
  model.addDumpField("residual");
  model.addDumpField("grad_u");
  model.addDumpField("damage");
  model.addDumpField("partitions");
  model.addDumpField("Sc");
  model.addDumpField("force");
  model.addDumpField("equivalent_stress");
  model.addDumpField("ultimate_strain");
  model.setBaseNameToDumper("igfem elements", "igfem elements");
  model.addDumpFieldToDumper("igfem elements", "material_index");
  model.addDumpFieldVectorToDumper("igfem elements", "displacement");
  ;
  model.addDumpFieldToDumper("igfem elements", "stress");
  model.addDumpFieldToDumper("igfem elements", "blocked_dofs");
  model.addDumpFieldToDumper("igfem elements", "residual");
  model.addDumpFieldToDumper("igfem elements", "grad_u");
  model.addDumpFieldToDumper("igfem elements", "damage");
  model.addDumpFieldToDumper("igfem elements", "partitions");
  model.addDumpFieldToDumper("igfem elements", "Sc");
  model.addDumpFieldToDumper("igfem elements", "force");
  model.addDumpFieldToDumper("igfem elements", "equivalent_stress");
  model.addDumpFieldToDumper("igfem elements", "ultimate_strain");

  model.dump();
  model.dump("igfem elements");

  /// get a reference to the damage materials
  MaterialDamageIterative<spatial_dimension> & material =
      dynamic_cast<MaterialDamageIterative<spatial_dimension> &>(
          model.getMaterial(0));
  MaterialIGFEMSawToothDamage<spatial_dimension> & igfem_material =
      dynamic_cast<MaterialIGFEMSawToothDamage<spatial_dimension> &>(
          model.getMaterial(2));

  Real error;
  bool converged = false;
  UInt nb_damaged_elements = 0;
  Real max_eq_stress_regular = 0;
  Real max_eq_stress_igfem = 0;

  /// solve the system
  // counter for the damage steps
  UInt s = 0;
  do {
    converged =
        model.solveStep<_scm_newton_raphson_tangent_modified,
                        SolveConvergenceCriteria::_increment>(1e-12, error, 2);

    if (converged == false) {
      std::cout << "The error is: " << error << std::endl;
      AKANTU_DEBUG_ASSERT(converged, "Did not converge");
    }

    /// compute damage
    max_eq_stress_regular = material.getNormMaxEquivalentStress();
    max_eq_stress_igfem = igfem_material.getNormMaxEquivalentStress();
    if (max_eq_stress_regular > max_eq_stress_igfem)
      nb_damaged_elements = material.updateDamage();
    else if (max_eq_stress_regular == max_eq_stress_igfem) {
      nb_damaged_elements = material.updateDamage();
      nb_damaged_elements += igfem_material.updateDamage();
    } else
      nb_damaged_elements = igfem_material.updateDamage();

    model.dump();
    model.dump("igfem elements");

    /// check the current damage state
    if (!checkDamageState(s, model, igfem_analysis)) {
      std::cout << "error in the damage compuation" << std::endl;
      finalize();
      return EXIT_FAILURE;
    }

    s++;
  } while (nb_damaged_elements);

  std::cout << action << " passed!!" << std::endl;
  finalize();

  return EXIT_SUCCESS;
}

/* -------------------------------------------------------------------------- */
bool checkDamageState(UInt step, const SolidMechanicsModelIGFEM & model,
                      bool igfem_analysis) {

  bool test_result = true;
  const UInt spatial_dimension = model.getSpatialDimension();
  const Mesh & mesh = model.getMesh();

  if (!igfem_analysis) {
    const ElementType element_type = _triangle_3;
    /// prepare output: compute barycenters for elements that can be damaged
    const Array<UInt> & element_filter =
        model.getMaterial(0).getElementFilter(element_type, _not_ghost);
    Array<Real> barycenters(element_filter.getSize(), spatial_dimension);
    Array<Real>::vector_iterator bary_it = barycenters.begin(spatial_dimension);
    for (UInt e = 0; e < element_filter.getSize(); ++e, ++bary_it) {
      UInt global_el_idx = element_filter(e);
      mesh.getBarycenter(global_el_idx, element_type, bary_it->storage(),
                         _not_ghost);
    }

    const Array<Real> & damage = model.getMaterial(0).getInternal<Real>(
        "damage")(element_type, _not_ghost);
    const Array<Real> & Sc =
        model.getMaterial(0).getInternal<Real>("Sc")(element_type, _not_ghost);

    std::ostringstream file_name;
    file_name << "step_" << std::setfill('0') << std::setw(3) << step << ".txt";

    std::ofstream file_output;
    file_output.open(file_name.str());
    file_output << std::setprecision(14);

    for (UInt e = 0; e < barycenters.getSize(); ++e)
      file_output << barycenters(e, 0) << " " << barycenters(e, 1) << " "
                  << damage(e) << " " << Sc(e) << std::endl;

  }

  else {

    /// read data
    Real default_tolerance = Math::getTolerance();
    Math::setTolerance(1.e-10);
    std::stringstream results_file;
    results_file << "step_" << std::setfill('0') << std::setw(3) << step
                 << ".txt";
    std::ifstream damage_input;
    damage_input.open(results_file.str().c_str());

    Array<Real> damage_result(0, 1);
    Array<Real> Sc_result(0, 1);
    Array<Real> bary_regular(0, spatial_dimension);
    Vector<Real> bary(spatial_dimension);
    Real dam = 0.;
    Real strength = 0;

    while (damage_input.good()) {
      damage_input >> bary(0) >> bary(1) >> dam >> strength;
      bary_regular.push_back(bary);
      damage_result.push_back(dam);
      Sc_result.push_back(strength);
    }

    /// compare the results
    Array<Real>::const_vector_iterator bary_it;

    Array<Real>::const_vector_iterator bary_begin =
        bary_regular.begin(spatial_dimension);
    Array<Real>::const_vector_iterator bary_end =
        bary_regular.end(spatial_dimension);
    /// compare the regular elements
    ElementType element_type = _triangle_3;
    const Array<UInt> & element_filter =
        model.getMaterial(0).getElementFilter(element_type, _not_ghost);
    const Array<Real> & damage_regular_el =
        model.getMaterial(0).getInternal<Real>("damage")(element_type,
                                                         _not_ghost);
    const Array<Real> & Sc_regular_el =
        model.getMaterial(0).getInternal<Real>("Sc")(element_type, _not_ghost);

    for (UInt e = 0; e < element_filter.getSize(); ++e) {
      UInt global_el_idx = element_filter(e);
      mesh.getBarycenter(global_el_idx, element_type, bary.storage(),
                         _not_ghost);
      /// find element
      for (bary_it = bary_begin; bary_it != bary_end; ++bary_it) {
        UInt matched_dim = 0;
        while (
            matched_dim < spatial_dimension &&
            Math::are_float_equal(bary(matched_dim), (*bary_it)(matched_dim)))
          ++matched_dim;
        if (matched_dim == spatial_dimension)
          break;
      }
      if (bary_it == bary_end) {
        std::cout << "Element barycenter not found!" << std::endl;
        return false;
      }

      UInt matched_el = bary_it - bary_begin;

      if (std::abs(damage_result(matched_el) - damage_regular_el(e)) > 1.e-12 ||
          std::abs(Sc_result(matched_el) - Sc_regular_el(e)) > 1.e-4) {
        test_result = false;
        break;
      }
    }
    /// compare the IGFEM elements
    UInt nb_sub_elements = 2;
    element_type = _igfem_triangle_4;
    const Array<UInt> & element_filter_igfem =
        model.getMaterial(2).getElementFilter(element_type, _not_ghost);
    const Array<Real> & damage_regular_el_igfem =
        model.getMaterial(2).getInternal<Real>("damage")(element_type,
                                                         _not_ghost);
    const Array<Real> & Sc_regular_el_igfem =
        model.getMaterial(2).getInternal<Real>("Sc")(element_type, _not_ghost);
    UInt * sub_el_ptr =
        model.getMaterial(2)
            .getInternal<UInt>("sub_material")(element_type, _not_ghost)
            .storage();

    for (UInt e = 0; e < element_filter_igfem.getSize(); ++e) {
      UInt global_el_idx = element_filter_igfem(e);
      for (UInt s = 0; s < nb_sub_elements; ++s, ++sub_el_ptr) {
        if (*sub_el_ptr)
          model.getSubElementBarycenter(global_el_idx, s, element_type, bary,
                                        _not_ghost);
        else
          continue;
        /// find element
        for (bary_it = bary_begin; bary_it != bary_end; ++bary_it) {
          UInt matched_dim = 0;
          while (
              matched_dim < spatial_dimension &&
              Math::are_float_equal(bary(matched_dim), (*bary_it)(matched_dim)))
            ++matched_dim;
          if (matched_dim == spatial_dimension)
            break;
        }
        if (bary_it == bary_end) {
          std::cout << "Element barycenter not found!" << std::endl;
          return false;
        }

        UInt matched_el = bary_it - bary_begin;

        if (std::abs(damage_result(matched_el) -
                     damage_regular_el_igfem(e * nb_sub_elements + s)) >
                1.e-12 ||
            std::abs(Sc_result(matched_el) -
                     Sc_regular_el_igfem(e * nb_sub_elements + s)) > 1.e-4) {
          test_result = false;
          break;
        }
      }
    }

    Math::setTolerance(default_tolerance);
  }
  return test_result;
}

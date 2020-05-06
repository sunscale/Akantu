/**
 * @file
 * test_material_igfem_iterative_stiffness_reduction_damage_step_transfer.cc
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @date   Thu Nov 26 12:20:15 2015
 *
 * @brief test the damage step transfer for the material iterative
 * stiffness reduction
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
#include "material_igfem_saw_tooth_damage.hh"
#include "material_iterative_stiffness_reduction.hh"
#include "solid_mechanics_model_igfem.hh"
/* -------------------------------------------------------------------------- */
using namespace akantu;

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

  initialize("material_stiffness_reduction_2.dat", argc, argv);

  const UInt spatial_dimension = 2;
  StaticCommunicator & comm =
      akantu::StaticCommunicator::getStaticCommunicator();
  Int psize = comm.getNbProc();
  Int prank = comm.whoAmI();

  /// read the mesh and partion it
  Mesh mesh(spatial_dimension);
  akantu::MeshPartition * partition = NULL;

  if (prank == 0) {

    mesh.read("test_damage_transfer.msh");

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
  Real radius_squared = (val - 0.6) * (val - 0.6);
  Vector<Real> center(spatial_dimension);
  center(0) = 0;
  center(1) = val;
  SK::Sphere_3 sphere(SK::Point_3(center(0), center(1), 0.), radius_squared);
  inclusions_list.push_back(sphere);
  TestMaterialSelector * mat_selector = new TestMaterialSelector(model);
  model.setMaterialSelector(*mat_selector);

  /// initialization of the model
  model.initFull();

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
  MaterialIterativeStiffnessReduction<spatial_dimension> & material =
      dynamic_cast<MaterialIterativeStiffnessReduction<spatial_dimension> &>(
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
  UInt regular_steps = 15;
  for (UInt s = 0; s < regular_steps; ++s) {
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

    if (!nb_damaged_elements)
      break;
    model.dump();
    model.dump("igfem elements");
  }

  const Array<UInt> & reduction_step_regular =
      material.getInternal<UInt>("damage_step")(_triangle_3, _not_ghost);
  UInt reduction_step_el_27 = reduction_step_regular(27);
  UInt reduction_step_el_19 = reduction_step_regular(19);

  /// create the interface
  Real new_radius = (val - 0.1);
  model.moveInterface(new_radius);
  model.dump();
  model.dump("igfem elements");

  /// check that the damage reduction step has been correctly computed
  /// regular element id -> igfem element id
  /// 27 -> 7; 19 -> 5
  const Array<UInt> & reduction_step_igfem = igfem_material.getInternal<UInt>(
      "damage_step")(_igfem_triangle_5, _not_ghost);
  Array<UInt>::const_scalar_iterator step_it = reduction_step_igfem.begin();

  /// check the igfem elements
  UInt nb_igfem_elements = mesh.getNbElement(_igfem_triangle_5, _not_ghost);
  UInt nb_quads = model.getFEEngine("IGFEMFEEngine")
                      .getNbIntegrationPoints(_igfem_triangle_5, _not_ghost);
  const Array<UInt> & sub_material = igfem_material.getInternal<UInt>(
      "sub_material")(_igfem_triangle_5, _not_ghost);
  Array<UInt>::const_scalar_iterator sub_it = sub_material.begin();
  for (UInt e = 0; e < nb_igfem_elements; ++e) {
    for (UInt q = 0; q < nb_quads; ++q, ++sub_it, ++step_it) {
      if (!*sub_it) {
        if (!Math::are_float_equal(*step_it, 0.)) {
          std::cout
              << "the reduction step for an elastic sub-element must be zero!!"
              << std::endl;
          finalize();
          return EXIT_FAILURE;
        }
      } else {
        if (e == 7) {
          if (!Math::are_float_equal(*step_it, reduction_step_el_27)) {
            std::cout << "error in computation of damage step!!" << std::endl;
            finalize();
            return EXIT_FAILURE;
          }
        } else if (e == 5) {
          if (!Math::are_float_equal(*step_it, reduction_step_el_19)) {
            std::cout << "error in computation of damage step!!" << std::endl;
            finalize();
            return EXIT_FAILURE;
          }
        } else {
          if (!Math::are_float_equal(*step_it, 0.)) {
            std::cout << "error in computation of damage step!!" << std::endl;
            finalize();
            return EXIT_FAILURE;
          }
        }
      }
    }
  }

  //// force the next damage event
  const Array<Real> & dam_igfem =
      igfem_material.getInternal<Real>("damage")(_igfem_triangle_5, _not_ghost);
  Array<Real> old_damage(dam_igfem);

  for (UInt s = 0; s < 1; ++s) {
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

    if (!nb_damaged_elements)
      break;
    model.dump();
    model.dump("igfem elements");
  }

  /// check that damage has been simultanously been updated on all the
  /// the integration points of one sub-element
  const Array<Real> & new_dam_igfem =
      igfem_material.getInternal<Real>("damage")(_igfem_triangle_5, _not_ghost);
  sub_it = sub_material.begin();
  Array<Real>::const_scalar_iterator new_dam_it = new_dam_igfem.begin();
  Array<Real>::const_scalar_iterator old_dam_it = old_damage.begin();
  step_it = reduction_step_igfem.begin();
  UInt reduction_constant = material.getParam<Real>("reduction_constant");

  for (UInt e = 0; e < nb_igfem_elements; ++e) {
    for (UInt q = 0; q < nb_quads;
         ++q, ++sub_it, ++step_it, ++new_dam_it, ++old_dam_it) {
      if (!*sub_it) {
        if (!Math::are_float_equal(*step_it, 0.) ||
            !Math::are_float_equal(*new_dam_it, 0.)) {
          std::cout << "the reduction step and damagefor an elastic "
                       "sub-element must be zero!!"
                    << std::endl;
          finalize();
          return EXIT_FAILURE;
        }
      } else {
        if (e == 7) {
          if (!Math::are_float_equal(*step_it, reduction_step_el_27 + 1) ||
              !Math::are_float_equal(
                  *new_dam_it, 1 - (1. / std::pow(reduction_constant,
                                                  reduction_step_el_27 + 1)))) {
            std::cout << "error in computation of damage step!!" << std::endl;
            finalize();
            return EXIT_FAILURE;
          }
        } else if (e == 5) {
          if (!Math::are_float_equal(*step_it, reduction_step_el_19) ||
              !Math::are_float_equal(*new_dam_it, *old_dam_it)) {
            std::cout << "error in computation of damage step!!" << std::endl;
            finalize();
            return EXIT_FAILURE;
          }
        } else {
          if (!Math::are_float_equal(*step_it, 0.) ||
              !Math::are_float_equal(*new_dam_it, 0.)) {
            std::cout << "error in computation of damage step!!" << std::endl;
            finalize();
            return EXIT_FAILURE;
          }
        }
      }
    }
  }

  finalize();

  return EXIT_SUCCESS;
}

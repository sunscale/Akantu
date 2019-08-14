/**
 * @file   test_ASR_damage.cc
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 *
 *
 * @brief  test the solidmechancis model for IGFEM analysis
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "solid_mechanics_model_igfem.hh"
/* -------------------------------------------------------------------------- */
#include "material_damage_iterative.hh"
#include "material_igfem_saw_tooth_damage.hh"
/* -------------------------------------------------------------------------- */
using namespace akantu;

/// function declaration
void applyBoundaryConditions(SolidMechanicsModelIGFEM & model);

class ASRMaterialSelector : public MaterialSelector {
public:
  ASRMaterialSelector(SolidMechanicsModelIGFEM & model) : model(model) {}

  UInt operator()(const Element & elem) {
    if (Mesh::getKind(elem.type) == _ek_igfem)
      /// choose IGFEM material
      return 2;

    const Mesh & mesh = model.getMesh();
    UInt spatial_dimension = model.getSpatialDimension();
    Vector<Real> barycenter(spatial_dimension);
    mesh.getBarycenter(elem, barycenter);
    if (model.isInside(barycenter))
      return 1;
    return 0;
  }

protected:
  SolidMechanicsModelIGFEM & model;
};

typedef Spherical SK;

int main(int argc, char * argv[]) {

  initialize("material_ASR.dat", argc, argv);

  /// problem dimension
  const UInt spatial_dimension = 2;
  StaticCommunicator & comm =
      akantu::StaticCommunicator::getStaticCommunicator();
  Int psize = comm.getNbProc();
  Int prank = comm.whoAmI();
  /// mesh creation
  Mesh mesh(spatial_dimension);
  akantu::MeshPartition * partition = NULL;
  if (prank == 0) {
    mesh.read("one_inclusion.msh");
    /// partition the mesh
    partition = new MeshPartitionScotch(mesh, spatial_dimension);
    partition->partitionate(psize);
  }

  /// model creation
  SolidMechanicsModelIGFEM model(mesh);
  model.initParallel(partition);
  delete partition;
  /// register the gel pocket list in the model
  std::list<SK::Sphere_3> gel_pocket_list;
  model.registerGeometryObject(gel_pocket_list, "gel");

  ASRMaterialSelector * mat_selector;

  mat_selector = new ASRMaterialSelector(model);
  model.setMaterialSelector(*mat_selector);
  model.initFull();

  /// add fields that should be dumped
  model.setBaseName("regular_elements");
  model.addDumpField("material_index");
  model.addDumpField("damage");
  model.addDumpField("Sc");
  model.addDumpField("partitions");
  model.addDumpField("eigen_grad_u");
  model.addDumpField("blocked_dofs");
  model.setBaseNameToDumper("igfem elements", "igfem elements");
  model.addDumpFieldToDumper("igfem elements", "material_index");
  model.addDumpFieldToDumper("igfem elements", "Sc");
  model.addDumpFieldToDumper("igfem elements", "damage");
  model.addDumpFieldToDumper("igfem elements", "lambda");
  model.addDumpFieldToDumper("igfem elements", "eigen_grad_u");
  model.addDumpFieldToDumper("igfem elements", "blocked_dofs");

  /// dump before the interface generation
  model.dump();
  model.dump("igfem elements");

  /// weaken one element to enforce damage there
  Array<Real> & Sc =
      model.getMaterial(0).getInternal<Real>("Sc")(_triangle_3, _not_ghost);
  Sc(11) = 1;
  /// create the gel pocket
  Real initial_gel_radius = 0.1;
  SK::Sphere_3 sphere_1(SK::Point_3(0., 0., 0.),
                        initial_gel_radius * initial_gel_radius);
  gel_pocket_list.push_back(sphere_1);

  /// create the interface
  model.update("gel");

  ///  apply eigenstrain the eigenstrain in the inclusions
  Matrix<Real> prestrain(spatial_dimension, spatial_dimension, 0.);
  for (UInt i = 0; i < spatial_dimension; ++i)
    prestrain(i, i) = 0.05;

  model.applyEigenGradU(prestrain, "gel", _not_ghost);
  applyBoundaryConditions(model);
  /// dump
  model.dump("igfem elements");
  model.dump();

  /// Instantiate objects of class MyDamageso
  MaterialDamageIterative<spatial_dimension> & mat_aggregate =
      dynamic_cast<MaterialDamageIterative<spatial_dimension> &>(
          model.getMaterial(0));
  MaterialIGFEMSawToothDamage<spatial_dimension> & mat_igfem =
      dynamic_cast<MaterialIGFEMSawToothDamage<spatial_dimension> &>(
          model.getMaterial(2));

  bool factorize = false;
  bool converged = false;
  Real error;

  UInt nb_damaged_elements = 0;
  Real max_eq_stress_aggregate = 0;
  Real max_igfem = 0;

  const Array<Real> & stress =
      model.getMaterial(2).getStress(_igfem_triangle_5, _not_ghost);
  Array<Real>::const_matrix_iterator stress_it =
      stress.begin(spatial_dimension, spatial_dimension);
  do {
    converged = model.solveStep<_scm_newton_raphson_tangent,
                                SolveConvergenceCriteria::_increment>(
        1e-6, error, 2, factorize);

    /// compute damage
    max_eq_stress_aggregate = mat_aggregate.getNormMaxEquivalentStress();
    max_igfem = mat_igfem.getNormMaxEquivalentStress();
    if (max_eq_stress_aggregate > max_igfem)
      nb_damaged_elements = mat_aggregate.updateDamage();
    else
      nb_damaged_elements = mat_igfem.updateDamage();
    std::cout << "damaged elements: " << nb_damaged_elements << std::endl;
    for (UInt i = 0; i < 5; ++i) {
      std::cout << *stress_it << std::endl;
      ++stress_it;
    }
    model.dump();
    model.dump("igfem elements");
  } while (nb_damaged_elements);

  model.dump();

  model.dump("igfem elements");

  finalize();
  return EXIT_SUCCESS;
}

/* -------------------------------------------------------------------------- */
void applyBoundaryConditions(SolidMechanicsModelIGFEM & model) {
  /// boundary conditions
  Mesh & mesh = model.getMesh();
  mesh.computeBoundingBox();
  const Vector<Real> & lowerBounds = mesh.getLowerBounds();
  const Vector<Real> & upperBounds = mesh.getUpperBounds();
  Real bottom = lowerBounds(1);
  Real top = upperBounds(1);
  Real left = lowerBounds(0);
  //  Real right = upperBounds(0);

  Real eps = std::abs((top - bottom) * 1e-12);
  const Array<Real> & pos = mesh.getNodes();
  Array<Real> & disp = model.getDisplacement();
  Array<bool> & boun = model.getBlockedDOFs();

  disp.clear();
  boun.clear();
  /// free expansion
  for (UInt i = 0; i < mesh.getNbNodes(); ++i) {

    if (std::abs(pos(i, 1) - bottom) < eps) {
      boun(i, 1) = true;
      disp(i, 1) = 0.0;
    }

    if (std::abs(pos(i, 0) - left) < eps) {
      boun(i, 0) = true;
      disp(i, 0) = 0.0;
    }
  }
}

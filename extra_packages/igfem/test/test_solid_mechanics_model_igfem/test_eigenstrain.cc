/**
 * @file   test_eigenstrain.cc
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 *
 *
 * @brief  test to that eigenstrain is only applied on one sub-material
 *
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "dumper_paraview.hh"
#include "material_igfem.hh"
#include "mesh_geom_common.hh"
#include "solid_mechanics_model_igfem.hh"
/* -------------------------------------------------------------------------- */
using namespace akantu;

class ASRMaterialSelector : public DefaultMaterialIGFEMSelector {
public:
  ASRMaterialSelector(SolidMechanicsModelIGFEM & model)
      : DefaultMaterialIGFEMSelector(model), model(model) {}

  UInt operator()(const Element & elem) {
    if (Mesh::getKind(elem.type) == _ek_igfem)
      /// choose IGFEM material
      return this->fallback_value_igfem;

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

  initialize("material_damage.dat", argc, argv);
  StaticCommunicator & comm =
      akantu::StaticCommunicator::getStaticCommunicator();
  Int psize = comm.getNbProc();
  Int prank = comm.whoAmI();

  /// problem dimension
  UInt spatial_dimension = 2;

  /// mesh creation and partioning
  Mesh mesh(spatial_dimension);
  akantu::MeshPartition * partition = NULL;

  if (prank == 0) {
    mesh.read("fine_mesh.msh");
    /// partition the mesh
    partition = new MeshPartitionScotch(mesh, spatial_dimension);
    partition->partitionate(psize);
  }

  /// model creation and initialization
  SolidMechanicsModelIGFEM model(mesh);
  model.initParallel(partition);
  delete partition;

  /// create the list to store the gel pockets
  std::list<SK::Sphere_3> gel_pocket_list;
  model.registerGeometryObject(gel_pocket_list, "mat_1");
  /// set the material selector
  ASRMaterialSelector * mat_selector = new ASRMaterialSelector(model);
  model.setMaterialSelector(*mat_selector);
  model.initFull();

  /// add fields that should be dumped
  model.setBaseName("regular_elements");
  model.addDumpField("material_index");
  model.addDumpField("damage");
  model.addDumpField("Sc");
  model.addDumpField("partitions");
  model.addDumpField("eigen_grad_u");
  model.setBaseNameToDumper("igfem elements", "igfem elements");
  model.addDumpFieldToDumper("igfem elements", "material_index");
  model.addDumpFieldToDumper("igfem elements", "Sc");
  model.addDumpFieldToDumper("igfem elements", "lambda");
  model.addDumpFieldToDumper("igfem elements", "partitions");
  model.addDumpFieldToDumper("igfem elements", "eigen_grad_u");

  /// dump
  model.dump("igfem elements");
  model.dump();

  /// create the inclusions
  SK::Sphere_3 sphere_1(SK::Point_3(0., 0., 0.), 0.13 * 0.13);
  SK::Sphere_3 sphere_2(SK::Point_3(0.5, 0.5, 0.), 0.4 * 0.4);
  SK::Sphere_3 sphere_3(SK::Point_3(-0.75, -0.75, 0.), 0.12 * 0.12);
  SK::Sphere_3 sphere_4(SK::Point_3(0.625, -0.625, 0.), 0.25 * 0.25);
  gel_pocket_list.push_back(sphere_1);
  gel_pocket_list.push_back(sphere_2);
  gel_pocket_list.push_back(sphere_3);
  gel_pocket_list.push_back(sphere_4);

  /// create the interface
  model.update("mat_1");
  if (mesh.getNbElement(_igfem_triangle_4, _not_ghost)) {
    /// something went wrong in the interface creation
    finalize();
    return EXIT_FAILURE;
  }

  ///  apply eigenstrain the eigenstrain in the inclusions
  Matrix<Real> prestrain(spatial_dimension, spatial_dimension, 0.);
  for (UInt i = 0; i < spatial_dimension; ++i)
    prestrain(i, i) = 0.07;

  model.applyEigenGradU(prestrain, "mat_1", _not_ghost);

  /// check that eigenstrain has been applied correctly
  /// check first the regular materials (the first two in the mat file)
  Real error = 0;
  Material & mat_1 = model.getMaterial(0);
  const Array<Real> & eigen_grad_u_1 =
      mat_1.getArray<Real>("eigen_grad_u", _triangle_3, _not_ghost);
  Array<Real>::const_matrix_iterator eigen_it =
      eigen_grad_u_1.begin(spatial_dimension, spatial_dimension);
  Array<Real>::const_matrix_iterator eigen_end =
      eigen_grad_u_1.end(spatial_dimension, spatial_dimension);
  for (; eigen_it != eigen_end; ++eigen_it) {
    const Matrix<Real> & eigen_grad_u = *eigen_it;
    error += eigen_grad_u.norm<L_2>();
  }

  Material & mat_2 = model.getMaterial(1);
  const Array<Real> & eigen_grad_u_2 =
      mat_2.getArray<Real>("eigen_grad_u", _triangle_3, _not_ghost);
  eigen_it = eigen_grad_u_2.begin(spatial_dimension, spatial_dimension);
  eigen_end = eigen_grad_u_2.end(spatial_dimension, spatial_dimension);
  for (; eigen_it != eigen_end; ++eigen_it) {
    const Matrix<Real> & eigen_grad_u = *eigen_it;
    Matrix<Real> diff = (prestrain - eigen_grad_u);
    error += diff.norm<L_2>();
  }

  MaterialIGFEM & mat_3 = dynamic_cast<MaterialIGFEM &>(model.getMaterial(2));
  const Array<Real> & eigen_grad_u_3 =
      mat_3.getArray<Real>("eigen_grad_u", _igfem_triangle_5, _not_ghost);
  UInt * sub_mat_ptr =
      mat_3.getArray<UInt>("sub_material", _igfem_triangle_5, _not_ghost)
          .storage();
  eigen_it = eigen_grad_u_3.begin(spatial_dimension, spatial_dimension);
  eigen_end = eigen_grad_u_3.end(spatial_dimension, spatial_dimension);
  for (; eigen_it != eigen_end; ++eigen_it, ++sub_mat_ptr) {
    const Matrix<Real> & eigen_grad_u = *eigen_it;
    if (!(*sub_mat_ptr)) {
      Matrix<Real> diff = (prestrain - eigen_grad_u);
      error += diff.norm<L_2>();
    } else
      error += eigen_grad_u.norm<L_2>();
  }

  std::cout << "The error in the prestrain is: " << error << std::endl;
  if (std::abs(error) > Math::getTolerance()) {
    std::cout << "The test failed!!!" << std::endl;
    finalize();
    return EXIT_FAILURE;
  }

  /// dump
  model.dump("igfem elements");
  model.dump();

  finalize();
  return EXIT_SUCCESS;
}

/**
 * @file   test_interface_position.cc
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 *
 *
 * @brief  patch test for interface close to standard nodes
 *
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

/* -------------------------------------------------------------------------- */
#include "material_elastic.hh"
#include "solid_mechanics_model_igfem.hh"
/* -------------------------------------------------------------------------- */
using namespace akantu;

int main(int argc, char * argv[]) {

  initialize("material.dat", argc, argv);
  StaticCommunicator & comm = StaticCommunicator::getStaticCommunicator();
  Int psize = comm.getNbProc();
  Int prank = comm.whoAmI();

  /// create a mesh and read the regular elements from the mesh file
  /// mesh creation
  const UInt spatial_dimension = 2;
  Mesh mesh(spatial_dimension);
  akantu::MeshPartition * partition = NULL;
  if (prank == 0) {
    mesh.read("test_igfem_element_orientation.msh");
    partition = new MeshPartitionScotch(mesh, spatial_dimension);
    partition->partitionate(psize);
  }

  /// model creation
  SolidMechanicsModelIGFEM model(mesh);
  model.initParallel(partition);
  delete partition;
  model.initFull();

  /// add fields that should be dumped
  model.setBaseName("regular_elements");
  model.setBaseNameToDumper("igfem elements", "igfem elements");
  model.addDumpField("material_index");
  model.addDumpField("partitions");
  model.addDumpFieldToDumper("igfem elements", "lambda");
  model.addDumpFieldToDumper("igfem elements", "material_index");
  model.addDumpFieldToDumper("igfem elements", "partitions");
  /// dump mesh before the IGFEM interface is created
  model.dump();
  model.dump("igfem elements");

  /// create the interace:
  std::list<SK::Sphere_3> sphere_list;
  Real radius = std::sqrt(0.2 * 0.2 + 0.2 * 0.2);
  Real tol = std::sqrt(2.) * 1e-11;
  radius += tol;
  SK::Sphere_3 sphere_1(SK::Point_3(-0.2, 0.3, 0.), radius * radius);
  sphere_list.push_back(sphere_1);
  model.registerGeometryObject(sphere_list, "gel");
  model.update();

  /// check that the igfem elements have been created correctly
  UInt nb_igfem_triangle_4 = mesh.getNbElement(_igfem_triangle_4, _not_ghost);
  UInt nb_igfem_triangle_5 = mesh.getNbElement(_igfem_triangle_5, _not_ghost);
  comm.allReduce(&nb_igfem_triangle_4, 1, _so_sum);
  comm.allReduce(&nb_igfem_triangle_5, 1, _so_sum);
  if (prank == 0) {
    if ((nb_igfem_triangle_4 != 2) || (nb_igfem_triangle_5 != 1)) {
      std::cout << "something went wrong in the interface creation"
                << std::endl;
      finalize();
      return EXIT_FAILURE;
    }
  }

  /// check that the igfem nodes have been created correctly
  UInt nb_global_nodes = 0;
  for (UInt n = 0; n < mesh.getNbNodes(); ++n) {
    if (mesh.isLocalOrMasterNode(n))
      nb_global_nodes += 1;
  }
  comm.allReduce(&nb_global_nodes, 1, _so_sum);
  if (prank == 0 && nb_global_nodes != 27) {
    std::cout << "something went wrong in the interface creation" << std::endl;
    finalize();
    return EXIT_FAILURE;
  }

  /// in this test the sphere interface is positioned in a way that it
  /// genereates two IGFEM elements (of type _igfem_triangle_4), that are
  /// enriched for
  /// compatibilty with their neighboring elements but they don't
  /// contain a material interface. Therefore, the same material has
  /// to be assigned to both sub-elements.

  /// this is to check that both sub-elements have the same material
  /// properties Note: This check works because all _igfem_triangle_4s
  /// in the mesh do not contain a material interface

  MaterialElastic<spatial_dimension> & aggregate_mat =
      dynamic_cast<MaterialElastic<spatial_dimension> &>(
          model.getMaterial("aggregate"));
  Real lambda_1 = aggregate_mat.getLambda();
  Real mu_1 = aggregate_mat.getMu();

  ElementType type = _igfem_triangle_4;
  GhostType ghost_type = _not_ghost;
  UInt nb_element = mesh.getNbElement(type, ghost_type);
  UInt nb_quads =
      model.getFEEngine("IGFEMFEEngine").getNbIntegrationPoints(type);
  Real * lambda = model.getMaterial("igfem_elastic")
                      .getArray<Real>("lambda", type, ghost_type)
                      .storage();
  Real * mu = model.getMaterial("igfem_elastic")
                  .getArray<Real>("mu", type, ghost_type)
                  .storage();

  Real error = 0;
  for (UInt e = 0; e < nb_element; ++e) {
    for (UInt q = 0; q < nb_quads; ++q, ++lambda, ++mu) {
      error += std::abs(lambda_1 - *lambda);
      error += std::abs(mu_1 - *mu);
    }
  }

  comm.allReduce(&error, 1, _so_sum);
  if (prank == 0) {
    std::cout << "The error is: " << error << std::endl;
    if (error > 1e-14) {
      std::cout << "The test failed!!!!" << std::endl;
      finalize();
      return EXIT_FAILURE;
    }
  }

  /// dump mesh after the IGFEM interface is created
  model.dump();
  model.dump("igfem elements");

  finalize();
  return EXIT_SUCCESS;
}

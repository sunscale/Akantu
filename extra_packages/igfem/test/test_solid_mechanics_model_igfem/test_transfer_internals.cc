/**
 * @file   test_transfer_internals.cc
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 *
 *
 * @brief  test to test the transfer of internals such as damage
 *
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "dumper_paraview.hh"
#include "mesh_geom_common.hh"
#include "solid_mechanics_model_igfem.hh"
/* -------------------------------------------------------------------------- */
using namespace akantu;

bool checkResults(Real & error, UInt & counter, SolidMechanicsModel & model,
                  Real diagonal);

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
  model.setBaseNameToDumper("igfem elements", "igfem elements");
  model.addDumpFieldToDumper("igfem elements", "material_index");
  model.addDumpFieldToDumper("igfem elements", "Sc");
  model.addDumpFieldToDumper("igfem elements", "damage");
  model.addDumpFieldToDumper("igfem elements", "lambda");
  model.addDumpFieldToDumper("igfem elements", "partitions");

  /// set damage state as a function of the position of the quadrature
  /// point the damage state is the absolute value of bary center
  /// position normalized by the half of the diagonal of the mesh (which is a
  /// square)

  /// compute the mesh diagonal
  mesh.computeBoundingBox();
  const Vector<Real> & lower_bounds = mesh.getLowerBounds();
  const Vector<Real> & upper_bounds = mesh.getUpperBounds();
  Real diagonal = upper_bounds.distance(lower_bounds);

  /// compute barycenters and set damage state
  GhostType ghost_type = _not_ghost;
  ElementType el_type = _triangle_3;
  UInt nb_element = mesh.getNbElement(el_type, ghost_type);
  Array<Real> barycenter(nb_element, spatial_dimension);

  Array<Real>::iterator<Vector<Real>> bary_it =
      barycenter.begin(spatial_dimension);
  for (UInt elem = 0; elem < nb_element; ++elem) {
    mesh.getBarycenter(elem, el_type, bary_it->storage(), ghost_type);
    UInt mat_index =
        model.getMaterialByElement(el_type, ghost_type).begin()[elem];
    UInt local_index =
        model.getMaterialLocalNumbering(el_type, ghost_type).begin()[elem];
    Material & mat = model.getMaterial(mat_index);
    Array<Real> & damage = mat.getArray<Real>("damage", el_type, ghost_type);
    Array<Real>::scalar_iterator damage_it = damage.begin();
    damage_it[local_index] = (*bary_it).norm() / (0.5 * diagonal);
    ++bary_it;
  }

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

  /// dump
  model.dump("igfem elements");
  model.dump();

  /// check that internals have been transferred correctly
  Real error = 0.;
  UInt counter = 0;
  bool check_passed = checkResults(error, counter, model, diagonal);
  if (!check_passed) {
    finalize();
    return EXIT_FAILURE;
  }

  comm.allReduce(&error, 1, _so_sum);
  comm.allReduce(&counter, 1, _so_sum);

  if (prank == 0) {
    std::cout << "The error is: " << error << std::endl;
    std::cout << "There are " << counter
              << " igfem quads with damage material in the mesh" << std::endl;
    std::cout
        << "There should be 156 igfem quads with damage material in the mesh"
        << std::endl;
  }
  if (error > 1e-14 || counter != 156) {
    finalize();
    return EXIT_FAILURE;
  }

  // /// grow two of the gel pockets (gel pocket 1 and 3) and repeat the test
  std::list<SK::Sphere_3> new_gel_pocket_list;
  /// grow sphere 2
  SK::Sphere_3 sphere_5(SK::Point_3(0., 0., 0.), 0.15 * 0.15);
  /// grow sphere 3
  SK::Sphere_3 sphere_6(SK::Point_3(-0.75, -0.75, 0.), 0.5 * 0.5);
  new_gel_pocket_list.push_back(sphere_5);
  new_gel_pocket_list.push_back(sphere_2);
  new_gel_pocket_list.push_back(sphere_6);
  new_gel_pocket_list.push_back(sphere_4);

  gel_pocket_list.clear();
  gel_pocket_list = new_gel_pocket_list;
  model.update("mat_1");

  /// check again that internals have been transferred correctly
  error = 0.;
  counter = 0;
  check_passed = checkResults(error, counter, model, diagonal);
  if (!check_passed) {
    finalize();
    return EXIT_FAILURE;
  }

  comm.allReduce(&error, 1, _so_sum);
  comm.allReduce(&counter, 1, _so_sum);

  if (prank == 0) {
    std::cout << "The error is: " << error << std::endl;
    std::cout << "There are " << counter
              << " igfem quads with damage material in the mesh" << std::endl;
    std::cout
        << "There should be 150 igfem quads with damage material in the mesh"
        << std::endl;
  }

  if (error > 1e-14 || counter != 150) {
    finalize();
    return EXIT_FAILURE;
  }

  /// dump
  model.dump("igfem elements");
  model.dump();

  finalize();
  return EXIT_SUCCESS;
}

/* -------------------------------------------------------------------------- */
bool checkResults(Real & error, UInt & counter, SolidMechanicsModel & model,
                  Real diagonal) {
  /// check that damage values have been correctly transferred
  FEEngine & fee = model.getFEEngine("IGFEMFEEngine");
  GhostType ghost_type = _not_ghost;
  Mesh & mesh = model.getMesh();
  UInt spatial_dimension = model.getSpatialDimension();
  bool check_passed = true;

  /// loop over all IGFEM elements of type _not_ghost
  Mesh::type_iterator it =
      mesh.firstType(spatial_dimension, ghost_type, _ek_igfem);
  Mesh::type_iterator last =
      mesh.lastType(spatial_dimension, ghost_type, _ek_igfem);
  for (; it != last; ++it) {
    ElementType igfem_el_type = *it;
    UInt nb_igfem_element = mesh.getNbElement(igfem_el_type);
    UInt nb_quads = fee.getNbIntegrationPoints(igfem_el_type, ghost_type);
    Array<Real> barycenter_igfem(nb_igfem_element, spatial_dimension);
    Array<Real>::vector_iterator bary_it =
        barycenter_igfem.begin(spatial_dimension);
    UInt * conn_val = mesh.getConnectivity(igfem_el_type, ghost_type).storage();
    Array<Real> & nodes = mesh.getNodes();
    UInt nb_parent_nodes = IGFEMHelper::getNbParentNodes(igfem_el_type);
    /// compute the bary center of the underlying parent element
    UInt nb_el_nodes = mesh.getNbNodesPerElement(igfem_el_type);
    for (UInt elem = 0; elem < nb_igfem_element; ++elem) {
      Real local_coord[spatial_dimension * nb_parent_nodes];
      UInt offset = elem * nb_el_nodes;
      for (UInt n = 0; n < nb_parent_nodes; ++n) {
        memcpy(local_coord + n * spatial_dimension,
               nodes.storage() + conn_val[offset + n] * spatial_dimension,
               spatial_dimension * sizeof(Real));
      }
      Math::barycenter(local_coord, nb_parent_nodes, spatial_dimension,
                       bary_it->storage());
      UInt mat_index =
          model.getMaterialByElement(igfem_el_type, ghost_type).begin()[elem];
      Material & mat = model.getMaterial(mat_index);
      UInt local_index =
          model.getMaterialLocalNumbering(igfem_el_type, ghost_type)
              .begin()[elem];
      Array<Real> & damage =
          mat.getArray<Real>("damage", igfem_el_type, ghost_type);
      Array<UInt> & sub_mat =
          mat.getArray<UInt>("sub_material", igfem_el_type, ghost_type);
      Array<Real> & strength =
          mat.getArray<Real>("Sc", igfem_el_type, ghost_type);
      Array<Real>::scalar_iterator damage_it = damage.begin();
      Array<UInt>::scalar_iterator sub_mat_it = sub_mat.begin();
      Array<Real>::scalar_iterator Sc_it = strength.begin();
      for (UInt q = 0; q < nb_quads; ++q) {
        UInt q_global = local_index * nb_quads + q;
        if (sub_mat_it[q_global] == 1) {
          if (std::abs(Sc_it[q_global] - 100) > 1e-15) {
            check_passed = false;
            return check_passed;
          }
          error += std::abs(
              (damage_it[q_global] - (*bary_it).norm() / (0.5 * diagonal)));
          ++counter;
        } else if ((std::abs(Sc_it[q_global]) > Math::getTolerance()) ||
                   (std::abs(damage_it[q_global]) > Math::getTolerance())) {
          check_passed = false;
          return check_passed;
        }
      }
      ++bary_it;
    }
  }
  return check_passed;
}

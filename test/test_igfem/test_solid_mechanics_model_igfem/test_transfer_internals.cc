/**
 * @file   test_transfer_internals.cc
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 *
 *
 * @brief  test to test the transfer of internals such as damage
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

/* -------------------------------------------------------------------------- */
#include "solid_mechanics_model_igfem.hh"
#include "aka_common.hh"
 #include "dumper_paraview.hh"
#include "mesh_geom_common.hh"
/* -------------------------------------------------------------------------- */
using namespace akantu;

class ASRMaterialSelector : public DefaultMaterialIGFEMSelector  {
public:
  ASRMaterialSelector(SolidMechanicsModelIGFEM & model) :
    DefaultMaterialIGFEMSelector(model),
    model(model) {} 

  UInt operator()(const Element & elem) {
    if(Mesh::getKind(elem.type) == _ek_igfem)
      /// choose IGFEM material
      return this->fallback_value_igfem;

    const Mesh & mesh = model.getMesh();
    UInt spatial_dimension = model.getSpatialDimension();
    Vector<Real> barycenter(spatial_dimension);
    mesh.getBarycenter(elem, barycenter);
    if(model.isInside(barycenter))
      return 1;
    return 0;
  }

protected:
  SolidMechanicsModelIGFEM & model;
};


typedef Spherical SK;

int main(int argc, char *argv[]) {

  initialize("material_damage.dat", argc, argv);
  StaticCommunicator & comm = akantu::StaticCommunicator::getStaticCommunicator();
  Int psize = comm.getNbProc();
  Int prank = comm.whoAmI();

  /// problem dimension
  UInt spatial_dimension = 2;  

  /// mesh creation and partioning
  Mesh mesh(spatial_dimension);
  akantu::MeshPartition * partition = NULL;

  if(prank == 0) {
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
  model.addDumpFieldToDumper("igfem elements", "partitions");

  /// set damage state as a function of the position of the quadrature
  /// point the damage state is the absolute value of bary center
  /// position normalized by the half of the diagonal of the mesh (which is a square)
  
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

  Array<Real>::iterator< Vector<Real> > bary_it = barycenter.begin(spatial_dimension);
  for (UInt elem = 0; elem < nb_element; ++elem) {
    mesh.getBarycenter(elem, el_type, bary_it->storage(), ghost_type);
    UInt mat_index = model.getMaterialByElement(el_type, ghost_type).begin()[elem];
    UInt local_index = model.getMaterialLocalNumbering(el_type, ghost_type).begin()[elem];
    Material & mat = model.getMaterial(mat_index);
    Array<Real> & damage = mat.getArray<Real>("damage", el_type, ghost_type); 
    Array<Real>::scalar_iterator damage_it = damage.begin();
    damage_it[local_index] = (*bary_it).norm()/(0.5 * diagonal);
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


  /// check that damage values have been correctly transferred 
  Real error = 0.;
  ElementType igfem_el_type = _igfem_triangle_5;
  UInt nb_igfem_element = mesh.getNbElement(igfem_el_type, ghost_type);
  if (mesh.getNbElement(_igfem_triangle_4)) {
    std::cout << "something went wrong in the interface creation" << std::endl;
    finalize();
    return EXIT_FAILURE;
  }
  UInt counter = 0;
  FEEngine & fee = model.getFEEngine("IGFEMFEEngine");
  UInt nb_quads = fee.getNbIntegrationPoints(igfem_el_type, ghost_type);
  Array<Real> barycenter_igfem(nb_igfem_element, spatial_dimension);
  bary_it = barycenter_igfem.begin(spatial_dimension);
  UInt * conn_val = mesh.getConnectivity(igfem_el_type, ghost_type).storage();
  Array<Real> & nodes = mesh.getNodes();
  UInt nb_parent_nodes = 3; /// compute the bary center of the underlying parent element
  UInt nb_el_nodes = mesh.getNbNodesPerElement(igfem_el_type);
  for (UInt elem = 0; elem < nb_igfem_element; ++elem) {
    Real local_coord[spatial_dimension * nb_parent_nodes];
    UInt offset = elem * nb_el_nodes;
    for (UInt n = 0; n < nb_parent_nodes; ++n) {
      memcpy(local_coord + n * spatial_dimension,
	     nodes.storage() + conn_val[offset + n] * spatial_dimension,
	     spatial_dimension*sizeof(Real));
    }
    Math::barycenter(local_coord, nb_parent_nodes, spatial_dimension, bary_it->storage());
    UInt mat_index = model.getMaterialByElement(igfem_el_type, ghost_type).begin()[elem];
    Material & mat = model.getMaterial(mat_index);
    UInt local_index = model.getMaterialLocalNumbering(igfem_el_type, ghost_type).begin()[elem];
    Array<Real> & damage = mat.getArray<Real>("damage", igfem_el_type, ghost_type); 
    Array<UInt> & sub_mat = mat.getArray<UInt>("sub_material", igfem_el_type, ghost_type);
    Array<Real> & strength = mat.getArray<Real>("Sc", igfem_el_type, ghost_type);
    Array<Real>::scalar_iterator damage_it = damage.begin();
    Array<UInt>::scalar_iterator sub_mat_it = sub_mat.begin();
    Array<Real>::scalar_iterator Sc_it = strength.begin();
    for (UInt q = 0; q < nb_quads; ++q) {
      UInt q_global = local_index * nb_quads + q;
      if (sub_mat_it[q_global] == 1) {
	if (std::abs(Sc_it[q_global] - 100) > 1e-15) {
	  finalize();
	  return EXIT_FAILURE;
	} 
	error += std::abs((damage_it[q_global] -(*bary_it).norm()/(0.5 * diagonal)));
	++counter;
      }
      else if ( (std::abs(Sc_it[q_global]) > Math::getTolerance()) || 
		(std::abs(damage_it[q_global]) > Math::getTolerance()) ) {
	finalize();
	return EXIT_FAILURE;
      }
    } 
    ++bary_it;
  }

  comm.allReduce(&error, 1, _so_sum);
  comm.allReduce(&counter, 1, _so_sum);

  if (prank == 0) {
    std::cout << "The error is: " << error << std::endl;
    std::cout << "There are " << counter << " igfem quads with damage material in the mesh" << std::endl;
    std::cout << "There should be 156 igfem quads with damage material in the mesh" << std::endl;
  }
  if (error > 1e-14 || counter != 156) {
    finalize();
    return EXIT_FAILURE;
  }
  finalize();
  return EXIT_SUCCESS;
}


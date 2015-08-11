/**
 * @file   test_solid_mechanics_model_igfem.cc
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
#include "solid_mechanics_model_igfem.hh"
#include "aka_common.hh"
// #include "mesh_segment_intersector.hh"
// #include "mesh_sphere_intersector.hh"
// #include "geom_helper_functions.hh"
// #include "mesh_geom_common.hh"
/* -------------------------------------------------------------------------- */
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <cmath>
#include <math.h> 
#include "dumper_paraview.hh"
#include "mesh_geom_common.hh"
//#include "dumpable_inline_impl.hh"
/* -------------------------------------------------------------------------- */
using namespace akantu;
void writeMesh(const std::string & filename, const Mesh & mesh);
class Sphere {
public:
  Sphere(const Vector<Real> & center, Real radius, Real tolerance = 0.) : center(center), radius(radius), tolerance(tolerance) {
  }

  bool isInside(const Vector<Real> & point) const {
    return (point.distance(center) < radius + tolerance);
  }

  const Vector<Real> & getCenter() const { return center; }
  Real & getRadius() { return radius; }

protected:
  Vector<Real> center;
  Real radius, tolerance;
};

void growGel(std::list<SK::Sphere_3> & query_list, Real new_radius) {
  std::list<SK::Sphere_3>::const_iterator query_it = query_list.begin(),
    query_end = query_list.end();
  std::list<SK::Sphere_3> sphere_list;
  for (; query_it != query_end ; ++query_it) {
    SK::Sphere_3 sphere(query_it->center(),
			new_radius * new_radius);
    sphere_list.push_back(sphere);
  }
  query_list.clear();
  query_list = sphere_list;
}

void applyBoundaryConditions(SolidMechanicsModelIGFEM & model, Real inner_radius) {
  /// boundary conditions
  Mesh & mesh = model.getMesh();
  mesh.computeBoundingBox();
  const Vector<Real> & lowerBounds = mesh.getLowerBounds();
  const Vector<Real> & upperBounds = mesh.getUpperBounds();
  Real bottom  = lowerBounds(1);
  Real top = upperBounds(1);
  Real left = lowerBounds(0);
  Real right = upperBounds(0);

  Real eps = std::abs((top - bottom) * 1e-12);
  const Array<Real> & pos = mesh.getNodes();
  Array<Real> & disp = model.getDisplacement();
  Array<bool> & boun = model.getBlockedDOFs();
  Real radius = 0;
  Real phi = 0;
  Real mu_2 = 10. / (2. * (1 + 0.3));
  Real mu_3 = 1. / (2. * (1 + 0.25));
  Real outer_radius = 2;

  Real lambda_2  = 2. * 0.3 * mu_2 / (1 - 2 * 0.3);
  Real lambda_3  = 2. * 0.25 * mu_3 / (1 - 2 * 0.25);
  Real alpha = (lambda_3 + mu_3 + mu_2) * outer_radius * outer_radius / ((lambda_2 + mu_2) * inner_radius * inner_radius
									 + (lambda_3 + mu_3) * (outer_radius * outer_radius - inner_radius * inner_radius)+(mu_2 * outer_radius * outer_radius));

  disp.clear();
  boun.clear();
  /// absolute confinement
  for (UInt i = 0; i < mesh.getNbNodes(); ++i) {
 
    if(std::abs(pos(i,0) - left) < eps) {
      radius = std::sqrt(pos(i,0)*pos(i,0) + pos(i,1)*pos(i,1));
      phi = std::atan2(pos(i,1), pos(i,0));
      boun(i,0) = true;
      disp(i,0) = cos(phi) * ( (radius - 4./radius) * alpha + 4./radius );
      boun(i,1) = true;
      disp(i,1) = sin(phi) * ( (radius - 4./radius) * alpha + 4./radius );
    }

    if(std::abs(pos(i,0) - right) < eps) {
      radius = std::sqrt(pos(i,0)*pos(i,0) + pos(i,1)*pos(i,1));
      phi = std::atan2(pos(i,1), pos(i,0));
      boun(i,0) = true;
      disp(i,0) = cos(phi) * ( (radius - 4./radius) * alpha + 4./radius );
      boun(i,1) = true;
      disp(i,1) = sin(phi) * ( (radius - 4./radius) * alpha + 4./radius );
    }

    if(std::abs(pos(i,1) - top) < eps) {
      radius = std::sqrt(pos(i,0)*pos(i,0) + pos(i,1)*pos(i,1));
      phi = std::atan2(pos(i,1), pos(i,0));
      boun(i,0) = true;
      disp(i,0) = cos(phi) * ( (radius - 4./radius) * alpha + 4./radius );
      boun(i,1) = true;
      disp(i,1) = sin(phi) * ( (radius - 4./radius) * alpha + 4./radius );
    }

    if(std::abs(pos(i,1) - bottom) < eps) {
      radius = std::sqrt(pos(i,0)*pos(i,0) + pos(i,1)*pos(i,1));
      phi = std::atan2(pos(i,1), pos(i,0));
      boun(i,0) = true;
      disp(i,0) = cos(phi) * ( (radius - 4./radius) * alpha + 4./radius );
      boun(i,1) = true;
      disp(i,1) = sin(phi) * ( (radius - 4./radius) * alpha + 4./radius );
    }

  }
}

class SphereMaterialSelector : public DefaultMaterialIGFEMSelector  {
public:
  SphereMaterialSelector(Real radius, Vector<Real> & center, SolidMechanicsModelIGFEM & model, Real tolerance = 1.e-12) :
    DefaultMaterialIGFEMSelector(model),
    model(model) {
    spheres.push_back(Sphere(center, radius, tolerance));   
  }

  UInt operator()(const Element & elem) {
    if(Mesh::getKind(elem.type) == _ek_igfem)
      return this->fallback_value_igfem;
    //  return 2;//2model.getMaterialIndex(2);
    const Mesh & mesh = model.getMesh();
    UInt spatial_dimension = model.getSpatialDimension();
    Vector<Real> barycenter(spatial_dimension);
    mesh.getBarycenter(elem, barycenter);
    std::vector<Sphere>::const_iterator iit = spheres.begin();
    std::vector<Sphere>::const_iterator eit = spheres.end();
    for (; iit != eit; ++iit) {
      const Sphere & sp = *iit;
      if(sp.isInside(barycenter)) {
	return 1;//model.getMaterialIndex("inside");;
      }
    }
    return 0;
    //return DefaultMaterialSelector::operator()(elem);
  }

  void update(Real new_radius) {
    std::vector<Sphere>::iterator iit = spheres.begin();
    std::vector<Sphere>::iterator eit = spheres.end();
    for (; iit != eit; ++iit) {
      Real & radius = iit->getRadius();
      radius = new_radius;
    } 
  }
  // UInt operator()(const Element & element) {
  //   if(Mesh::getKind(element.type) == _ek_igfem) 
  //     return this->fallback_value_igfem;
  //   else
  //     return DefaultMaterialSelector::operator()(element);
  // }

protected:
  SolidMechanicsModelIGFEM & model;
  std::vector<Sphere> spheres;
};


typedef Spherical SK;

int main(int argc, char *argv[]) {

  initialize("material.dat", argc, argv);

  /// problem dimension
  UInt spatial_dimension = 2;  

  /// mesh creation
  Mesh mesh(spatial_dimension);
  mesh.read(std::string(argv[1]));
  Math::setTolerance(1e-14);
  /// dumper for the mesh
  // DumperParaview dumper_igfem("mesh_igfem");
  // dumper_igfem.registerMesh(mesh, spatial_dimension, _not_ghost, _ek_igfem);
  // DumperParaview dumper_regular("mesh_regular");
  // dumper_regular.registerMesh(mesh, spatial_dimension, _not_ghost, _ek_regular);
  // dumper_regular.dump();

  /// geometry of inclusion
  Real radius_inclusion = 0.401;
  // Spherical kernel testing the addition of nodes
  SK::Sphere_3 sphere(SK::Point_3(0, 0, 0), radius_inclusion * radius_inclusion);
  std::list<SK::Sphere_3> sphere_list;
  sphere_list.push_back(sphere);
  ID domain_name = "gel";

  /// model creation
  SolidMechanicsModelIGFEM model(mesh);

  SphereMaterialSelector * mat_selector;

  Vector<Real> c_sphere(spatial_dimension);
  c_sphere.clear();

  mat_selector = new SphereMaterialSelector(0.401, c_sphere, model);
  model.setMaterialSelector(*mat_selector);
  model.initFull();

  /// register the sphere list in the model

  /// add fields that should be dumped
  model.setBaseName("regular_elements");
  model.addDumpField("material_index");
  model.addDumpFieldVector("displacement");
  model.addDumpField("blocked_dofs");
  model.addDumpField("stress");
  model.addDumpField("damage");
  model.addDumpField("Sc");

  /// set damage at two points

  // set damage at two points
  // GhostType ghost_type = _not_ghost;
  // ElementType element_type = _triangle_3;
  // Array<Real> & damage = const_cast<Array<Real> &>(model.getMaterial("agg_inside").getInternal("eigenstrain")(element_type, ghost_type));
  // damage(0) = 0.1;
  // damage(1) = 0.2;
  // damage(3) = 0.5;
  // model.dump();
  //  dumper_igfem.getDumper().setMode(iohelper::TEXT);
  /// intersect the mesh with the sphere

  // dumper_igfem.dump();
  model.registerGeometryObject(sphere_list, domain_name);
  model.update(domain_name);
  Real inner_radius = 0.401;
  applyBoundaryConditions(model,inner_radius);


  model.setBaseNameToDumper("igfem elements", "igfem elements");
  model.addDumpFieldToDumper("igfem elements", "lambda");
  model.addDumpFieldVectorToDumper("igfem elements", "real_displacement");
  model.addDumpFieldVectorToDumper("igfem elements", "displacement");
  model.addDumpFieldToDumper("igfem elements", "material_index");
  model.addDumpFieldToDumper("igfem elements", "stress");
  model.addDumpFieldToDumper("igfem elements", "Sc");
  model.addDumpFieldToDumper("igfem elements", "damage");
  model.addDumpFieldToDumper("igfem elements", "eigenstrain");
  model.dump("igfem elements");
  model.dump();


  /// apply eigenstrain
  // GhostType ghost_type = _not_ghost;
  // ElementType element_type = _igfem_triangle_5;
  // Array<Real> & eigenstrain = const_cast<Array<Real> &>(model.getMaterial("igfem_elastic").getInternal("eigenstrain")(element_type, ghost_type));
  // const Array<UInt> & is_inside = model.getMaterial("igfem_elastic").getInternalUInt("sub_material")(element_type, ghost_type);

  // Vector<Real> prestrain(4);
  // prestrain.clear();
  // prestrain(0) = 4.3;
  // prestrain(3) = 4.3;
  // Array<Real>::vector_iterator eig_it = eigenstrain.begin(4);
  // Array<UInt>::const_scalar_iterator sub_mat_it = is_inside.begin();
  // Array<Real>::vector_iterator eig_end = eigenstrain.end(4);
  // for(; eig_it != eig_end; ++eig_it, ++sub_mat_it) {
  //   if((*sub_mat_it) == 0)
  //     *eig_it = prestrain;
  // }
  // bool factorize = false;
  // bool converged = false;
  // Real error; 
  // converged = model.solveStep<_scm_newton_raphson_tangent, _scc_increment>(1e-4, error, 2, factorize);
  model.dump("igfem elements");
  model.dump();
  std::cout << "initial inclusion" << std::endl;

  // Math::setTolerance(1e-14);
  // ///dumper_igfem.dump();
  // Array<UInt> & connectivity_igfem = const_cast<Array<UInt> &>(mesh.getConnectivity(_igfem_triangle_5));


  

  growGel(sphere_list, 0.72);
  mat_selector->update(0.72);
  model.update(domain_name);
  inner_radius = 0.72;
  applyBoundaryConditions(model,inner_radius);
  model.dump("igfem elements");
  model.dump();

  std::cout << "first successful grow" << std::endl;


  growGel(sphere_list, 0.84);
  mat_selector->update(0.84);
  model.update(domain_name);
  inner_radius = 0.84;
  applyBoundaryConditions(model,inner_radius);
  model.dump("igfem elements");
  model.dump();

  std::cout << "second successful grow" << std::endl;


  growGel(sphere_list, 0.96);
  mat_selector->update(0.96);
  model.update(domain_name);
  inner_radius = 0.96;
  applyBoundaryConditions(model,inner_radius);
  model.dump("igfem elements");
  model.dump();

  std::cout << "third successful grow" << std::endl;

  


  finalize();
  return EXIT_SUCCESS;
}

/* -------------------------------------------------------------------------- */
void writeMesh(const std::string & filename, const Mesh & mesh) {
  std::ofstream outfile;
  const Array<Real> & nodes = mesh.getNodes();
  UInt nb_nodes = nodes.getSize();
  outfile.open(filename.c_str());

  outfile << "$MeshFormat" << std::endl;
  outfile << "2.1 0 8" << std::endl;;
  outfile << "$EndMeshFormat" << std::endl;;

  outfile << std::setprecision(std::numeric_limits<Real>::digits10);
  outfile << "$Nodes" << std::endl;;
  outfile << nodes.getSize() << std::endl;

  outfile << std::uppercase;
  for(UInt i = 0; i < nodes.getSize(); ++i) {
    Int offset = i * nodes.getNbComponent();
    outfile << i+1;
    for(UInt j = 0; j < nodes.getNbComponent(); ++j) {
      outfile << " " << nodes.storage()[offset + j];
    }

    for (UInt p = nodes.getNbComponent(); p < 3; ++p)
      outfile << " " << 0.;
    outfile << std::endl;;
  }

  outfile << std::nouppercase;
  outfile << "$EndNodes" << std::endl;;


  outfile << "$Elements" << std::endl;;

  Mesh::type_iterator it  = mesh.firstType(_all_dimensions, _not_ghost, _ek_regular);
  Mesh::type_iterator end = mesh.lastType(_all_dimensions, _not_ghost, _ek_regular);

  Int nb_elements = 0;
  for(; it != end; ++it) {
    const Array<UInt> & connectivity = mesh.getConnectivity(*it, _not_ghost);
    nb_elements += connectivity.getSize();
  }

  it  = mesh.firstType(_all_dimensions, _not_ghost, _ek_igfem);
  end = mesh.lastType(_all_dimensions, _not_ghost, _ek_igfem);

  for(; it != end; ++it) {
    const Array<UInt> & connectivity = mesh.getConnectivity(*it, _not_ghost);
    nb_elements += (2 * connectivity.getSize());
  }


  outfile << nb_elements << std::endl;
  end = mesh.lastType(_all_dimensions, _not_ghost, _ek_regular);
  UInt element_idx = 1;
  for(it  = mesh.firstType(_all_dimensions, _not_ghost, _ek_regular); it != end; ++it) {
    ElementType type = *it;
    const Array<UInt> & connectivity = mesh.getConnectivity(type, _not_ghost);

    UInt * tag[2] = {NULL, NULL};
    try {
      const Array<UInt> & data_tag_0 = mesh.getData<UInt>("tag_0", type, _not_ghost);
      tag[0] = data_tag_0.storage();
    } catch(...) { tag[0] = NULL; }

    try {
      const Array<UInt> & data_tag_1 = mesh.getData<UInt>("tag_1", type, _not_ghost);
      tag[1] = data_tag_1.storage();
    } catch(...) { tag[1] = NULL; }

    for(UInt i = 0; i < connectivity.getSize(); ++i) {
      UInt offset = i * connectivity.getNbComponent();
      outfile << element_idx << " " << 2 << " 2";

      /// \todo write the real data in the file
      for (UInt t = 0; t < 2; ++t)
        if(tag[t]) outfile << " " << tag[t][i];
        else outfile << " 0";

      for(UInt j = 0; j < connectivity.getNbComponent(); ++j) {
        outfile << " " << connectivity.storage()[offset + j] + 1;
      }
      outfile << std::endl;
      element_idx++;
    }
  }
  end = mesh.lastType(_all_dimensions, _not_ghost, _ek_igfem);
  for(it  = mesh.firstType(_all_dimensions, _not_ghost, _ek_igfem); it != end; ++it) {
    ElementType type = *it;
    const Array<UInt> & connectivity = mesh.getConnectivity(type, _not_ghost);

    UInt * tag[2] = {NULL, NULL};
    try {
      const Array<UInt> & data_tag_0 = mesh.getData<UInt>("tag_0", type, _not_ghost);
      tag[0] = data_tag_0.storage();
    } catch(...) { tag[0] = NULL; }

    try {
      const Array<UInt> & data_tag_1 = mesh.getData<UInt>("tag_1", type, _not_ghost);
      tag[1] = data_tag_1.storage();
    } catch(...) { tag[1] = NULL; }

    for(UInt i = 0; i < connectivity.getSize(); ++i) {
      UInt offset = i * connectivity.getNbComponent();

      /// write first sub-element
      outfile << element_idx << " " << 2 << " 2";

      /// \todo write the real data in the file
      for (UInt t = 0; t < 2; ++t)
        if(tag[t]) outfile << " " << tag[t][i];
        else outfile << " 0";

      
      outfile << " " << connectivity.storage()[offset + 0] + 1;
      outfile << " " << connectivity.storage()[offset + 3] + 1;
      outfile << " " << connectivity.storage()[offset + 4] + 1;

      outfile << std::endl;
      element_idx++;
      // /// write second sub-element
      // outfile << element_idx << " " << 3 << " 2";

      // /// \todo write the real data in the file
      // for (UInt t = 0; t < 2; ++t)
      //   if(tag[t]) outfile << " " << tag[t][i];
      //   else outfile << " 0";

      
      // outfile << " " << connectivity.storage()[offset + 3] + 1;
      // outfile << " " << connectivity.storage()[offset + 1] + 1;
      // outfile << " " << connectivity.storage()[offset + 2] + 1;
      // outfile << " " << connectivity.storage()[offset + 4] + 1;
	
      // outfile << std::endl;
      // element_idx++;
    }
  

    for(UInt i = 0; i < connectivity.getSize(); ++i) {
      UInt offset = i * connectivity.getNbComponent();

      /// write second sub-element
      outfile << element_idx << " " << 3 << " 2";

      /// \todo write the real data in the file
      for (UInt t = 0; t < 2; ++t)
        if(tag[t]) outfile << " " << tag[t][i];
        else outfile << " 0";

      
      outfile << " " << connectivity.storage()[offset + 3] + 1;
      outfile << " " << connectivity.storage()[offset + 1] + 1;
      outfile << " " << connectivity.storage()[offset + 2] + 1;
      outfile << " " << connectivity.storage()[offset + 4] + 1;
	
      outfile << std::endl;
      element_idx++;
    }
  }

  outfile << "$EndElements" << std::endl;;

  outfile.close();
}


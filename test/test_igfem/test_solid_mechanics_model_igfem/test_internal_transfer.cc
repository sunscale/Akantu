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
  model.dump();

  /// set damage at two points
  GhostType ghost_type = _not_ghost;
  ElementType element_type = _triangle_3;
  Array<Real> & damage = const_cast<Array<Real> &>(model.getMaterial("agg_inside").getInternal<Real>("damage")(element_type, ghost_type));
  damage(0) = 0.1;
  damage(1) = 0.2;
  damage(3) = 0.5;
  model.dump();

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

  finalize();
  return EXIT_SUCCESS;
}


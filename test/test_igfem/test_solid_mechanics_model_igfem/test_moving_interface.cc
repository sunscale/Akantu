/**
 * @file   test_solid_mechanics_model_igfem.cc
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 *
 *
 * @brief  test the solidmechanics model for a moving interface
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
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <cmath>
#include <math.h> 
#include "dumper_paraview.hh"
#include "mesh_geom_common.hh"

/* -------------------------------------------------------------------------- */
using namespace akantu;

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

Real computeAlpha(Real inner_radius, Real outer_radius, const Vector<Real> & lambda, const Vector<Real> & mu) {

  Real alpha = (lambda(1) + mu(1) + mu(0)) * outer_radius * outer_radius / ((lambda(0) + mu(0)) * inner_radius * inner_radius
									 + (lambda(1) + mu(1)) * (outer_radius * outer_radius - inner_radius * inner_radius)+(mu(0) * outer_radius * outer_radius));

  return alpha;
}

void applyBoundaryConditions(SolidMechanicsModelIGFEM & model, Real inner_radius, Real outer_radius, const Vector<Real> & lambda, const Vector<Real> & mu) {
  /// boundary conditions for circular inclusion:
  Real alpha = computeAlpha(inner_radius, outer_radius, lambda, mu);

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
  SphereMaterialSelector(std::vector<Sphere> & sphere_list, SolidMechanicsModelIGFEM & model) :
    DefaultMaterialIGFEMSelector(model),
    model(model),
    spheres(sphere_list) { 
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


protected:
  SolidMechanicsModelIGFEM & model;
  std::vector<Sphere> spheres;
};


typedef Spherical SK;


/// the following modeling problem is explained in:
/// T.-P. Fries "A corrected XFEM approximation without problems in blending elements", 2008 
int main(int argc, char *argv[]) {

  initialize("material.dat", argc, argv);

  /// problem dimension
  const UInt spatial_dimension = 2;  

  /// mesh creation
  Mesh mesh(spatial_dimension);
  mesh.read("plate.msh");
  Math::setTolerance(1e-14);

  /// geometry of IGFEM interface: circular inclusion
  Real radius_inclusion = 0.401;
  Vector<Real> center(spatial_dimension, 0.);
  /// @todo: Simplify this: need to create two type of spheres:
  /// one for the geometry and one for the material selector
  SK::Sphere_3 sphere(SK::Point_3(center(0), center(1), 0), radius_inclusion * radius_inclusion);
  std::list<SK::Sphere_3> sphere_list;
  sphere_list.push_back(sphere);
  ID domain_name = "gel";

  /// model creation
  SolidMechanicsModelIGFEM model(mesh);
  SphereMaterialSelector * mat_selector;
  
  /// set material selector and initialize the model completely
  std::vector<Sphere> spheres;
  spheres.push_back(Sphere(center, radius_inclusion, 1.e-12));
  mat_selector = new SphereMaterialSelector(spheres, model);
  model.setMaterialSelector(*mat_selector);
  model.initFull();

  /// register the sphere list in the model
  model.registerGeometryObject(sphere_list, domain_name);

  /// add fields that should be dumped
  model.setBaseName("regular_elements");
  model.setBaseNameToDumper("igfem elements", "igfem elements");
  model.addDumpField("material_index");
  model.addDumpFieldVector("displacement");
  model.addDumpField("blocked_dofs");
  model.addDumpField("stress");
  model.addDumpFieldToDumper("igfem elements", "lambda");
  model.addDumpFieldVectorToDumper("igfem elements", "real_displacement");
  model.addDumpFieldVectorToDumper("igfem elements", "displacement");
  model.addDumpFieldToDumper("igfem elements", "material_index");
  model.addDumpFieldToDumper("igfem elements", "stress");


  /// dump mesh before the IGFEM interface is created
  model.dump();
  model.dump("igfem elements");

  /// create the interface
  model.update(domain_name);
  /// dump the mesh after the IGFEM interface has been created
  model.dump();
  model.dump("igfem elements");


  /// move interface
  growGel(sphere_list, 1.0);
  mat_selector->update(1.0);
  model.update(domain_name);

  /// dump after the interface has moved
  model.dump();
  model.dump("igfem elements");

  finalize();
  return EXIT_SUCCESS;
}



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
#include "material_damage_iterative.hh"
#include "material_igfem_saw_tooth_damage.hh"
//#include "dumpable_inline_impl.hh"
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

void applyBoundaryConditions(SolidMechanicsModelIGFEM & model) {
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

  disp.clear();
  boun.clear();
  /// free expansion
  for (UInt i = 0; i < mesh.getNbNodes(); ++i) {
 
    if(std::abs(pos(i,1) - bottom) < eps){
      boun(i,1) = true;
      disp(i,1) = 0.0;
    }

    if(std::abs(pos(i,0) - left) < eps){
      boun(i,0) = true;
      disp(i,0) = 0.0;
    }

  }
}

void applyEigenstrain(SolidMechanicsModelIGFEM & model) {
  ///  apply eigenstrain
  /// igfem elements
  GhostType ghost_type = _not_ghost;
  ElementType element_type = _igfem_triangle_5;
  Array<Real> & eigenstrain = const_cast<Array<Real> &>(model.getMaterial("igfem_damage").getInternal<Real>("eigen_grad_u")(element_type, ghost_type));
  const Array<UInt> & is_inside = model.getMaterial("igfem_damage").getInternal<UInt>("sub_material")(element_type, ghost_type);

  Vector<Real> prestrain(4);
  prestrain.clear();
  prestrain(0) = 0.07;
  prestrain(3) = 0.07;
  Array<Real>::vector_iterator eig_it = eigenstrain.begin(4);
  Array<UInt>::const_scalar_iterator sub_mat_it = is_inside.begin();
  Array<Real>::vector_iterator eig_end = eigenstrain.end(4);
  for(; eig_it != eig_end; ++eig_it, ++sub_mat_it) {
    if((*sub_mat_it) == 0)
      *eig_it = prestrain;
  }


  /// regular elements
  element_type = _triangle_3;
  Array<Real> & eigenstrain_reg = const_cast<Array<Real> &>(model.getMaterial("gel").getInternal<Real>("eigen_grad_u")(element_type, ghost_type));

  eig_it = eigenstrain_reg.begin(4);
  eig_end = eigenstrain_reg.end(4);
  for(; eig_it != eig_end; ++eig_it, ++sub_mat_it) {
    *eig_it = prestrain;
  }

}

class SphereMaterialSelector : public DefaultMaterialIGFEMSelector  {
public:
  SphereMaterialSelector(Real radius, Vector<Real> & center, SolidMechanicsModelIGFEM & model, Real aggregate_radius, Real tolerance = 1.e-12) :
    DefaultMaterialIGFEMSelector(model),
    model(model),
    aggregate(center, aggregate_radius, tolerance) {
    gel_pockets.push_back(Sphere(center, radius, tolerance));   
  }

  UInt operator()(const Element & elem) {
    if(Mesh::getKind(elem.type) == _ek_igfem)
      return 3;
    //  return 2;//2model.getMaterialIndex(2);
    const Mesh & mesh = model.getMesh();
    UInt spatial_dimension = model.getSpatialDimension();
    Vector<Real> barycenter(spatial_dimension);
    mesh.getBarycenter(elem, barycenter);
    std::vector<Sphere>::const_iterator iit = gel_pockets.begin();
    std::vector<Sphere>::const_iterator eit = gel_pockets.end();
    for (; iit != eit; ++iit) {
      const Sphere & sp = *iit;
      if(sp.isInside(barycenter)) {
	return 2; //elastic gel material
	//return 0; /// aggregate
      }
    }
    if (aggregate.isInside(barycenter))
      return 0;
    ///return 2;
    return 1;
  }

  void update(Real new_radius) {
    std::vector<Sphere>::iterator iit = gel_pockets.begin();
    std::vector<Sphere>::iterator eit = gel_pockets.end();
    for (; iit != eit; ++iit) {
      Real & radius = iit->getRadius();
      radius = new_radius;
    } 
  }

protected:
  SolidMechanicsModelIGFEM & model;
  std::vector<Sphere> gel_pockets;
  Sphere aggregate;
};


typedef Spherical SK;

int main(int argc, char *argv[]) {

  initialize("material_ASR.dat", argc, argv);

  /// problem dimension
  const UInt spatial_dimension = 2;  

  /// mesh creation
  Mesh mesh(spatial_dimension);
  mesh.read(std::string(argv[1]));
  Math::setTolerance(1e-14);

  Real aggregate_radius = 0.016;
  /// geometry of inclusion
  Real radius_inclusion = 0.0002;///0.0158;//0.0002;
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

  mat_selector = new SphereMaterialSelector(0., c_sphere, model, aggregate_radius);
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
  model.addDumpField("eigen_grad_u");
  // model.dump();

  model.registerGeometryObject(sphere_list, domain_name);
  mat_selector->update(radius_inclusion);
  model.update(domain_name);
  applyBoundaryConditions(model);
  applyEigenstrain(model);
  model.setBaseNameToDumper("igfem elements", "igfem elements");
  model.addDumpFieldToDumper("igfem elements", "lambda");
  model.addDumpFieldVectorToDumper("igfem elements", "real_displacement");
  model.addDumpFieldVectorToDumper("igfem elements", "displacement");
  model.addDumpFieldToDumper("igfem elements", "material_index");
  model.addDumpFieldToDumper("igfem elements", "stress");
  model.addDumpFieldToDumper("igfem elements", "Sc");
  model.addDumpFieldToDumper("igfem elements", "damage");
  model.addDumpFieldToDumper("igfem elements", "eigen_grad_u");


  /// Instantiate objects of class MyDamage
  MaterialDamageIterative<spatial_dimension> & mat_paste = dynamic_cast<MaterialDamageIterative<spatial_dimension> & >(model.getMaterial(1));
  MaterialDamageIterative<spatial_dimension> & mat_aggregate = dynamic_cast<MaterialDamageIterative<spatial_dimension> & >(model.getMaterial(0)); 
  MaterialIGFEMSawToothDamage<spatial_dimension> & mat_igfem = dynamic_cast<MaterialIGFEMSawToothDamage<spatial_dimension> & >(model.getMaterial(3)); 


  bool factorize = false;
  bool converged = false;
  Real error; 

  UInt nb_damaged_elements = 0;
  Real max_eq_stress_aggregate = 0;
  Real max_eq_stress_paste = 0;
  Real max_igfem = 0;

  do {
    converged = model.solveStep<_scm_newton_raphson_tangent, _scc_increment>(1e-4, error, 2, factorize);
    /// compute damage 
    max_eq_stress_aggregate = mat_aggregate.getNormMaxEquivalentStress();
    max_eq_stress_paste = mat_paste.getNormMaxEquivalentStress();
    max_igfem = mat_igfem.getNormMaxEquivalentStress();
    if (max_eq_stress_aggregate > max_eq_stress_paste)
      if (max_eq_stress_aggregate > max_igfem)
	nb_damaged_elements = mat_aggregate.updateDamage();
      else
	nb_damaged_elements = mat_igfem.updateDamage();
    else
      if (max_eq_stress_paste > max_igfem)
	nb_damaged_elements = mat_paste.updateDamage();
      else
	nb_damaged_elements = mat_igfem.updateDamage();
  } while (nb_damaged_elements);

  model.dump("igfem elements");
  model.dump();




  /// grow the gel i = 12
  for (UInt i = 1; i < 12; ++i) {
    // Real new_radius = radius_inclusion - i * 0.00015;
    Real new_radius = radius_inclusion + i * 0.0002;
    growGel(sphere_list, new_radius);
    mat_selector->update(new_radius);
    model.update(domain_name);

    applyBoundaryConditions(model);
    applyEigenstrain(model);
    factorize = false;
    converged = false; 
    nb_damaged_elements = 0;
    do {
      converged = model.solveStep<_scm_newton_raphson_tangent, _scc_increment>(1e-4, error, 2, factorize);
      /// compute damage 
      max_eq_stress_aggregate = mat_aggregate.getNormMaxEquivalentStress();
      max_eq_stress_paste = mat_paste.getNormMaxEquivalentStress();
      max_igfem = mat_igfem.getNormMaxEquivalentStress();
      if (max_eq_stress_aggregate > max_eq_stress_paste)
	if (max_eq_stress_aggregate > max_igfem)
	  nb_damaged_elements = mat_aggregate.updateDamage();
	else
	  nb_damaged_elements = mat_igfem.updateDamage();
      else
	if (max_eq_stress_paste > max_igfem)
	  nb_damaged_elements = mat_paste.updateDamage();
	else
	  nb_damaged_elements = mat_igfem.updateDamage();
  
     // if (max_eq_stress_aggregate > max_eq_stress_paste)
      // 	nb_damaged_elements = mat_aggregate.updateDamage();
      // else
      // 	nb_damaged_elements = mat_paste.updateDamage();
    } while (nb_damaged_elements);

 
    model.dump();
    model.dump("igfem elements");

    std::cout << "the step is " << i <<std::endl;
  }





  model.dump("igfem elements");
  model.dump();


  // // Math::setTolerance(1e-14);
  // // ///dumper_igfem.dump();
  // // Array<UInt> & connectivity_igfem = const_cast<Array<UInt> &>(mesh.getConnectivity(_igfem_triangle_5));

  
  // growGel(sphere_list, 0.72);
  // mat_selector->update(0.72);
  // model.update(domain_name);
  // inner_radius = 0.72;
  // applyBoundaryConditions(model,inner_radius);
  // factorize = false;
  // converged = false;
  // model.dump("igfem elements");
  // model.dump();
  // converged = model.solveStep<_scm_newton_raphson_tangent, _scc_increment>(1e-4, error, 2, factorize);
  // model.dump("igfem elements");
  // model.dump();
  

  finalize();
  return EXIT_SUCCESS;
}


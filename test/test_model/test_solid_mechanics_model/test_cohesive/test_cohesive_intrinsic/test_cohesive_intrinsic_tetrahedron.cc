/**
 * @file   test_cohesive_intrinsic_tetrahedron.cc
 *
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Tue Aug 27 2013
 * @date last modification: Fri Sep 19 2014
 *
 * @brief  Test for cohesive elements
 *
 * @section LICENSE
 *
 * Copyright (©) 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include <limits>
#include <fstream>
#include <iostream>

/* -------------------------------------------------------------------------- */
#include "solid_mechanics_model_cohesive.hh"
#include "material_cohesive.hh"
/* -------------------------------------------------------------------------- */

using namespace akantu;

void updateDisplacement(SolidMechanicsModelCohesive &,
			Array<UInt> &,
			ElementType,
			Vector<Real> &);

bool checkTractions(SolidMechanicsModelCohesive & model,
		    ElementType type,
		    Vector<Real> & opening,
		    Vector<Real> & theoretical_traction,
		    Matrix<Real> & rotation);

void findNodesToCheck(const Mesh & mesh,
		      const Array<UInt> & elements,
		      ElementType type,
		      Array<UInt> & nodes_to_check);

bool checkEquilibrium(const Array<Real> & residual);

bool checkResidual(const Array<Real> & residual,
		   const Vector<Real> & traction,
		   const Array<UInt> & nodes_to_check,
		   const Matrix<Real> & rotation);

int main(int argc, char *argv[]) {
  initialize("material_tetrahedron.dat", argc, argv);

  //  debug::setDebugLevel(dblDump);

  const UInt spatial_dimension = 3;
  const UInt max_steps = 60;
  const Real increment_constant = 0.01;
  Math::setTolerance(1.e-12);

  const ElementType type = _tetrahedron_10;

  Mesh mesh(spatial_dimension);
  mesh.read("tetrahedron.msh");

  SolidMechanicsModelCohesive model(mesh);

  /// model initialization
  model.initFull();

  model.limitInsertion(_x, -0.01, 0.01);
  model.insertIntrinsicElements();

  Array<bool> & boundary = model.getBlockedDOFs();
  boundary.set(true);

  UInt nb_element = mesh.getNbElement(type);

  model.updateResidual();

  model.setBaseName("intrinsic_tetrahedron");
  model.addDumpFieldVector("displacement");
  model.addDumpField("residual");
  model.dump();

  model.setBaseNameToDumper("cohesive elements", "cohesive_elements_tetrahedron");
  model.addDumpFieldVectorToDumper("cohesive elements", "displacement");
  model.addDumpFieldToDumper("cohesive elements", "damage");
  model.dump("cohesive elements");

  /// find elements to displace
  Array<UInt> elements;
  Real * bary = new Real[spatial_dimension];
  for (UInt el = 0; el < nb_element; ++el) {
    mesh.getBarycenter(el, type, bary);
    if (bary[0] > 0.01) elements.push_back(el);
  }
  delete[] bary;

  /// find nodes to check
  Array<UInt> nodes_to_check;
  findNodesToCheck(mesh, elements, type, nodes_to_check);

  /// rotate mesh
  Real angle = 1.;

  Matrix<Real> rotation(spatial_dimension, spatial_dimension);
  rotation.clear();
  rotation(0, 0) = std::cos(angle);
  rotation(0, 1) = std::sin(angle) * -1.;
  rotation(1, 0) = std::sin(angle);
  rotation(1, 1) = std::cos(angle);
  rotation(2, 2) = 1.;


  Vector<Real> increment_tmp(spatial_dimension);
  for (UInt dim = 0; dim < spatial_dimension; ++dim) {
    increment_tmp(dim) = (dim + 1) * increment_constant;
  }

  Vector<Real> increment(spatial_dimension);
  increment.mul<false>(rotation, increment_tmp);

  Array<Real> & position = mesh.getNodes();
  Array<Real> position_tmp(position);

  Array<Real>::iterator<Vector<Real> > position_it = position.begin(spatial_dimension);
  Array<Real>::iterator<Vector<Real> > position_end = position.end(spatial_dimension);
  Array<Real>::iterator<Vector<Real> > position_tmp_it
    = position_tmp.begin(spatial_dimension);

  for (; position_it != position_end; ++position_it, ++position_tmp_it)
    position_it->mul<false>(rotation, *position_tmp_it);

  model.dump();
  model.dump("cohesive elements");

  updateDisplacement(model, elements, type, increment);


  Real theoretical_Ed = 0;

  Vector<Real> opening(spatial_dimension);
  Vector<Real> traction(spatial_dimension);
  Vector<Real> opening_old(spatial_dimension);
  Vector<Real> traction_old(spatial_dimension);

  opening.clear();
  traction.clear();
  opening_old.clear();
  traction_old.clear();

  Vector<Real> Dt(spatial_dimension);
  Vector<Real> Do(spatial_dimension);

  const Array<Real> & residual = model.getResidual();

  /// Main loop
  for (UInt s = 1; s <= max_steps; ++s) {

    model.updateResidual();

    opening += increment_tmp;
    if (checkTractions(model, type, opening, traction, rotation) ||
	checkEquilibrium(residual) ||
	checkResidual(residual, traction, nodes_to_check, rotation)) {
      finalize();
      return EXIT_FAILURE;
    }

    /// compute energy
    Do = opening;
    Do -= opening_old;

    Dt = traction_old;
    Dt += traction;

    theoretical_Ed += .5 * Do.dot(Dt);

    opening_old = opening;
    traction_old = traction;


    updateDisplacement(model, elements, type, increment);

    if(s % 10 == 0) {
      std::cout << "passing step " << s << "/" << max_steps << std::endl;
      model.dump();
      model.dump("cohesive elements");
    }
  }

  model.dump();
  model.dump("cohesive elements");

  Real Ed = model.getEnergy("dissipated");
  theoretical_Ed *= 4.;

  std::cout << Ed << " " << theoretical_Ed << std::endl;

  if (!Math::are_float_equal(Ed, theoretical_Ed) || std::isnan(Ed)) {
    std::cout << "The dissipated energy is incorrect" << std::endl;
    finalize();
    return EXIT_FAILURE;
  }

  finalize();

  std::cout << "OK: test_cohesive_intrinsic_tetrahedron was passed!" << std::endl;
  return EXIT_SUCCESS;
}

/* -------------------------------------------------------------------------- */

void updateDisplacement(SolidMechanicsModelCohesive & model,
			Array<UInt> & elements,
			ElementType type,
			Vector<Real> & increment) {

  UInt spatial_dimension = model.getSpatialDimension();
  Mesh & mesh = model.getFEEngine().getMesh();
  UInt nb_element = elements.getSize();
  UInt nb_nodes = mesh.getNbNodes();
  UInt nb_nodes_per_element = mesh.getNbNodesPerElement(type);

  const Array<UInt> & connectivity = mesh.getConnectivity(type);
  Array<Real> & displacement = model.getDisplacement();
  Array<bool> update(nb_nodes);
  update.clear();

  for (UInt el = 0; el < nb_element; ++el) {
    for (UInt n = 0; n < nb_nodes_per_element; ++n) {
      UInt node = connectivity(elements(el), n);
      if (!update(node)) {
	Vector<Real> node_disp(displacement.storage() + node * spatial_dimension,
			       spatial_dimension);
	node_disp += increment;
	update(node) = true;
      }
    }
  }
}

/* -------------------------------------------------------------------------- */

bool checkTractions(SolidMechanicsModelCohesive & model,
		    ElementType type,
		    Vector<Real> & opening,
		    Vector<Real> & theoretical_traction,
		    Matrix<Real> & rotation) {
  UInt spatial_dimension = model.getSpatialDimension();

  const MaterialCohesive & mat_cohesive
    = dynamic_cast < const MaterialCohesive & > (model.getMaterial(1));

  Real sigma_c = mat_cohesive.getParam< RandomInternalField<Real, FacetInternalField> >("sigma_c");
  const Real beta = mat_cohesive.getParam<Real>("beta");
  const Real G_cI = mat_cohesive.getParam<Real>("G_cI");
  //  Real G_cII = mat_cohesive.getParam<Real>("G_cII");
  const Real delta_0 = mat_cohesive.getParam<Real>("delta_0");
  const Real kappa = mat_cohesive.getParam<Real>("kappa");
  Real delta_c = 2 * G_cI / sigma_c;
  sigma_c *= delta_c / (delta_c - delta_0);

  ElementType type_facet = Mesh::getFacetType(type);
  ElementType type_cohesive = FEEngine::getCohesiveElementType(type_facet);

  const Array<Real> & traction = mat_cohesive.getTraction(type_cohesive);
  const Array<Real> & damage = mat_cohesive.getDamage(type_cohesive);

  UInt nb_quad_per_el
    = model.getFEEngine("CohesiveFEEngine").getNbQuadraturePoints(type_cohesive);
  UInt nb_element = model.getMesh().getNbElement(type_cohesive);
  UInt tot_nb_quad = nb_element * nb_quad_per_el;

  Vector<Real> normal_opening(spatial_dimension);
  normal_opening.clear();
  normal_opening(0) = opening(0);
  Real normal_opening_norm = normal_opening.norm();

  Vector<Real> tangential_opening(spatial_dimension);
  tangential_opening.clear();
  for (UInt dim = 1; dim < spatial_dimension; ++dim)
    tangential_opening(dim) = opening(dim);

  Real tangential_opening_norm = tangential_opening.norm();

  Real beta2_kappa2 = beta * beta/kappa/kappa;
  Real beta2_kappa  = beta * beta/kappa;

  Real delta = std::sqrt(tangential_opening_norm * tangential_opening_norm
			 * beta2_kappa2 +
			 normal_opening_norm * normal_opening_norm);

  delta = std::max(delta, delta_0);

  Real theoretical_damage = std::min(delta / delta_c, 1.);

  if (Math::are_float_equal(theoretical_damage, 1.))
    theoretical_traction.clear();
  else {
    theoretical_traction  = tangential_opening;
    theoretical_traction *= beta2_kappa;
    theoretical_traction += normal_opening;
    theoretical_traction *= sigma_c / delta * (1. - theoretical_damage);
  }

  // adjust damage
  theoretical_damage = std::max((delta - delta_0) / (delta_c - delta_0), 0.);
  theoretical_damage = std::min(theoretical_damage, 1.);

  Vector<Real> theoretical_traction_rotated(spatial_dimension);
  theoretical_traction_rotated.mul<false>(rotation, theoretical_traction);

  for (UInt q = 0; q < tot_nb_quad; ++q) {
    for (UInt dim = 0; dim < spatial_dimension; ++dim) {
      if (!Math::are_float_equal(theoretical_traction_rotated(dim),
				 traction(q, dim))) {
	std::cout << "Tractions are incorrect" << std::endl;
	return 1;
      }
    }

    if (!Math::are_float_equal(theoretical_damage, damage(q))) {
      std::cout << "Damage is incorrect" << std::endl;
      return 1;
    }
  }

  return 0;
}

/* -------------------------------------------------------------------------- */

void findNodesToCheck(const Mesh & mesh,
		      const Array<UInt> & elements,
		      ElementType type,
		      Array<UInt> & nodes_to_check) {

  const Array<UInt> & connectivity = mesh.getConnectivity(type);
  const Array<Real> & position = mesh.getNodes();

  UInt nb_nodes = position.getSize();
  UInt nb_nodes_per_elem = connectivity.getNbComponent();

  Array<bool> checked_nodes(nb_nodes);
  checked_nodes.clear();

  for (UInt el = 0; el < elements.getSize(); ++el) {

    UInt element = elements(el);
    Vector<UInt> conn_el(connectivity.storage() + nb_nodes_per_elem * element,
			 nb_nodes_per_elem);

    for (UInt n = 0; n < nb_nodes_per_elem; ++n) {
      UInt node = conn_el(n);
      if (Math::are_float_equal(position(node, 0), 0.)
	  && checked_nodes(node) == false) {
	checked_nodes(node) = true;
	nodes_to_check.push_back(node);
      }
    }
  }
}

/* -------------------------------------------------------------------------- */

bool checkEquilibrium(const Array<Real> & residual) {

  UInt spatial_dimension = residual.getNbComponent();

  Vector<Real> residual_sum(spatial_dimension);
  residual_sum.clear();

  Array<Real>::const_iterator<Vector<Real> > res_it
    = residual.begin(spatial_dimension);
  Array<Real>::const_iterator<Vector<Real> > res_end
    = residual.end(spatial_dimension);

  for (; res_it != res_end; ++res_it)
    residual_sum += *res_it;

  for (UInt s = 0; s < spatial_dimension; ++s) {
    if (!Math::are_float_equal(residual_sum(s), 0.)) {
      std::cout << "System is not in equilibrium!" << std::endl;
      return 1;
    }
  }

  return 0;
}

/* -------------------------------------------------------------------------- */

bool checkResidual(const Array<Real> & residual,
		   const Vector<Real> & traction,
		   const Array<UInt> & nodes_to_check,
		   const Matrix<Real> & rotation) {

  UInt spatial_dimension = residual.getNbComponent();

  Vector<Real> total_force(spatial_dimension);
  total_force.clear();

  for (UInt n = 0; n < nodes_to_check.getSize(); ++n) {
    UInt node = nodes_to_check(n);

    Vector<Real> res(residual.storage() + node * spatial_dimension,
		     spatial_dimension);

    total_force += res;
  }

  Vector<Real> theoretical_total_force(spatial_dimension);
  theoretical_total_force.mul<false>(rotation, traction);
  theoretical_total_force *= -1 * 2 * 2;

  for (UInt s = 0; s < spatial_dimension; ++s) {
    if (!Math::are_float_equal(total_force(s), theoretical_total_force(s))) {
      std::cout << "Total force isn't correct!" << std::endl;
      return 1;
    }
  }

  return 0;
}

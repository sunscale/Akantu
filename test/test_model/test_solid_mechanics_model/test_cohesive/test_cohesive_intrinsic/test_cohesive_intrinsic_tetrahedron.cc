/**
 * @file   test_cohesive_intrinsic_tetrahedron.cc
 *
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Tue Aug 27 2013
 * @date last modification: Mon Dec 18 2017
 *
 * @brief  Test for cohesive elements
 *
 *
 * Copyright (©) 2014-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#include <fstream>
#include <iostream>
#include <limits>

/* -------------------------------------------------------------------------- */
#include "material_cohesive.hh"
#include "solid_mechanics_model_cohesive.hh"
/* -------------------------------------------------------------------------- */

using namespace akantu;

class Checker {
public:
  Checker(const SolidMechanicsModelCohesive & model,
          const Array<UInt> & elements, ElementType type);

  void check(const Vector<Real> & opening, const Matrix<Real> & rotation) {
    checkTractions(opening, rotation);
    checkEquilibrium();
    computeEnergy(opening);
  }

  void updateDisplacement(const Vector<Real> & increment);

protected:
  void checkTractions(const Vector<Real> & opening,
                      const Matrix<Real> & rotation);
  void checkEquilibrium();
  void checkResidual(const Matrix<Real> & rotation);
  void computeEnergy(const Vector<Real> & opening);

private:
  std::set<UInt> nodes_to_check;
  const SolidMechanicsModelCohesive & model;
  ElementType type;
  // const Array<UInt> & elements;

  const Material & mat_cohesive;

  Real sigma_c;
  const Real beta;
  const Real G_c;
  const Real delta_0;
  const Real kappa;
  Real delta_c;

  const UInt spatial_dimension;

  const ElementType type_facet;
  const ElementType type_cohesive;
  const Array<Real> & traction;
  const Array<Real> & damage;
  const UInt nb_quad_per_el;
  const UInt nb_element;

  const Real beta2_kappa2;
  const Real beta2_kappa;

  Vector<Real> theoretical_traction;
  Vector<Real> traction_old;
  Vector<Real> opening_old;

  Real Ed;
};

/* -------------------------------------------------------------------------- */
int main(int argc, char * argv[]) {
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
  model.getElementInserter().setLimit(_x, -0.01, 0.01);

  /// model initialization
  model.initFull();

  Array<bool> & boundary = model.getBlockedDOFs();
  boundary.set(true);
  UInt nb_element = mesh.getNbElement(type);

  model.setBaseName("intrinsic_tetrahedron");
  model.addDumpFieldVector("displacement");
  model.addDumpField("internal_force");
  model.dump();

  model.setBaseNameToDumper("cohesive elements",
                            "cohesive_elements_tetrahedron");
  model.addDumpFieldVectorToDumper("cohesive elements", "displacement");
  model.addDumpFieldToDumper("cohesive elements", "damage");
  model.dump("cohesive elements");

  /// find elements to displace
  Array<UInt> elements;
  Vector<Real> bary(spatial_dimension);
  for (UInt el = 0; el < nb_element; ++el) {
    mesh.getBarycenter({type, el, _not_ghost}, bary);
    if (bary(_x) > 0.01)
      elements.push_back(el);
  }

  /// rotate mesh
  Real angle = 1.;

  // clang-format off
  Matrix<Real> rotation{
    {std::cos(angle), std::sin(angle) * -1., 0.},
    {std::sin(angle), std::cos(angle),       0.},
    {0.,              0.,                    1.}};
  // clang-format on

  Vector<Real> increment_tmp{increment_constant, 2. * increment_constant,
                             3. * increment_constant};
  Vector<Real> increment = rotation * increment_tmp;

  auto & position = mesh.getNodes();
  auto position_it = position.begin(spatial_dimension);
  auto position_end = position.end(spatial_dimension);

  for (; position_it != position_end; ++position_it) {
    auto & pos = *position_it;
    pos = rotation * pos;
  }

  model.dump();
  model.dump("cohesive elements");

  /// find nodes to check
  Checker checker(model, elements, type);

  checker.updateDisplacement(increment);

  Real theoretical_Ed = 0;

  Vector<Real> opening(spatial_dimension, 0.);
  Vector<Real> opening_old(spatial_dimension, 0.);

  /// Main loop
  for (UInt s = 1; s <= max_steps; ++s) {
    model.solveStep();

    model.dump();
    model.dump("cohesive elements");

    opening += increment_tmp;
    checker.check(opening, rotation);

    checker.updateDisplacement(increment);
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

  std::cout << "OK: test_cohesive_intrinsic_tetrahedron was passed!"
            << std::endl;
  return EXIT_SUCCESS;
}

/* -------------------------------------------------------------------------- */
void Checker::updateDisplacement(const Vector<Real> & increment) {
  Mesh & mesh = model.getFEEngine().getMesh();
  const auto & connectivity = mesh.getConnectivity(type);
  auto & displacement = model.getDisplacement();
  Array<bool> update(displacement.size());
  update.clear();

  auto conn_it = connectivity.begin(connectivity.getNbComponent());
  auto conn_end = connectivity.begin(connectivity.getNbComponent());
  for (; conn_it != conn_end; ++conn_it) {
    const auto & conn = *conn_it;
    for (UInt n = 0; n < conn.size(); ++n) {
      UInt node = conn(n);
      if (!update(node)) {
        Vector<Real> node_disp(displacement.storage() +
                                   node * spatial_dimension,
                               spatial_dimension);
        node_disp += increment;
        update(node) = true;
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
Checker::Checker(const SolidMechanicsModelCohesive & model,
                 const Array<UInt> & elements, ElementType type)
    : model(model), type(std::move(type)), // elements(elements),
      mat_cohesive(model.getMaterial(1)), sigma_c(mat_cohesive.get("sigma_c")),
      beta(mat_cohesive.get("beta")), G_c(mat_cohesive.get("G_c")),
      delta_0(mat_cohesive.get("delta_0")), kappa(mat_cohesive.get("kappa")),
      spatial_dimension(model.getSpatialDimension()),
      type_facet(Mesh::getFacetType(type)),
      type_cohesive(FEEngine::getCohesiveElementType(type_facet)),
      traction(mat_cohesive.getArray<Real>("tractions", type_cohesive)),
      damage(mat_cohesive.getArray<Real>("damage", type_cohesive)),
      nb_quad_per_el(model.getFEEngine("CohesiveFEEngine")
                         .getNbIntegrationPoints(type_cohesive)),
      nb_element(model.getMesh().getNbElement(type_cohesive)),
      beta2_kappa2(beta * beta / kappa / kappa),
      beta2_kappa(beta * beta / kappa) {
  const Mesh & mesh = model.getMesh();
  const auto & connectivity = mesh.getConnectivity(type);
  const auto & position = mesh.getNodes();
  auto conn_it = connectivity.begin(connectivity.getNbComponent());
  for (const auto & element : elements) {
    Vector<UInt> conn_el(conn_it[element]);
    for (UInt n = 0; n < conn_el.size(); ++n) {
      UInt node = conn_el(n);
      if (Math::are_float_equal(position(node, _x), 0.))
        nodes_to_check.insert(node);
    }
  }

  delta_c = 2 * G_c / sigma_c;
  sigma_c *= delta_c / (delta_c - delta_0);
}

/* -------------------------------------------------------------------------- */
void Checker::checkTractions(const Vector<Real> & opening,
                             const Matrix<Real> & rotation) {
  auto normal_opening = opening * Vector<Real>{1., 0., 0.};
  auto tangential_opening = opening - normal_opening;
  const Real normal_opening_norm = normal_opening.norm();
  const Real tangential_opening_norm = tangential_opening.norm();
  const Real delta =
      std::max(std::sqrt(tangential_opening_norm * tangential_opening_norm *
                             beta2_kappa2 +
                         normal_opening_norm * normal_opening_norm),
               0.);

  Real theoretical_damage = std::min(delta / delta_c, 1.);

  theoretical_traction = (tangential_opening * beta2_kappa + normal_opening) *
                         sigma_c / delta * (1. - theoretical_damage);
  // adjust damage
  theoretical_damage = std::max((delta - delta_0) / (delta_c - delta_0), 0.);
  theoretical_damage = std::min(theoretical_damage, 1.);

  Vector<Real> theoretical_traction_rotated = rotation * theoretical_traction;

  std::for_each(
      traction.begin(spatial_dimension), traction.end(spatial_dimension),
      [&theoretical_traction_rotated](auto && traction) {
        Real diff =
            Vector<Real>(theoretical_traction_rotated - traction).norm<L_inf>();
        if (diff > 1e-14)
          throw std::domain_error("Tractions are incorrect");
      });

  std::for_each(damage.begin(), damage.end(),
                [&theoretical_damage](auto && damage) {
                  if (not Math::are_float_equal(theoretical_damage, damage))
                    throw std::domain_error("Damage is incorrect");
                });
}

/* -------------------------------------------------------------------------- */
void Checker::computeEnergy(const Vector<Real> & opening) {
  /// compute energy
  auto Do = opening - opening_old;
  auto Dt = traction_old + theoretical_traction;

  Ed += .5 * Do.dot(Dt);

  opening_old = opening;
  traction_old = theoretical_traction;
}

/* -------------------------------------------------------------------------- */
void Checker::checkEquilibrium() {
  Vector<Real> residual_sum(spatial_dimension, 0.);
  const auto & residual = model.getInternalForce();
  auto res_it = residual.begin(spatial_dimension);
  auto res_end = residual.end(spatial_dimension);
  for (; res_it != res_end; ++res_it)
    residual_sum += *res_it;

  if (!Math::are_float_equal(residual_sum.norm<L_inf>(), 0.))
    throw std::domain_error("System is not in equilibrium!");
}

/* -------------------------------------------------------------------------- */
void Checker::checkResidual(const Matrix<Real> & rotation) {
  Vector<Real> total_force(spatial_dimension, 0.);
  const auto & residual = model.getInternalForce();
  for (auto node : nodes_to_check) {
    Vector<Real> res(residual.begin(spatial_dimension)[node]);
    total_force += res;
  }

  Vector<Real> theoretical_total_force(spatial_dimension);
  theoretical_total_force.mul<false>(rotation, theoretical_traction);
  theoretical_total_force *= -1 * 2 * 2;

  for (UInt s = 0; s < spatial_dimension; ++s) {
    if (!Math::are_float_equal(total_force(s), theoretical_total_force(s))) {
      std::cout << "Total force isn't correct!" << std::endl;
      std::terminate();
    }
  }
}

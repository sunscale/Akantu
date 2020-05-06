/**
 * @file   test_cohesive_extrinsic_fatigue.cc
 *
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Fri Feb 20 2015
 * @date last modification: Tue Feb 20 2018
 *
 * @brief  Test for the linear fatigue cohesive law
 *
 *
 * Copyright (©) 2015-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "material_cohesive_linear_fatigue.hh"
#include "solid_mechanics_model_cohesive.hh"
#include <limits>

/* -------------------------------------------------------------------------- */

using namespace akantu;

// the following class contains an implementation of the 1D linear
// fatigue cohesive law
class MaterialFatigue {
public:
  MaterialFatigue(Real delta_f, Real sigma_c, Real delta_c)
      : delta_f(delta_f), sigma_c(sigma_c), delta_c(delta_c), delta_prec(0),
        traction(sigma_c), delta_max(0),
        stiff_plus(std::numeric_limits<Real>::max()),
        tolerance(Math::getTolerance()){};

  Real computeTraction(Real delta) {
    if (delta - delta_c > -tolerance)
      traction = 0;
    else if (delta_max < tolerance && delta < tolerance)
      traction = sigma_c;
    else {
      Real delta_dot = delta - delta_prec;

      if (delta_dot > -tolerance) {
        stiff_plus *= 1 - delta_dot / delta_f;
        traction += stiff_plus * delta_dot;
        Real max_traction = sigma_c * (1 - delta / delta_c);

        if (traction - max_traction > -tolerance || delta_max < tolerance) {
          traction = max_traction;
          stiff_plus = traction / delta;
        }
      } else {
        Real stiff_minus = traction / delta_prec;
        stiff_plus += (stiff_plus - stiff_minus) * delta_dot / delta_f;
        traction += stiff_minus * delta_dot;
      }
    }

    delta_prec = delta;
    delta_max = std::max(delta, delta_max);
    return traction;
  }

private:
  const Real delta_f;
  const Real sigma_c;
  const Real delta_c;
  Real delta_prec;
  Real traction;
  Real delta_max;
  Real stiff_plus;
  const Real tolerance;
};

void imposeOpening(SolidMechanicsModelCohesive &, Real);
void arange(Array<Real> &, Real, Real, Real);

/* -------------------------------------------------------------------------- */
int main(int argc, char * argv[]) {
  initialize("material_fatigue.dat", argc, argv);

  Math::setTolerance(1e-13);

  const UInt spatial_dimension = 2;
  const ElementType type = _quadrangle_4;

  Mesh mesh(spatial_dimension);
  mesh.read("fatigue.msh");

  // init stuff
  const ElementType type_facet = Mesh::getFacetType(type);
  const ElementType type_cohesive =
      FEEngine::getCohesiveElementType(type_facet);

  SolidMechanicsModelCohesive model(mesh);
  model.initFull(
      SolidMechanicsModelCohesiveOptions(_explicit_lumped_mass, true));

  MaterialCohesiveLinearFatigue<2> & numerical_material =
      dynamic_cast<MaterialCohesiveLinearFatigue<2> &>(
          model.getMaterial("cohesive"));

  Real delta_f = numerical_material.getParam("delta_f");
  Real delta_c = numerical_material.getParam("delta_c");
  Real sigma_c = 1;

  const Array<Real> & traction_array =
      numerical_material.getTraction(type_cohesive);

  MaterialFatigue theoretical_material(delta_f, sigma_c, delta_c);

  // model.setBaseName("fatigue");
  // model.addDumpFieldVector("displacement");
  // model.addDumpField("stress");
  // model.dump();

  // stretch material
  Real strain = 1;
  Array<Real> & displacement = model.getDisplacement();
  const Array<Real> & position = mesh.getNodes();

  for (UInt n = 0; n < mesh.getNbNodes(); ++n)
    displacement(n, 0) = position(n, 0) * strain;

  model.assembleInternalForces();
  // model.dump();

  // insert cohesive elements
  model.checkCohesiveStress();

  // create the displacement sequence
  Real increment = 0.01;

  Array<Real> openings;

  arange(openings, 0, 0.5, increment);
  arange(openings, 0.5, 0.1, increment);
  arange(openings, 0.1, 0.7, increment);
  arange(openings, 0.7, 0.3, increment);
  arange(openings, 0.3, 0.6, increment);
  arange(openings, 0.6, 0.3, increment);
  arange(openings, 0.3, 0.7, increment);
  arange(openings, 0.7, 1.3, increment);

  const Array<UInt> & switches = numerical_material.getSwitches(type_cohesive);

  // std::ofstream edis("fatigue_edis.txt");

  // impose openings
  for (UInt i = 0; i < openings.size(); ++i) {

    // compute numerical traction
    imposeOpening(model, openings(i));
    model.assembleInternalForces();
    // model.dump();
    Real numerical_traction = traction_array(0, 0);

    // compute theoretical traction
    Real theoretical_traction =
        theoretical_material.computeTraction(openings(i));

    // test traction
    if (std::abs(numerical_traction - theoretical_traction) > 1e-13)
      AKANTU_ERROR("The numerical traction "
                   << numerical_traction << " and theoretical traction "
                   << theoretical_traction << " are not coincident");

    // edis << model.getEnergy("dissipated") << std::endl;
  }

  if (switches(0) != 7)
    AKANTU_ERROR("The number of switches is wrong");

  std::cout << "OK: the test_cohesive_extrinsic_fatigue passed." << std::endl;
  return 0;
}

/* -------------------------------------------------------------------------- */

void imposeOpening(SolidMechanicsModelCohesive & model, Real opening) {

  UInt spatial_dimension = model.getSpatialDimension();
  Mesh & mesh = model.getFEEngine().getMesh();
  Array<Real> & position = mesh.getNodes();
  Array<Real> & displacement = model.getDisplacement();
  UInt nb_nodes = mesh.getNbNodes();

  Array<bool> update(nb_nodes);
  update.clear();

  Mesh::type_iterator it = mesh.firstType(spatial_dimension);
  Mesh::type_iterator end = mesh.lastType(spatial_dimension);

  for (; it != end; ++it) {
    ElementType type = *it;
    UInt nb_element = mesh.getNbElement(type);
    UInt nb_nodes_per_element = mesh.getNbNodesPerElement(type);

    const Array<UInt> & connectivity = mesh.getConnectivity(type);
    Vector<Real> barycenter(spatial_dimension);

    for (UInt el = 0; el < nb_element; ++el) {
      mesh.getBarycenter({type, el, _not_ghost}, barycenter);
      if (barycenter(0) > 1) {
        for (UInt n = 0; n < nb_nodes_per_element; ++n) {
          UInt node = connectivity(el, n);
          if (!update(node)) {
            displacement(node, 0) = opening + position(node, 0);
            update(node) = true;
          }
        }
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
void arange(Array<Real> & openings, Real begin, Real end, Real increment) {
  if (begin < end) {
    for (Real opening = begin; opening < end - increment / 2.;
         opening += increment)
      openings.push_back(opening);
  } else {
    for (Real opening = begin; opening > end + increment / 2.;
         opening -= increment)
      openings.push_back(opening);
  }
}

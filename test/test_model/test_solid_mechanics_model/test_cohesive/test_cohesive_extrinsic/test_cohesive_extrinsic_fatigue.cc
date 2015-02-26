/**
 * @file   test_cohesive_intrinsic_fatigue.cc
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 * @date   Fri Feb 20 10:13:14 2015
 *
 * @brief  Test for the linear fatigue cohesive law
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "solid_mechanics_model_cohesive.hh"
#include "material_cohesive.hh"

/* -------------------------------------------------------------------------- */

using namespace akantu;

// the following class contains an implementation of the 1D linear
// fatigue cohesive law
class MaterialFatigue {
public:
  MaterialFatigue(Real delta_f, Real sigma_c, Real delta_c) :
    delta_f(delta_f), sigma_c(sigma_c), delta_c(delta_c),
    tolerance(Math::getTolerance()) {};

  Real computeTraction(Real delta) {
    if (delta - delta_c > -tolerance)
      traction = 0;
    else if (delta_max < tolerance && delta < tolerance)
      traction = sigma_c;
    else {
      Real delta_dot = delta - delta_prec;

      if (delta_dot > -tolerance) {
	Real new_stiffness = stiff_plus * (1 - delta_dot / delta_f);
	Real new_traction = traction + new_stiffness * delta_dot;
	Real max_traction = sigma_c * (1 - delta / delta_c);

	if (new_traction - max_traction > -tolerance || delta_max < tolerance) {
	  stiff_plus = max_traction / delta;
	  traction = stiff_plus * delta;
	} else {
	  stiff_plus = new_stiffness;
	  traction = new_traction;
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
int main(int argc, char *argv[]) {
  initialize("material_fatigue.dat", argc, argv);

  const UInt spatial_dimension = 2;
  const ElementType type = _quadrangle_4;

  Mesh mesh(spatial_dimension);
  mesh.read("fatigue.msh");

  // init stuff
  const ElementType type_facet = Mesh::getFacetType(type);
  const ElementType type_cohesive = FEEngine::getCohesiveElementType(type_facet);

  SolidMechanicsModelCohesive model(mesh);
  model.initFull(SolidMechanicsModelCohesiveOptions(_explicit_lumped_mass, true));

  MaterialCohesive & numerical_material
    = dynamic_cast<MaterialCohesive &>(model.getMaterial("cohesive"));

  Real delta_f = numerical_material.getParam<Real>("delta_f");
  Real delta_c = numerical_material.getParam<Real>("delta_c");
  Real sigma_c = 1;

  const Array<Real> & traction_array = numerical_material.getTraction(type_cohesive);

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

  model.updateResidual();
  //  model.dump();

  // insert cohesive elements
  model.checkCohesiveStress();

  // create the displacement sequence
  Real increment = 0.01;

  Array<Real> openings;

  arange(openings,   0, 0.5, increment);
  arange(openings, 0.5, 0.1, increment);
  arange(openings, 0.1, 0.7, increment);
  arange(openings, 0.7, 0.3, increment);
  arange(openings, 0.3, 0.6, increment);
  arange(openings, 0.6, 0.3, increment);
  arange(openings, 0.3, 0.7, increment);
  arange(openings, 0.7, 1.3, increment);

  // impose openings
  for (UInt i = 0; i < openings.getSize(); ++i) {

    // compute numerical traction
    imposeOpening(model, openings(i));
    model.updateResidual();
    Real numerical_traction = traction_array(0, 0);

    // compute theoretical traction
    Real theoretical_traction = theoretical_material.computeTraction(openings(i));

    // test traction
    if (std::abs(numerical_traction - theoretical_traction) > 1e-13)
      AKANTU_DEBUG_ERROR("The numerical traction " << numerical_traction
			 << " and theoretical traction " << theoretical_traction
			 << " are not coincident");
  }

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
      mesh.getBarycenter(el, type, barycenter.storage());
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
    for (Real opening = begin; opening < end - 1e-13; opening += increment)
      openings.push_back(opening);
  } else {
    for (Real opening = begin; opening > end + 1e-13; opening -= increment)
      openings.push_back(opening);
  }
}

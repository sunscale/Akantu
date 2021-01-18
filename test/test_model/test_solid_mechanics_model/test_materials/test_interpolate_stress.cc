/**
 * @file   test_interpolate_stress.cc
 *
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Thu Jun 07 2012
 * @date last modification: Tue Nov 07 2017
 *
 * @brief  Test for the stress interpolation function
 *
 *
 * Copyright (©)  2010-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "integrator_gauss.hh"
#include "mesh_utils.hh"
#include "shape_lagrange.hh"
#include "solid_mechanics_model.hh"
/* -------------------------------------------------------------------------- */

using namespace akantu;

Real function(Real x, Real y, Real z) { return 100. + 2. * x + 3. * y + 4 * z; }

int main(int argc, char * argv[]) {
  initialize("material_interpolate.dat", argc, argv);

  debug::setDebugLevel(dblWarning);

  const UInt spatial_dimension = 3;
  const ElementType type = _tetrahedron_10;

  Mesh mesh(spatial_dimension);
  mesh.read("interpolation.msh");
  const ElementType type_facet = mesh.getFacetType(type);

  Mesh & mesh_facets = mesh.initMeshFacets("mesh_facets");
  MeshUtils::buildAllFacets(mesh, mesh_facets);

  SolidMechanicsModel model(mesh);

  /// model initialization
  model.initFull();

  Array<Real> & position = mesh.getNodes();
  UInt nb_facet = mesh_facets.getNbElement(type_facet);
  UInt nb_element = mesh.getNbElement(type);

  /// compute quadrature points positions on facets
  typedef FEEngineTemplate<IntegratorGauss, ShapeLagrange> MyFEEngineType;

  model.registerFEEngineObject<MyFEEngineType>("FacetsFEEngine", mesh_facets,
                                               spatial_dimension - 1);
  model.getFEEngine("FacetsFEEngine").initShapeFunctions();

  UInt nb_quad_per_facet =
      model.getFEEngine("FacetsFEEngine").getNbIntegrationPoints(type_facet);
  UInt nb_tot_quad = nb_quad_per_facet * nb_facet;

  Array<Real> quad_facets(nb_tot_quad, spatial_dimension);

  model.getFEEngine("FacetsFEEngine")
      .interpolateOnIntegrationPoints(position, quad_facets, spatial_dimension,
                                      type_facet);

  auto && facet_to_element = mesh_facets.getSubelementToElement(type);
  UInt nb_facet_per_elem = facet_to_element.getNbComponent();

  ElementTypeMapArray<Real> element_quad_facet;
  element_quad_facet.alloc(nb_element * nb_facet_per_elem * nb_quad_per_facet,
                           spatial_dimension, type);

  ElementTypeMapArray<Real> interpolated_stress("interpolated_stress", "");
  interpolated_stress.initialize(
      mesh, _nb_component = spatial_dimension * spatial_dimension,
      _spatial_dimension = spatial_dimension);

  Array<Real> & interp_stress = interpolated_stress(type);
  interp_stress.resize(nb_element * nb_facet_per_elem * nb_quad_per_facet);

  Array<Real> & el_q_facet = element_quad_facet(type);

  for (UInt el = 0; el < nb_element; ++el) {
    for (UInt f = 0; f < nb_facet_per_elem; ++f) {
      UInt global_facet = facet_to_element(el, f).element;

      for (UInt q = 0; q < nb_quad_per_facet; ++q) {
        for (UInt s = 0; s < spatial_dimension; ++s) {
          el_q_facet(el * nb_facet_per_elem * nb_quad_per_facet +
                         f * nb_quad_per_facet + q,
                     s) = quad_facets(global_facet * nb_quad_per_facet + q, s);
        }
      }
    }
  }

  /// compute quadrature points position of the elements
  UInt nb_quad_per_element = model.getFEEngine().getNbIntegrationPoints(type);
  UInt nb_tot_quad_el = nb_quad_per_element * nb_element;

  Array<Real> quad_elements(nb_tot_quad_el, spatial_dimension);

  model.getFEEngine().interpolateOnIntegrationPoints(position, quad_elements,
                                                     spatial_dimension, type);

  /// assign some values to stresses
  Array<Real> & stress =
      const_cast<Array<Real> &>(model.getMaterial(0).getStress(type));

  for (UInt q = 0; q < nb_tot_quad_el; ++q) {
    for (UInt s = 0; s < spatial_dimension * spatial_dimension; ++s) {
      stress(q, s) = s * function(quad_elements(q, 0), quad_elements(q, 1),
                                  quad_elements(q, 2));
    }
  }

  /// interpolate stresses on facets' quadrature points
  model.getMaterial(0).initElementalFieldInterpolation(element_quad_facet);
  model.getMaterial(0).interpolateStress(interpolated_stress);

  Real tolerance = 1.e-10;

  /// check results
  for (UInt el = 0; el < nb_element; ++el) {
    for (UInt f = 0; f < nb_facet_per_elem; ++f) {

      for (UInt q = 0; q < nb_quad_per_facet; ++q) {
        for (UInt s = 0; s < spatial_dimension * spatial_dimension; ++s) {

          Real x = el_q_facet(el * nb_facet_per_elem * nb_quad_per_facet +
                                  f * nb_quad_per_facet + q,
                              0);
          Real y = el_q_facet(el * nb_facet_per_elem * nb_quad_per_facet +
                                  f * nb_quad_per_facet + q,
                              1);
          Real z = el_q_facet(el * nb_facet_per_elem * nb_quad_per_facet +
                                  f * nb_quad_per_facet + q,
                              2);

          Real theoretical = s * function(x, y, z);

          Real numerical =
              interp_stress(el * nb_facet_per_elem * nb_quad_per_facet +
                                f * nb_quad_per_facet + q,
                            s);

          if (std::abs(theoretical - numerical) > tolerance) {
            std::cout << "Theoretical and numerical values aren't coincident!"
                      << std::endl;
            return EXIT_FAILURE;
          }
        }
      }
    }
  }

  std::cout << "OK: Stress interpolation test passed." << std::endl;
  return EXIT_SUCCESS;
}

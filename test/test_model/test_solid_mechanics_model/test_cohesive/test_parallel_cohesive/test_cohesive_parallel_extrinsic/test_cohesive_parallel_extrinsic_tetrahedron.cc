/**
 * @file   test_cohesive_parallel_extrinsic_tetrahedron.cc
 *
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Wed Nov 05 2014
 * @date last modification: Wed Nov 08 2017
 *
 * @brief  3D extrinsic cohesive elements test
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
#include "material_cohesive_linear.hh"
#include "solid_mechanics_model_cohesive.hh"

/* -------------------------------------------------------------------------- */
using namespace akantu;

Real function(Real constant, Real x, Real y, Real z) {
  return constant + 2. * x + 3. * y + 4 * z;
}

int main(int argc, char * argv[]) {
  initialize("material.dat", argc, argv);

  debug::setDebugLevel(dblWarning);

  // const UInt max_steps = 1000;
  // Real increment = 0.005;
  const UInt spatial_dimension = 3;
  Math::setTolerance(1.e-12);

  ElementType type = _tetrahedron_10;
  ElementType type_facet = Mesh::getFacetType(type);
  ElementType type_cohesive = FEEngine::getCohesiveElementType(type_facet);

  Mesh mesh(spatial_dimension);

  const auto & comm = Communicator::getStaticCommunicator();
  Int psize = comm.getNbProc();
  Int prank = comm.whoAmI();

  akantu::MeshPartition * partition = NULL;
  if (prank == 0) {
    // Read the mesh
    mesh.read("tetrahedron.msh");

    /// partition the mesh
    partition = new MeshPartitionScotch(mesh, spatial_dimension);
    partition->partitionate(psize);
  }

  SolidMechanicsModelCohesive model(mesh);

  /// model initialization
  model.initParallel(partition, NULL, true);
  model.initFull(
      SolidMechanicsModelCohesiveOptions(_explicit_lumped_mass, true));

  const MaterialCohesiveLinear<3> & mat_cohesive =
      dynamic_cast<const MaterialCohesiveLinear<3> &>(model.getMaterial(1));

  const Real sigma_c =
      mat_cohesive.getParam<RandomInternalField<Real, FacetInternalField>>(
          "sigma_c");
  const Real beta = mat_cohesive.getParam<Real>("beta");
  //  const Real G_cI = mat_cohesive.getParam<Real>("G_cI");

  Array<Real> & position = mesh.getNodes();

  /* ------------------------------------------------------------------------ */
  /* Facet part                                                               */
  /* ------------------------------------------------------------------------ */

  /// compute quadrature points positions on facets
  const Mesh & mesh_facets = model.getMeshFacets();
  UInt nb_facet = mesh_facets.getNbElement(type_facet);
  UInt nb_quad_per_facet =
      model.getFEEngine("FacetsFEEngine").getNbIntegrationPoints(type_facet);
  UInt nb_tot_quad = nb_quad_per_facet * nb_facet;

  Array<Real> quad_facets(nb_tot_quad, spatial_dimension);

  model.getFEEngine("FacetsFEEngine")
      .interpolateOnIntegrationPoints(position, quad_facets, spatial_dimension,
                                      type_facet);

  /* ------------------------------------------------------------------------ */
  /* End of facet part                                                        */
  /* ------------------------------------------------------------------------ */

  /// compute quadrature points position of the elements
  UInt nb_quad_per_element = model.getFEEngine().getNbIntegrationPoints(type);
  UInt nb_element = mesh.getNbElement(type);
  UInt nb_tot_quad_el = nb_quad_per_element * nb_element;

  Array<Real> quad_elements(nb_tot_quad_el, spatial_dimension);

  model.getFEEngine().interpolateOnIntegrationPoints(position, quad_elements,
                                                     spatial_dimension, type);

  /// assign some values to stresses
  Array<Real> & stress =
      const_cast<Array<Real> &>(model.getMaterial(0).getStress(type));

  Array<Real>::iterator<Matrix<Real>> stress_it =
      stress.begin(spatial_dimension, spatial_dimension);

  for (UInt q = 0; q < nb_tot_quad_el; ++q, ++stress_it) {

    /// compute values
    for (UInt i = 0; i < spatial_dimension; ++i) {
      for (UInt j = i; j < spatial_dimension; ++j) {
        UInt index = i * spatial_dimension + j;
        (*stress_it)(i, j) =
            index * function(sigma_c * 5, quad_elements(q, 0),
                             quad_elements(q, 1), quad_elements(q, 2));
      }
    }

    /// fill symmetrical part
    for (UInt i = 0; i < spatial_dimension; ++i) {
      for (UInt j = 0; j < i; ++j) {
        (*stress_it)(i, j) = (*stress_it)(j, i);
      }
    }
  }

  /// compute stress on facet quads
  Array<Real> stress_facets(nb_tot_quad, spatial_dimension * spatial_dimension);

  Array<Real>::iterator<Matrix<Real>> stress_facets_it =
      stress_facets.begin(spatial_dimension, spatial_dimension);

  for (UInt q = 0; q < nb_tot_quad; ++q, ++stress_facets_it) {
    /// compute values
    for (UInt i = 0; i < spatial_dimension; ++i) {
      for (UInt j = i; j < spatial_dimension; ++j) {
        UInt index = i * spatial_dimension + j;
        (*stress_facets_it)(i, j) =
            index * function(sigma_c * 5, quad_facets(q, 0), quad_facets(q, 1),
                             quad_facets(q, 2));
      }
    }

    /// fill symmetrical part
    for (UInt i = 0; i < spatial_dimension; ++i) {
      for (UInt j = 0; j < i; ++j) {
        (*stress_facets_it)(i, j) = (*stress_facets_it)(j, i);
      }
    }
  }

  /// insert cohesive elements
  model.checkCohesiveStress();

  /// check insertion stress
  const Array<Real> & normals = model.getFEEngine("FacetsFEEngine")
                                    .getNormalsOnIntegrationPoints(type_facet);
  const Array<Real> & tangents = model.getTangents(type_facet);
  const Array<Real> & sigma_c_eff =
      mat_cohesive.getInsertionTraction(type_cohesive);

  Vector<Real> normal_stress(spatial_dimension);

  const Array<std::vector<Element>> & coh_element_to_facet =
      mesh_facets.getElementToSubelement(type_facet);

  Array<Real>::iterator<Matrix<Real>> quad_facet_stress =
      stress_facets.begin(spatial_dimension, spatial_dimension);

  Array<Real>::const_iterator<Vector<Real>> quad_normal =
      normals.begin(spatial_dimension);

  Array<Real>::const_iterator<Vector<Real>> quad_tangents =
      tangents.begin(tangents.getNbComponent());

  for (UInt f = 0; f < nb_facet; ++f) {
    const Element & cohesive_element = coh_element_to_facet(f)[1];

    for (UInt q = 0; q < nb_quad_per_facet;
         ++q, ++quad_facet_stress, ++quad_normal, ++quad_tangents) {
      if (cohesive_element == ElementNull)
        continue;

      normal_stress.mul<false>(*quad_facet_stress, *quad_normal);

      Real normal_contrib = normal_stress.dot(*quad_normal);

      Real first_tangent_contrib = 0;

      for (UInt dim = 0; dim < spatial_dimension; ++dim)
        first_tangent_contrib += normal_stress(dim) * (*quad_tangents)(dim);

      Real second_tangent_contrib = 0;

      for (UInt dim = 0; dim < spatial_dimension; ++dim)
        second_tangent_contrib +=
            normal_stress(dim) * (*quad_tangents)(dim + spatial_dimension);

      Real tangent_contrib =
          std::sqrt(first_tangent_contrib * first_tangent_contrib +
                    second_tangent_contrib * second_tangent_contrib);

      normal_contrib = std::max(0., normal_contrib);

      Real effective_norm =
          std::sqrt(normal_contrib * normal_contrib +
                    tangent_contrib * tangent_contrib / beta / beta);

      if (effective_norm < sigma_c)
        continue;

      if (!Math::are_float_equal(
              effective_norm,
              sigma_c_eff(cohesive_element.element * nb_quad_per_facet + q))) {
        std::cout << "Insertion tractions do not match" << std::endl;
        finalize();
        return EXIT_FAILURE;
      }
    }
  }
  finalize();
  return EXIT_SUCCESS;
}

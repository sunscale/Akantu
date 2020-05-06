/**
 * @file   test_cohesive_facet_stress_synchronizer.cc
 *
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Wed Nov 05 2014
 * @date last modification: Wed Nov 08 2017
 *
 * @brief  Test for facet stress synchronizer
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
#include <fstream>
#include <iostream>
#include <limits>

/* -------------------------------------------------------------------------- */
#include "solid_mechanics_model_cohesive.hh"
/* -------------------------------------------------------------------------- */

using namespace akantu;

Real function(Real constant, Real x, Real y, Real z) {
  return constant + 2. * x + 3. * y + 4 * z;
}

int main(int argc, char * argv[]) {
  initialize("material.dat", argc, argv);

  const UInt spatial_dimension = 3;

  ElementType type = _tetrahedron_10;
  ElementType type_facet = Mesh::getFacetType(type);
  Math::setTolerance(1.e-11);

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
  model.initParallel(partition, NULL, true);
  model.initFull(
      SolidMechanicsModelCohesiveOptions(_explicit_lumped_mass, true));

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
        (*stress_it)(i, j) = function(index, quad_elements(q, 0),
                                      quad_elements(q, 1), quad_elements(q, 2));
      }
    }

    /// fill symmetrical part
    for (UInt i = 0; i < spatial_dimension; ++i) {
      for (UInt j = 0; j < i; ++j) {
        (*stress_it)(i, j) = (*stress_it)(j, i);
      }
    }

    // stress_it->clear();
    // for (UInt i = 0; i < spatial_dimension; ++i)
    //   (*stress_it)(i, i) = sigma_c * 5;
  }

  /// compute and communicate stress on facets
  model.checkCohesiveStress();

  /* ------------------------------------------------------------------------ */
  /* Check facet stress                                                       */
  /* ------------------------------------------------------------------------ */

  const Array<Real> & facet_stress = model.getStressOnFacets(type_facet);
  const Array<bool> & facet_check =
      model.getElementInserter().getCheckFacets(type_facet);
  const Array<std::vector<Element>> & elements_to_facet =
      model.getMeshFacets().getElementToSubelement(type_facet);

  Array<Real>::iterator<Vector<Real>> quad_facet_it =
      quad_facets.begin(spatial_dimension);
  Array<Real>::const_iterator<Matrix<Real>> facet_stress_it =
      facet_stress.begin(spatial_dimension, spatial_dimension * 2);

  Matrix<Real> current_stress(spatial_dimension, spatial_dimension);

  for (UInt f = 0; f < nb_facet; ++f) {

    if (!facet_check(f) || (elements_to_facet(f)[0].ghost_type == _not_ghost &&
                            elements_to_facet(f)[1].ghost_type == _not_ghost)) {
      quad_facet_it += nb_quad_per_facet;
      facet_stress_it += nb_quad_per_facet;
      continue;
    }

    for (UInt q = 0; q < nb_quad_per_facet;
         ++q, ++quad_facet_it, ++facet_stress_it) {
      /// compute current_stress
      for (UInt i = 0; i < spatial_dimension; ++i) {
        for (UInt j = i; j < spatial_dimension; ++j) {
          UInt index = i * spatial_dimension + j;
          current_stress(i, j) =
              function(index, (*quad_facet_it)(0), (*quad_facet_it)(1),
                       (*quad_facet_it)(2));
        }
      }

      /// fill symmetrical part
      for (UInt i = 0; i < spatial_dimension; ++i) {
        for (UInt j = 0; j < i; ++j) {
          current_stress(i, j) = current_stress(j, i);
        }
      }

      /// compare it to interpolated one
      for (UInt s = 0; s < 2; ++s) {
        Matrix<Real> stress_to_check(facet_stress_it->storage() +
                                         s * spatial_dimension *
                                             spatial_dimension,
                                     spatial_dimension, spatial_dimension);

        for (UInt i = 0; i < spatial_dimension; ++i) {
          for (UInt j = 0; j < spatial_dimension; ++j) {
            if (!Math::are_float_equal(stress_to_check(i, j),
                                       current_stress(i, j))) {
              std::cout << "Stress doesn't match" << std::endl;
              finalize();
              return EXIT_FAILURE;
            }
          }
        }
      }
    }
  }

  finalize();
  if (prank == 0)
    std::cout << "OK: test_cohesive_facet_stress_synchronizer passed!"
              << std::endl;
  return EXIT_SUCCESS;
}

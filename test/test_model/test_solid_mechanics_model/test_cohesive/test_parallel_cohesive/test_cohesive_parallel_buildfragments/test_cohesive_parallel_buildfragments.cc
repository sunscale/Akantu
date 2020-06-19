/**
 * @file   test_cohesive_parallel_buildfragments.cc
 *
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Wed Nov 05 2014
 * @date last modification: Tue Feb 20 2018
 *
 * @brief  Test to build fragments in parallel
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
#include <algorithm>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
/* -------------------------------------------------------------------------- */
#include "fragment_manager.hh"
#include "material_cohesive.hh"
#include "solid_mechanics_model_cohesive.hh"
/* -------------------------------------------------------------------------- */

using namespace akantu;

void verticalInsertionLimit(SolidMechanicsModelCohesive &);
void displaceElements(SolidMechanicsModelCohesive &, const Real, const Real);
bool isInertiaEqual(const Vector<Real> &, const Vector<Real> &);
void rotateArray(Array<Real> & array, Real angle);
UInt getNbBigFragments(FragmentManager &, UInt);

const UInt spatial_dimension = 3;
const UInt total_nb_fragment = 4;
const Real rotation_angle = M_PI / 4.;
const Real global_tolerance = 1.e-9;

int main(int argc, char * argv[]) {
  initialize("material.dat", argc, argv);

  Math::setTolerance(global_tolerance);

  Mesh mesh(spatial_dimension);

  const auto & comm = Communicator::getStaticCommunicator();
  Int psize = comm.getNbProc();
  Int prank = comm.whoAmI();

  akantu::MeshPartition * partition = NULL;
  if (prank == 0) {
    // Read the mesh
    mesh.read("mesh.msh");

    /// partition the mesh
    MeshUtils::purifyMesh(mesh);
    partition = new MeshPartitionScotch(mesh, spatial_dimension);
    partition->partitionate(psize);
  }

  SolidMechanicsModelCohesive model(mesh);
  model.initParallel(partition, NULL, true);

  delete partition;

  /// model initialization
  model.initFull(
      SolidMechanicsModelCohesiveOptions(_explicit_lumped_mass, true));

  mesh.computeBoundingBox();
  Real L = mesh.getUpperBounds()(0) - mesh.getLowerBounds()(0);
  Real h = mesh.getUpperBounds()(1) - mesh.getLowerBounds()(1);
  Real rho = model.getMaterial("bulk").getParam<Real>("rho");

  Real theoretical_mass = L * h * h * rho;
  Real frag_theo_mass = theoretical_mass / total_nb_fragment;

  UInt nb_element =
      mesh.getNbElement(spatial_dimension, _not_ghost, _ek_regular);
  comm.allReduce(&nb_element, 1, _so_sum);
  UInt nb_element_per_fragment = nb_element / total_nb_fragment;

  FragmentManager fragment_manager(model);

  fragment_manager.computeAllData();
  getNbBigFragments(fragment_manager, nb_element_per_fragment + 1);

  model.setBaseName("extrinsic");
  model.addDumpFieldVector("displacement");
  model.addDumpField("velocity");
  model.addDumpField("stress");
  model.addDumpField("partitions");
  model.addDumpField("fragments");
  model.addDumpField("fragments mass");
  model.addDumpField("moments of inertia");
  model.addDumpField("elements per fragment");
  model.dump();

  model.setBaseNameToDumper("cohesive elements", "cohesive_elements_test");
  model.addDumpFieldVectorToDumper("cohesive elements", "displacement");
  model.addDumpFieldToDumper("cohesive elements", "damage");
  model.dump("cohesive elements");

  /// set check facets
  verticalInsertionLimit(model);

  model.assembleMassLumped();
  model.synchronizeBoundaries();

  /// impose initial displacement
  Array<Real> & displacement = model.getDisplacement();
  Array<Real> & velocity = model.getVelocity();
  const Array<Real> & position = mesh.getNodes();
  UInt nb_nodes = mesh.getNbNodes();

  for (UInt n = 0; n < nb_nodes; ++n) {
    displacement(n, 0) = position(n, 0) * 0.1;
    velocity(n, 0) = position(n, 0);
  }

  rotateArray(mesh.getNodes(), rotation_angle);
  //  rotateArray(displacement, rotation_angle);
  //  rotateArray(velocity, rotation_angle);

  model.updateResidual();
  model.checkCohesiveStress();
  model.dump();
  model.dump("cohesive elements");

  const Array<Real> & fragment_mass = fragment_manager.getMass();
  const Array<Real> & fragment_center = fragment_manager.getCenterOfMass();

  Real el_size = L / total_nb_fragment;
  Real lim = -L / 2 + el_size * 0.99;

  /// define theoretical inertia moments
  Vector<Real> small_frag_inertia(spatial_dimension);
  small_frag_inertia(0) = frag_theo_mass * (h * h + h * h) / 12.;
  small_frag_inertia(1) = frag_theo_mass * (el_size * el_size + h * h) / 12.;
  small_frag_inertia(2) = frag_theo_mass * (el_size * el_size + h * h) / 12.;

  std::sort(small_frag_inertia.storage(),
            small_frag_inertia.storage() + spatial_dimension,
            std::greater<Real>());

  const Array<Real> & inertia_moments = fragment_manager.getMomentsOfInertia();
  Array<Real>::const_iterator<Vector<Real>> inertia_moments_begin =
      inertia_moments.begin(spatial_dimension);

  /// displace one fragment each time
  for (UInt frag = 1; frag <= total_nb_fragment; ++frag) {
    if (prank == 0)
      std::cout << "Generating fragment: " << frag << std::endl;

    fragment_manager.computeAllData();

    /// check number of big fragments
    UInt nb_big_fragment =
        getNbBigFragments(fragment_manager, nb_element_per_fragment + 1);

    model.dump();
    model.dump("cohesive elements");

    if (frag < total_nb_fragment) {
      if (nb_big_fragment != 1) {
        AKANTU_ERROR(
            "The number of big fragments is wrong: " << nb_big_fragment);
        return EXIT_FAILURE;
      }
    } else {
      if (nb_big_fragment != 0) {
        AKANTU_ERROR(
            "The number of big fragments is wrong: " << nb_big_fragment);
        return EXIT_FAILURE;
      }
    }

    /// check number of fragments
    UInt nb_fragment_num = fragment_manager.getNbFragment();

    if (nb_fragment_num != frag) {
      AKANTU_ERROR("The number of fragments is wrong! Numerical: "
                   << nb_fragment_num << " Theoretical: " << frag);
      return EXIT_FAILURE;
    }

    /// check mass computation
    if (frag < total_nb_fragment) {
      Real total_mass = 0.;
      UInt small_fragments = 0;

      for (UInt f = 0; f < nb_fragment_num; ++f) {
        const Vector<Real> & current_inertia = inertia_moments_begin[f];

        if (Math::are_float_equal(fragment_mass(f, 0), frag_theo_mass)) {

          /// check center of mass
          if (fragment_center(f, 0) > (L * frag / total_nb_fragment - L / 2)) {
            AKANTU_ERROR("Fragment center is wrong!");
            return EXIT_FAILURE;
          }

          /// check moment of inertia
          if (!isInertiaEqual(current_inertia, small_frag_inertia)) {
            AKANTU_ERROR("Inertia moments are wrong");
            return EXIT_FAILURE;
          }

          ++small_fragments;
          total_mass += frag_theo_mass;
        } else {
          /// check the moment of inertia for the biggest fragment
          Real big_frag_mass = frag_theo_mass * (total_nb_fragment - frag + 1);
          Real big_frag_size = el_size * (total_nb_fragment - frag + 1);

          Vector<Real> big_frag_inertia(spatial_dimension);
          big_frag_inertia(0) = big_frag_mass * (h * h + h * h) / 12.;
          big_frag_inertia(1) =
              big_frag_mass * (big_frag_size * big_frag_size + h * h) / 12.;
          big_frag_inertia(2) =
              big_frag_mass * (big_frag_size * big_frag_size + h * h) / 12.;

          std::sort(big_frag_inertia.storage(),
                    big_frag_inertia.storage() + spatial_dimension,
                    std::greater<Real>());

          if (!isInertiaEqual(current_inertia, big_frag_inertia)) {
            AKANTU_ERROR("Inertia moments are wrong");
            return EXIT_FAILURE;
          }
        }
      }

      if (small_fragments != nb_fragment_num - 1) {
        AKANTU_ERROR("The number of small fragments is wrong!");
        return EXIT_FAILURE;
      }

      if (!Math::are_float_equal(total_mass,
                                 small_fragments * frag_theo_mass)) {
        AKANTU_ERROR("The mass of small fragments is wrong!");
        return EXIT_FAILURE;
      }
    }

    /// displace fragments
    rotateArray(mesh.getNodes(), -rotation_angle);
    //    rotateArray(displacement, -rotation_angle);
    displaceElements(model, lim, el_size * 2);
    rotateArray(mesh.getNodes(), rotation_angle);
    //    rotateArray(displacement, rotation_angle);

    model.updateResidual();

    lim += el_size;
  }

  model.dump();
  model.dump("cohesive elements");

  /// check centers
  const Array<Real> & fragment_velocity = fragment_manager.getVelocity();

  Real initial_position = -L / 2. + el_size / 2.;

  for (UInt frag = 0; frag < total_nb_fragment; ++frag) {
    Real theoretical_center = initial_position + el_size * frag;

    UInt f_index = 0;
    while (
        f_index < total_nb_fragment &&
        !Math::are_float_equal(fragment_center(f_index, 0), theoretical_center))
      ++f_index;

    if (f_index == total_nb_fragment) {
      AKANTU_ERROR("The fragments' center is wrong!");
      return EXIT_FAILURE;
    }

    f_index = 0;
    while (f_index < total_nb_fragment &&
           !Math::are_float_equal(fragment_velocity(f_index, 0),
                                  theoretical_center))
      ++f_index;

    if (f_index == total_nb_fragment) {
      AKANTU_ERROR("The fragments' velocity is wrong!");
      return EXIT_FAILURE;
    }
  }

  finalize();

  if (prank == 0)
    std::cout << "OK: test_cohesive_buildfragments was passed!" << std::endl;
  return EXIT_SUCCESS;
}

void verticalInsertionLimit(SolidMechanicsModelCohesive & model) {
  UInt spatial_dimension = model.getSpatialDimension();
  const Mesh & mesh_facets = model.getMeshFacets();
  const Array<Real> & position = mesh_facets.getNodes();

  for (ghost_type_t::iterator gt = ghost_type_t::begin();
       gt != ghost_type_t::end(); ++gt) {
    GhostType ghost_type = *gt;

    Mesh::type_iterator it =
        mesh_facets.firstType(spatial_dimension - 1, ghost_type);
    Mesh::type_iterator end =
        mesh_facets.lastType(spatial_dimension - 1, ghost_type);
    for (; it != end; ++it) {
      ElementType type = *it;

      Array<bool> & check_facets =
          model.getElementInserter().getCheckFacets(type, ghost_type);

      const Array<UInt> & connectivity =
          mesh_facets.getConnectivity(type, ghost_type);
      UInt nb_nodes_per_facet = connectivity.getNbComponent();

      for (UInt f = 0; f < check_facets.getSize(); ++f) {
        if (!check_facets(f))
          continue;

        UInt nb_aligned_nodes = 1;
        Real first_node_pos = position(connectivity(f, 0), 0);

        for (; nb_aligned_nodes < nb_nodes_per_facet; ++nb_aligned_nodes) {
          Real other_node_pos = position(connectivity(f, nb_aligned_nodes), 0);

          if (!Math::are_float_equal(first_node_pos, other_node_pos))
            break;
        }

        if (nb_aligned_nodes != nb_nodes_per_facet) {
          check_facets(f) = false;
        }
      }
    }
  }
}

void displaceElements(SolidMechanicsModelCohesive & model, const Real lim,
                      const Real amount) {
  UInt spatial_dimension = model.getSpatialDimension();
  Array<Real> & displacement = model.getDisplacement();
  Mesh & mesh = model.getMesh();
  UInt nb_nodes = mesh.getNbNodes();
  Array<bool> displaced(nb_nodes);
  displaced.clear();
  Vector<Real> barycenter(spatial_dimension);

  for (ghost_type_t::iterator gt = ghost_type_t::begin();
       gt != ghost_type_t::end(); ++gt) {
    GhostType ghost_type = *gt;

    Mesh::type_iterator it = mesh.firstType(spatial_dimension, ghost_type);
    Mesh::type_iterator end = mesh.lastType(spatial_dimension, ghost_type);
    for (; it != end; ++it) {
      ElementType type = *it;
      const Array<UInt> & connectivity = mesh.getConnectivity(type, ghost_type);
      UInt nb_element = connectivity.getSize();
      UInt nb_nodes_per_element = connectivity.getNbComponent();

      Array<UInt>::const_vector_iterator conn_el =
          connectivity.begin(nb_nodes_per_element);

      for (UInt el = 0; el < nb_element; ++el) {
        mesh.getBarycenter(el, type, barycenter.storage(), ghost_type);

        if (barycenter(0) < lim) {
          const Vector<UInt> & conn = conn_el[el];

          for (UInt n = 0; n < nb_nodes_per_element; ++n) {
            UInt node = conn(n);
            if (!displaced(node)) {
              displacement(node, 0) -= amount;
              displaced(node) = true;
            }
          }
        }
      }
    }
  }
}

bool isInertiaEqual(const Vector<Real> & a, const Vector<Real> & b) {
  UInt nb_terms = a.size();
  UInt equal_terms = 0;

  while (equal_terms < nb_terms &&
         std::abs(a(equal_terms) - b(equal_terms)) / a(equal_terms) <
             Math::getTolerance())
    ++equal_terms;

  return equal_terms == nb_terms;
}

void rotateArray(Array<Real> & array, Real angle) {
  UInt spatial_dimension = array.getNbComponent();

  Real rotation_values[] = {std::cos(angle),
                            std::sin(angle),
                            0,
                            -std::sin(angle),
                            std::cos(angle),
                            0,
                            0,
                            0,
                            1};
  Matrix<Real> rotation(rotation_values, spatial_dimension, spatial_dimension);

  RVector displaced_node(spatial_dimension);
  auto node_it = array.begin(spatial_dimension);
  auto node_end = array.end(spatial_dimension);

  for (; node_it != node_end; ++node_it) {
    displaced_node.mul<false>(rotation, *node_it);
    *node_it = displaced_node;
  }
}

UInt getNbBigFragments(FragmentManager & fragment_manager,
                       UInt minimum_nb_elements) {
  fragment_manager.computeNbElementsPerFragment();
  const Array<UInt> & nb_elements_per_fragment =
      fragment_manager.getNbElementsPerFragment();
  UInt nb_fragment = fragment_manager.getNbFragment();
  UInt nb_big_fragment = 0;

  for (UInt frag = 0; frag < nb_fragment; ++frag) {
    if (nb_elements_per_fragment(frag) >= minimum_nb_elements) {
      ++nb_big_fragment;
    }
  }

  return nb_big_fragment;
}

/**
 * @file   cohesive_element_inserter.cc
 *
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Wed Dec 04 2013
 * @date last modification: Mon Feb 19 2018
 *
 * @brief  Cohesive element inserter functions
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
#include "cohesive_element_inserter.hh"
#include "cohesive_element_inserter_helper.hh"
#include "communicator.hh"
#include "element_group.hh"
#include "element_synchronizer.hh"
#include "global_ids_updater.hh"
#include "mesh_accessor.hh"
#include "mesh_iterators.hh"
/* -------------------------------------------------------------------------- */
#include <algorithm>
#include <limits>
/* -------------------------------------------------------------------------- */

namespace akantu {

CohesiveElementInserter::CohesiveElementInserter(Mesh & mesh, const ID & id)
    : Parsable(ParserType::_cohesive_inserter), id(id), mesh(mesh),
      mesh_facets(mesh.initMeshFacets()),
      insertion_facets("insertion_facets", id),
      insertion_limits(mesh.getSpatialDimension(), 2),
      check_facets("check_facets", id) {

  this->registerParam("cohesive_surfaces", physical_surfaces, _pat_parsable,
                      "List of groups to consider for insertion");
  this->registerParam("cohesive_zones", physical_zones, _pat_parsable,
                      "List of groups to consider for insertion");
  this->registerParam("bounding_box", insertion_limits, _pat_parsable,
                      "Global limit for insertion");

  UInt spatial_dimension = mesh.getSpatialDimension();

  /// init insertion limits
  for (UInt dim = 0; dim < spatial_dimension; ++dim) {
    insertion_limits(dim, 0) = std::numeric_limits<Real>::max() * Real(-1.);
    insertion_limits(dim, 1) = std::numeric_limits<Real>::max();
  }

  insertion_facets.initialize(mesh_facets,
                              _spatial_dimension = spatial_dimension - 1,
                              _with_nb_element = true, _default_value = false);
}

/* -------------------------------------------------------------------------- */
CohesiveElementInserter::~CohesiveElementInserter() = default;

/* -------------------------------------------------------------------------- */
void CohesiveElementInserter::parseSection(const ParserSection & section) {
  Parsable::parseSection(section);

  if (is_extrinsic) {
    limitCheckFacets(this->check_facets);
  }
}

/* -------------------------------------------------------------------------- */
void CohesiveElementInserter::limitCheckFacets() {
  limitCheckFacets(this->check_facets);
}

/* -------------------------------------------------------------------------- */
void CohesiveElementInserter::setLimit(SpatialDirection axis, Real first_limit,
                                       Real second_limit) {
  AKANTU_DEBUG_ASSERT(
      axis < mesh.getSpatialDimension(),
      "You are trying to limit insertion in a direction that doesn't exist");

  insertion_limits(axis, 0) = std::min(first_limit, second_limit);
  insertion_limits(axis, 1) = std::max(first_limit, second_limit);
}

/* -------------------------------------------------------------------------- */
UInt CohesiveElementInserter::insertIntrinsicElements() {
  limitCheckFacets(insertion_facets);
  return insertElements();
}

/* -------------------------------------------------------------------------- */
void CohesiveElementInserter::limitCheckFacets(
    ElementTypeMapArray<bool> & check_facets) {
  AKANTU_DEBUG_IN();

  UInt spatial_dimension = mesh.getSpatialDimension();

  check_facets.initialize(mesh_facets,
                          _spatial_dimension = spatial_dimension - 1,
                          _with_nb_element = true, _default_value = true);
  check_facets.set(true);

  // remove the pure ghost elements
  for_each_element(
      mesh_facets,
      [&](auto && facet) {
        const auto & element_to_facet = mesh_facets.getElementToSubelement(
            facet.type, facet.ghost_type)(facet.element);
        auto & left = element_to_facet[0];
        auto & right = element_to_facet[1];
        if (right == ElementNull ||
            (left.ghost_type == _ghost && right.ghost_type == _ghost)) {
          check_facets(facet) = false;
          return;
        }
#ifndef AKANTU_NDEBUG
        if (left == ElementNull) {
          AKANTU_DEBUG_WARNING("By convention element should not have "
                               "ElementNull on there first side: "
                               << facet);
        }
#endif

        if (left.kind() == _ek_cohesive or right.kind() == _ek_cohesive) {
          check_facets(facet) = false;
        }
      },
      _spatial_dimension = spatial_dimension - 1);

  auto tolerance = Math::getTolerance();
  Vector<Real> bary_facet(spatial_dimension);
  // set the limits to the bounding box
  for_each_element(
      mesh_facets,
      [&](auto && facet) {
        auto & need_check = check_facets(facet, 0);
        if (not need_check) {
          return;
        }

        mesh_facets.getBarycenter(facet, bary_facet);
        UInt coord_in_limit = 0;

        while (coord_in_limit < spatial_dimension and
               bary_facet(coord_in_limit) >
                   (insertion_limits(coord_in_limit, 0) - tolerance) and
               bary_facet(coord_in_limit) <
                   (insertion_limits(coord_in_limit, 1) + tolerance)) {
          ++coord_in_limit;
        }

        if (coord_in_limit != spatial_dimension) {
          need_check = false;
        }
      },
      _spatial_dimension = spatial_dimension - 1);

  // remove the physical zones
  if (mesh.hasData("physical_names") and not physical_zones.empty()) {
    auto && physical_names = mesh.getData<std::string>("physical_names");
    for_each_element(
        mesh_facets,
        [&](auto && facet) {
          const auto & element_to_facet = mesh_facets.getElementToSubelement(
              facet.type, facet.ghost_type)(facet.element);
          auto count = 0;
          for (auto i : arange(2)) {
            const auto & element = element_to_facet[i];
            if (element == ElementNull) {
              continue;
            }
            const auto & name = physical_names(element);
            count += find(physical_zones.begin(), physical_zones.end(), name) !=
                     physical_zones.end();
          }

          if (count != 2) {
            check_facets(facet) = false;
          }
        },
        _spatial_dimension = spatial_dimension - 1);
  }

  if (physical_surfaces.empty()) {
    AKANTU_DEBUG_OUT();
    return;
  }

  if (not mesh_facets.hasData("physical_names")) {
    AKANTU_DEBUG_ASSERT(
        physical_surfaces.empty(),
        "No physical names in the mesh but insertion limited to a group");
    AKANTU_DEBUG_OUT();
    return;
  }

  const auto & physical_ids =
      mesh_facets.getData<std::string>("physical_names");

  // set the limits to the physical surfaces
  for_each_element(
      mesh_facets,
      [&](auto && facet) {
        auto & need_check = check_facets(facet, 0);
        if (not need_check) {
          return;
        }

        const auto & physical_id = physical_ids(facet);
        auto it = find(physical_surfaces.begin(), physical_surfaces.end(),
                       physical_id);

        need_check = (it != physical_surfaces.end());
      },
      _spatial_dimension = spatial_dimension - 1);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
UInt CohesiveElementInserter::insertElements(bool only_double_facets) {
  CohesiveNewNodesEvent node_event(AKANTU_CURRENT_FUNCTION);
  NewElementsEvent element_event(AKANTU_CURRENT_FUNCTION);

  if (mesh_facets.isDistributed()) {
    mesh_facets.getElementSynchronizer().synchronizeOnce(
        *this, SynchronizationTag::_ce_groups);
  }

  CohesiveElementInserterHelper cohesive_element_inserter_helper(
      mesh, insertion_facets);

  UInt nb_new_elements{0};
  if (only_double_facets) {
    nb_new_elements = cohesive_element_inserter_helper.insertFacetsOnly();
  } else {
    nb_new_elements = cohesive_element_inserter_helper.insertCohesiveElement();
    element_event.getList().copy(
        cohesive_element_inserter_helper.getNewElements());
  }

  auto && doubled_nodes = cohesive_element_inserter_helper.getDoubledNodes();
  auto nb_new_nodes = doubled_nodes.size();

  node_event.getList().reserve(nb_new_nodes);
  node_event.getOldNodesList().reserve(nb_new_nodes);

  for (auto && doubled_node : make_view(doubled_nodes, 2)) {
    node_event.getList().push_back(doubled_node(1));
    node_event.getOldNodesList().push_back(doubled_node(0));
  }

  if (nb_new_elements > 0) {
    updateInsertionFacets();
  }

  MeshAccessor mesh_accessor(mesh);
  std::tie(nb_new_nodes, nb_new_elements) =
      mesh_accessor.updateGlobalData(node_event, element_event);

  return nb_new_elements;
}

/* -------------------------------------------------------------------------- */
void CohesiveElementInserter::updateInsertionFacets() {
  AKANTU_DEBUG_IN();

  UInt spatial_dimension = mesh.getSpatialDimension();

  for (auto && facet_gt : ghost_types) {
    for (auto && facet_type :
         mesh_facets.elementTypes(spatial_dimension - 1, facet_gt)) {
      auto & ins_facets = insertion_facets(facet_type, facet_gt);

      // this is the intrinsic case
      if (not is_extrinsic) {
        continue;
      }

      auto & f_check = check_facets(facet_type, facet_gt);
      for (auto && pair : zip(ins_facets, f_check)) {
        bool & ins = std::get<0>(pair);
        bool & check = std::get<1>(pair);
        if (ins) {
          ins = check = false;
        }
      }
    }
  }

  // resize for the newly added facets
  insertion_facets.initialize(mesh_facets,
                              _spatial_dimension = spatial_dimension - 1,
                              _with_nb_element = true, _default_value = false);

  // resize for the newly added facets
  if (is_extrinsic) {
    check_facets.initialize(mesh_facets,
                            _spatial_dimension = spatial_dimension - 1,
                            _with_nb_element = true, _default_value = false);
  } else {
    insertion_facets.set(false);
  }

  AKANTU_DEBUG_OUT();
}

} // namespace akantu

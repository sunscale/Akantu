/**
 * @file   cohesive_element_inserter.cc
 *
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Wed Dec 04 2013
 * @date last modification: Sun Oct 04 2015
 *
 * @brief  Cohesive element inserter functions
 *
 * @section LICENSE
 *
 * Copyright  (©)  2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de Lausanne)
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
#include "cohesive_element_inserter.hh"
#include "communicator.hh"
#include "element_group.hh"
#include "global_ids_updater.hh"
#include "mesh_iterators.hh"
/* -------------------------------------------------------------------------- */
#include <algorithm>
#include <limits>
/* -------------------------------------------------------------------------- */

namespace akantu {

CohesiveElementInserter::CohesiveElementInserter(Mesh & mesh, bool is_extrinsic,
                                                 const ID & id)
    : Parsable(ParserType::_cohesive_inserter), id(id), mesh(mesh),
      mesh_facets(mesh.initMeshFacets()),
      insertion_facets("insertion_facets", id),
      insertion_limits(mesh.getSpatialDimension(), 2),
      check_facets("check_facets", id) {

  this->registerParam("groups", physical_groups, _pat_parsable,
                      "List of groups to consider for insertion");
  this->registerParam("bounding_box", insertion_limits, _pat_parsable,
                      "Global limit for insertion");

  init(is_extrinsic);
}

/* -------------------------------------------------------------------------- */
CohesiveElementInserter::~CohesiveElementInserter() = default;

/* -------------------------------------------------------------------------- */
void CohesiveElementInserter::init(bool is_extrinsic) {
  AKANTU_DEBUG_IN();

  UInt spatial_dimension = mesh.getSpatialDimension();

  MeshUtils::resetFacetToDouble(mesh_facets);

  /// initialize facet insertion array
  insertion_facets.initialize(mesh_facets,
                              _spatial_dimension = (spatial_dimension - 1),
                              _with_nb_element = true);

  /// init insertion limits
  for (UInt dim = 0; dim < spatial_dimension; ++dim) {
    insertion_limits(dim, 0) = std::numeric_limits<Real>::max() * (-1.);
    insertion_limits(dim, 1) = std::numeric_limits<Real>::max();
  }

  if (is_extrinsic) {
    initFacetsCheck();
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void CohesiveElementInserter::initFacetsCheck() {
  AKANTU_DEBUG_IN();

  UInt spatial_dimension = mesh.getSpatialDimension();

  check_facets.initialize(mesh_facets,
                          _spatial_dimension = spatial_dimension - 1,
                          _with_nb_element = true,
                          _default_value = true);

  for_each_elements(
      mesh_facets,
      [&](auto && facet) {
        const auto & element_to_facet = mesh_facets.getElementToSubelement(
            facet.type, facet.ghost_type)(facet.element);
        auto & left = element_to_facet[0];
        auto & right = element_to_facet[1];
        if (right == ElementNull ||
            (left.ghost_type == _ghost && right.ghost_type == _ghost)) {
          check_facets(facet) = false;
        }
      },
      _spatial_dimension = spatial_dimension - 1);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void CohesiveElementInserter::limitCheckFacets() {
  AKANTU_DEBUG_IN();

  UInt spatial_dimension = mesh.getSpatialDimension();
  Vector<Real> bary_facet(spatial_dimension);

  for_each_elements(mesh_facets,
                    [&](auto && facet) {
                      auto & need_check = check_facets(facet);
                      if (not need_check)
                        return;

                      mesh_facets.getBarycenter(facet, bary_facet);
                      UInt coord_in_limit = 0;

                      while (coord_in_limit < spatial_dimension &&
                             bary_facet(coord_in_limit) >
                                 insertion_limits(coord_in_limit, 0) &&
                             bary_facet(coord_in_limit) <
                                 insertion_limits(coord_in_limit, 1))
                        ++coord_in_limit;

                      if (coord_in_limit != spatial_dimension)
                        need_check = false;
                    },
                    _spatial_dimension = spatial_dimension - 1);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void CohesiveElementInserter::setLimit(SpacialDirection axis, Real first_limit,
                                       Real second_limit) {
  AKANTU_DEBUG_ASSERT(
      axis < mesh.getSpatialDimension(),
      "You are trying to limit insertion in a direction that doesn't exist");

  insertion_limits(axis, 0) = std::min(first_limit, second_limit);
  insertion_limits(axis, 1) = std::max(first_limit, second_limit);
}

/* -------------------------------------------------------------------------- */
UInt CohesiveElementInserter::insertIntrinsicElements() {
  AKANTU_DEBUG_IN();

  UInt spatial_dimension = mesh.getSpatialDimension();

  Vector<Real> bary_facet(spatial_dimension);

  for_each_elements(
      mesh_facets,
      [&](auto && facet) {
        const auto & element_to_facet = mesh_facets.getElementToSubelement(
            facet.type, facet.ghost_type)(facet.element);
        if (element_to_facet[1] == ElementNull)
          return;

        mesh_facets.getBarycenter(facet, bary_facet);
        UInt coord_in_limit = 0;

        while (coord_in_limit < spatial_dimension &&
               bary_facet(coord_in_limit) >
                   insertion_limits(coord_in_limit, 0) &&
               bary_facet(coord_in_limit) < insertion_limits(coord_in_limit, 1))
          ++coord_in_limit;

        if (coord_in_limit == spatial_dimension)
          insertion_facets(facet) = true;
      },
      _spatial_dimension = spatial_dimension - 1);

  AKANTU_DEBUG_OUT();

  return insertElements();
}

/* -------------------------------------------------------------------------- */
void CohesiveElementInserter::insertIntrinsicElements(
    const std::string & physname, UInt material_index) {
  AKANTU_DEBUG_IN();

  UInt spatial_dimension = mesh.getSpatialDimension();
  ElementTypeMapArray<UInt> * phys_data;
  try {
    phys_data = &(mesh_facets.getData<UInt>("physical_names"));
  } catch (...) {
    phys_data = &(mesh_facets.registerData<UInt>("physical_names"));
    phys_data->initialize(mesh_facets,
                          _spatial_dimension = spatial_dimension - 1,
                          _with_nb_element = true);
    // mesh_facets.initElementTypeMapArray(*phys_data, 1, spatial_dimension -
    // 1,
    //                                     false, _ek_regular, true);
  }
  Vector<Real> bary_facet(spatial_dimension);
  mesh_facets.createElementGroup(physname);

  GhostType ghost_type = _not_ghost;

  for (auto && type_facet :
       mesh_facets.elementTypes(spatial_dimension - 1, ghost_type)) {
    auto & f_insertion = insertion_facets(type_facet, ghost_type);
    auto & element_to_facet =
        mesh_facets.getElementToSubelement(type_facet, ghost_type);

    UInt nb_facet = mesh_facets.getNbElement(type_facet, ghost_type);
    UInt coord_in_limit = 0;

    auto & group = mesh.getElementGroup(physname);
    auto & group_facet = mesh_facets.getElementGroup(physname);

    Vector<Real> bary_physgroup(spatial_dimension);
    Real norm_bary;
    for (auto && e : group.getElements(type_facet, ghost_type)) {
      mesh.getBarycenter(Element{type_facet, e, ghost_type}, bary_physgroup);
#ifndef AKANTU_NDEBUG
      bool find_a_partner = false;
#endif
      norm_bary = bary_physgroup.norm();
      Array<UInt> & material_id = (*phys_data)(type_facet, ghost_type);

      for (UInt f = 0; f < nb_facet; ++f) {

        if (element_to_facet(f)[1] == ElementNull)
          continue;

        mesh_facets.getBarycenter(Element{type_facet, f, ghost_type},
                                  bary_facet);

        coord_in_limit = 0;

        while (coord_in_limit < spatial_dimension &&
               (std::abs(bary_facet(coord_in_limit) -
                         bary_physgroup(coord_in_limit)) /
                    norm_bary <
                Math::getTolerance()))
          ++coord_in_limit;

        if (coord_in_limit == spatial_dimension) {
          f_insertion(f) = true;
#ifndef AKANTU_NDEBUG
          find_a_partner = true;
#endif
          group_facet.add(type_facet, f, ghost_type, false);
          material_id(f) = material_index;
          break;
        }
      }
      AKANTU_DEBUG_ASSERT(find_a_partner,
                          "The element nO "
                              << e << " of physical group " << physname
                              << " did not find its associated facet!"
                              << " Try to decrease math tolerance. "
                              << std::endl);
    }
  }
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
UInt CohesiveElementInserter::insertElements(bool only_double_facets) {

  CohesiveNewNodesEvent node_event;
  NewElementsEvent element_event;

  Array<UInt> new_pairs(0, 2);

  UInt nb_new_elements = MeshUtils::insertCohesiveElements(
      mesh, mesh_facets, insertion_facets, new_pairs, element_event.getList(),
      only_double_facets);

  UInt nb_new_nodes = new_pairs.size();
  node_event.getList().reserve(nb_new_nodes);
  node_event.getOldNodesList().reserve(nb_new_nodes);
  for (UInt n = 0; n < nb_new_nodes; ++n) {
    node_event.getList().push_back(new_pairs(n, 1));
    node_event.getOldNodesList().push_back(new_pairs(n, 0));
  }

  if (mesh.isDistributed()) {

    /// update nodes type
    updateNodesType(mesh, node_event);
    updateNodesType(mesh_facets, node_event);

    /// update global ids
    nb_new_nodes = this->updateGlobalIDs(node_event);

    /// compute total number of new elements
    const auto & comm = mesh.getCommunicator();
    comm.allReduce(nb_new_elements, SynchronizerOperation::_sum);
  }

  if (nb_new_nodes > 0)
    mesh.sendEvent(node_event);

  if (nb_new_elements > 0) {
    updateInsertionFacets();
    // mesh.updateTypesOffsets(_not_ghost);
    mesh.sendEvent(element_event);
    MeshUtils::resetFacetToDouble(mesh_facets);
  }

  if (mesh.isDistributed()) {
    /// update global ids
    this->synchronizeGlobalIDs(node_event);
  }

  return nb_new_elements;
}

/* -------------------------------------------------------------------------- */
void CohesiveElementInserter::updateInsertionFacets() {
  AKANTU_DEBUG_IN();

  UInt spatial_dimension = mesh.getSpatialDimension();

  for (ghost_type_t::iterator gt = ghost_type_t::begin();
       gt != ghost_type_t::end(); ++gt) {

    GhostType facet_gt = *gt;
    Mesh::type_iterator it =
        mesh_facets.firstType(spatial_dimension - 1, facet_gt);
    Mesh::type_iterator last =
        mesh_facets.lastType(spatial_dimension - 1, facet_gt);

    for (; it != last; ++it) {
      ElementType facet_type = *it;

      Array<bool> & ins_facets = insertion_facets(facet_type, facet_gt);

      // this is the extrinsic case
      if (check_facets.exists(facet_type, facet_gt)) {
        Array<bool> & f_check = check_facets(facet_type, facet_gt);

        UInt nb_facets = f_check.size();

        for (UInt f = 0; f < ins_facets.size(); ++f) {
          if (ins_facets(f)) {
            ++nb_facets;
            ins_facets(f) = false;
            f_check(f) = false;
          }
        }

        f_check.resize(nb_facets);
      } else { // and this the intrinsic one
        ins_facets.resize(mesh_facets.getNbElement(facet_type, facet_gt));
        ins_facets.set(false);
      }
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void CohesiveElementInserter::printself(std::ostream & stream,
                                        int indent) const {
  std::string space;
  for (Int i = 0; i < indent; i++, space += AKANTU_INDENT)
    ;

  stream << space << "CohesiveElementInserter [" << std::endl;

  stream << space << " + mesh [" << std::endl;
  mesh.printself(stream, indent + 2);
  stream << space << AKANTU_INDENT << "]" << std::endl;

  stream << space << " + mesh_facets [" << std::endl;
  mesh_facets.printself(stream, indent + 2);
  stream << space << AKANTU_INDENT << "]" << std::endl;

  stream << space << "]" << std::endl;
}

} // namespace akantu

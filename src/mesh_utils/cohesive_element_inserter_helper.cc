/**
 * @file   cohesive_element_inserter_helper.cc
 *
 * @author Nicolas Richart
 *
 * @date creation  jeu sep 03 2020
 *
 * @brief A Documented file.
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
#include "cohesive_element_inserter_helper.hh"
/* -------------------------------------------------------------------------- */
#include "data_accessor.hh"
#include "element_synchronizer.hh"
#include "fe_engine.hh"
#include "mesh_accessor.hh"
/* -------------------------------------------------------------------------- */
#include <queue>
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
CohesiveElementInserterHelper::CohesiveElementInserterHelper(
    Mesh & mesh, const ElementTypeMapArray<bool> & facet_insertion)
    : doubled_nodes(0, 2, "doubled_nodes"), mesh(mesh),
      mesh_facets(mesh.getMeshFacets()) {

  auto spatial_dimension = mesh_facets.getSpatialDimension();

  for (auto gt : ghost_types) {
    for (auto type : mesh_facets.elementTypes(_ghost_type = gt)) {
      nb_new_facets(type, gt) = mesh_facets.getNbElement(type, gt);
    }
  }

  std::array<Int, 2> nb_facet_to_insert{0, 0};
  // creates a vector of the facets to insert
  std::vector<Element> potential_facets_to_double;
  for (auto gt_facet : ghost_types) {
    for (auto type_facet :
         mesh_facets.elementTypes(spatial_dimension - 1, gt_facet)) {
      const auto & f_insertion = facet_insertion(type_facet, gt_facet);
      auto & counter = nb_facet_to_insert[gt_facet == _not_ghost ? 0 : 1];

      for (auto && data : enumerate(f_insertion)) {
        if (std::get<1>(data)) {
          UInt el = std::get<0>(data);
          potential_facets_to_double.push_back({type_facet, el, gt_facet});
          ++counter;
        }
      }
    }
  }

  // Defines a global order of insertion oof cohesive element to ensure node
  // doubling appens in the smae order, this is necessary for the global node
  // numbering
  if (mesh_facets.isDistributed()) {
    const auto & comm = mesh_facets.getCommunicator();

    ElementTypeMapArray<Int> global_orderings;

    global_orderings.initialize(mesh_facets,
                                _spatial_dimension = spatial_dimension - 1,
                                _with_nb_element = true, _default_value = -1);

    auto starting_index = nb_facet_to_insert[0];
    comm.exclusiveScan(starting_index);

    // define the local numbers for facet to insert
    for (auto gt_facet : ghost_types) {
      for (auto type_facet :
           mesh_facets.elementTypes(spatial_dimension - 1, gt_facet)) {
        for (auto data : zip(facet_insertion(type_facet, gt_facet),
                             global_orderings(type_facet, gt_facet))) {
          if (std::get<0>(data)) {
            std::get<1>(data) = starting_index;
            ++starting_index;
          }
        }
      }
    }

    // retreives the oorder number from neighoring processors
    auto && synchronizer = mesh_facets.getElementSynchronizer();
    SimpleElementDataAccessor<Int> data_accessor(
        global_orderings, SynchronizationTag::_ce_insertion_order);
    synchronizer.synchronizeOnce(data_accessor,
                                 SynchronizationTag::_ce_insertion_order);

    // sort the facets to double based on this global ordering
    std::sort(potential_facets_to_double.begin(),
              potential_facets_to_double.end(),
              [&global_orderings](auto && el1, auto && el2) {
                return global_orderings(el1) < global_orderings(el2);
              });
  }

  for (auto d : arange(spatial_dimension)) {
    facets_to_double_by_dim[d] = std::make_unique<Array<Element>>(
        0, 2, "facets_to_double_" + std::to_string(d));
  }

  auto & facets_to_double = *facets_to_double_by_dim[spatial_dimension - 1];

  MeshAccessor mesh_accessor(mesh_facets);
  auto & elements_to_subelements = mesh_accessor.getElementToSubelement();

  for (auto && facet_to_double : potential_facets_to_double) {
    auto gt_facet = facet_to_double.ghost_type;
    auto type_facet = facet_to_double.type;
    auto & elements_to_facets = elements_to_subelements(type_facet, gt_facet);

    auto & elements_to_facet = elements_to_facets(facet_to_double.element);
    if (elements_to_facet[1].type == _not_defined
#if defined(AKANTU_COHESIVE_ELEMENT)
        || elements_to_facet[1].kind() == _ek_cohesive
#endif
    ) {
      AKANTU_DEBUG_WARNING("attempt to double a facet on the boundary");
      continue;
    }

    auto new_facet = nb_new_facets(type_facet, gt_facet)++;

    facets_to_double.push_back(Vector<Element>{
        facet_to_double, Element{type_facet, new_facet, gt_facet}});

    /// update facet_to_element vector
    auto & element_to_update = elements_to_facet[1];
    auto el = element_to_update.element;

    auto & facets_to_elements = mesh_facets.getSubelementToElement(
        element_to_update.type, element_to_update.ghost_type);
    auto facets_to_element = Vector<Element>(
        make_view(facets_to_elements, facets_to_elements.getNbComponent())
            .begin()[el]);
    auto facet_to_update = std::find(facets_to_element.begin(),
                                     facets_to_element.end(), facet_to_double);

    AKANTU_DEBUG_ASSERT(facet_to_update != facets_to_element.end(),
                        "Facet not found");

    facet_to_update->element = new_facet;

    /// update elements connected to facet
    const auto & first_facet_list = elements_to_facet;
    elements_to_facets.push_back(first_facet_list);

    /// set new and original facets as boundary facets
    elements_to_facets(new_facet)[0] = elements_to_facets(new_facet)[1];
    elements_to_facets(new_facet)[1] = ElementNull;

    elements_to_facets(facet_to_double.element)[1] = ElementNull;
  }
}

/* -------------------------------------------------------------------------- */
inline bool
CohesiveElementInserterHelper::hasElement(const Vector<UInt> & nodes_element,
                                          const Vector<UInt> & nodes) {
  // if one of the nodes of "nodes" is not in "nodes_element" it stops
  auto it = std::mismatch(nodes.begin(), nodes.end(), nodes_element.begin(),
                          [&](auto && node, auto && /*node2*/) -> bool {
                            auto it = std::find(nodes_element.begin(),
                                                nodes_element.end(), node);
                            return (it != nodes_element.end());
                          });

  // true if all "nodes" where found in "nodes_element"
  return (it.first == nodes.end());
}

/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
inline bool CohesiveElementInserterHelper::removeElementsInVector(
    const std::vector<Element> & elem_to_remove,
    std::vector<Element> & elem_list) {
  if (elem_list.size() <= elem_to_remove.size())
    return true;

  auto el_it = elem_to_remove.begin();
  auto el_last = elem_to_remove.end();
  std::vector<Element>::iterator el_del;

  UInt deletions = 0;

  for (; el_it != el_last; ++el_it) {
    el_del = std::find(elem_list.begin(), elem_list.end(), *el_it);

    if (el_del != elem_list.end()) {
      elem_list.erase(el_del);
      ++deletions;
    }
  }

  AKANTU_DEBUG_ASSERT(deletions == 0 || deletions == elem_to_remove.size(),
                      "Not all elements have been erased");

  return deletions == 0;
}

/* -------------------------------------------------------------------------- */
void CohesiveElementInserterHelper::updateElementalConnectivity(
    Mesh & mesh, UInt old_node, UInt new_node,
    const std::vector<Element> & element_list,
    const std::vector<Element> * facet_list) {
  AKANTU_DEBUG_IN();

  for (auto & element : element_list) {
    if (element.type == _not_defined)
      continue;

    Vector<UInt> connectivity = mesh.getConnectivity(element);

    if (element.kind() == _ek_cohesive) {
      AKANTU_DEBUG_ASSERT(
          facet_list != nullptr,
          "Provide a facet list in order to update cohesive elements");

      const Vector<Element> facets =
          mesh_facets.getSubelementToElement(element);

      auto facet_nb_nodes = connectivity.size() / 2;

      /// loop over cohesive element's facets
      for (const auto & facet : enumerate(facets)) {
        /// skip facets if not present in the list
        if (std::find(facet_list->begin(), facet_list->end(),
                      std::get<1>(facet)) == facet_list->end()) {
          continue;
        }

        auto n = std::get<0>(facet);

        auto begin = connectivity.begin() + n * facet_nb_nodes;
        auto end = begin + facet_nb_nodes;

        auto it = std::find(begin, end, old_node);
        AKANTU_DEBUG_ASSERT(it != end, "Node not found in current element");

        *it = new_node;
      }
    } else {
      auto it = std::find(connectivity.begin(), connectivity.end(), old_node);
      AKANTU_DEBUG_ASSERT(it != connectivity.end(),
                          "Node not found in current element");

      /// update connectivity
      *it = new_node;
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void CohesiveElementInserterHelper::updateSubelementToElement(UInt dim,
                                                              bool facet_mode) {
  auto & facets_to_double = *facets_to_double_by_dim[dim];
  auto & facets_to_subfacets =
      elementsOfDimToElementsOfDim(dim + facet_mode, dim);

  for (auto && data :
       zip(make_view(facets_to_double, 2), facets_to_subfacets)) {
    const auto & old_subfacet = std::get<0>(data)[0];
    const auto & new_subfacet = std::get<0>(data)[1];
    auto & facet_to_subfacets = std::get<1>(data);

    MeshAccessor mesh_accessor(mesh_facets);
    for (auto & facet : facet_to_subfacets) {
      Vector<Element> subfacets = mesh_accessor.getSubelementToElement(facet);

      auto && subfacet_to_update_it =
          std::find(subfacets.begin(), subfacets.end(), old_subfacet);

      AKANTU_DEBUG_ASSERT(subfacet_to_update_it != subfacets.end(),
                          "Subfacet not found");

      *subfacet_to_update_it = new_subfacet;
    }
  }
}

/* -------------------------------------------------------------------------- */
void CohesiveElementInserterHelper::updateElementToSubelement(UInt dim,
                                                              bool facet_mode) {
  auto & facets_to_double = *facets_to_double_by_dim[dim];
  auto & facets_to_subfacets =
      elementsOfDimToElementsOfDim(dim + facet_mode, dim);

  MeshAccessor mesh_accessor(mesh_facets);
  // resize the arrays
  mesh_accessor.getElementToSubelement().initialize(
      mesh_facets, _spatial_dimension = dim, _with_nb_element = true);

  for (auto && data :
       zip(make_view(facets_to_double, 2), facets_to_subfacets)) {
    const auto & new_facet = std::get<0>(data)[1];
    mesh_accessor.getElementToSubelement(new_facet) = std::get<1>(data);
  }
}

/* -------------------------------------------------------------------------- */
void CohesiveElementInserterHelper::updateQuadraticSegments(UInt dim) {
  AKANTU_DEBUG_IN();
  auto spatial_dimension = mesh.getSpatialDimension();
  auto & facets_to_double = *facets_to_double_by_dim[dim];

  MeshAccessor mesh_accessor(mesh_facets);
  auto & connectivities = mesh_accessor.getConnectivities();

  /// this ones matter only for segments in 3D
  Array<std::vector<Element>> * element_to_subfacet_double = nullptr;
  Array<std::vector<Element>> * facet_to_subfacet_double = nullptr;

  if (dim == spatial_dimension - 2) {
    element_to_subfacet_double = &elementsOfDimToElementsOfDim(dim + 2, dim);
    facet_to_subfacet_double = &elementsOfDimToElementsOfDim(dim + 1, dim);
  }

  auto & element_to_subelement = mesh_facets.getElementToSubelement();
  std::vector<UInt> middle_nodes;

  for (auto && facet_to_double : make_view(facets_to_double, 2)) {
    auto & old_facet = facet_to_double[0];
    if (old_facet.type != _segment_3)
      continue;

    auto node = connectivities(old_facet, 2);
    if (not mesh.isPureGhostNode(node))
      middle_nodes.push_back(node);
  }

  auto n = doubled_nodes.size();
  doubleNodes(middle_nodes);

  for (auto && data : enumerate(make_view(facets_to_double, 2))) {
    auto facet = std::get<0>(data);
    auto & old_facet = std::get<1>(data)[0];
    if (old_facet.type != _segment_3)
      continue;

    auto old_node = connectivities(old_facet, 2);

    if (mesh.isPureGhostNode(old_node))
      continue;

    auto new_node = doubled_nodes(n, 1);

    auto & new_facet = std::get<1>(data)[1];
    connectivities(new_facet, 2) = new_node;

    if (dim == spatial_dimension - 2) {
      updateElementalConnectivity(mesh_facets, old_node, new_node,
                                  element_to_subelement(new_facet, 0));

      updateElementalConnectivity(mesh, old_node, new_node,
                                  (*element_to_subfacet_double)(facet),
                                  &(*facet_to_subfacet_double)(facet));
    } else {
      updateElementalConnectivity(mesh, old_node, new_node,
                                  element_to_subelement(new_facet, 0));
    }

    ++n;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
UInt CohesiveElementInserterHelper::insertCohesiveElement() {
  auto nb_new_facets = insertFacetsOnly();
  if (nb_new_facets == 0)
    return 0;

  updateCohesiveData();
  return nb_new_facets;
}

/* -------------------------------------------------------------------------- */
UInt CohesiveElementInserterHelper::insertFacetsOnly() {
  UInt spatial_dimension = mesh.getSpatialDimension();

  if (facets_to_double_by_dim[spatial_dimension - 1]->size() == 0)
    return 0;

  if (spatial_dimension == 1) {
    doublePointFacet();
  } else if (spatial_dimension == 2) {
    doubleFacets<1>();
    findSubfacetToDouble<1>();

    doubleSubfacet<2>();

  } else if (spatial_dimension == 3) {
    doubleFacets<2>();
    findSubfacetToDouble<2>();

    doubleFacets<1>();
    findSubfacetToDouble<1>();

    doubleSubfacet<3>();
  }

  return facets_to_double_by_dim[spatial_dimension - 1]->size();
}

/* -------------------------------------------------------------------------- */
template <UInt dim> void CohesiveElementInserterHelper::doubleFacets() {
  AKANTU_DEBUG_IN();
  NewElementsEvent new_facets;
  auto spatial_dimension = mesh_facets.getSpatialDimension();

  auto & facets_to_double = *facets_to_double_by_dim[dim];
  MeshAccessor mesh_accessor(mesh_facets);
  for (auto && facet_to_double : make_view(facets_to_double, 2)) {
    auto && old_facet = facet_to_double[0];
    auto && new_facet = facet_to_double[1];

    auto & facets_connectivities =
        mesh_accessor.getConnectivity(old_facet.type, old_facet.ghost_type);
    facets_connectivities.resize(
        nb_new_facets(old_facet.type, old_facet.ghost_type));

    auto facets_connectivities_begin =
        make_view(facets_connectivities, facets_connectivities.getNbComponent())
            .begin();

    // copy the connectivities
    Vector<UInt> new_conn(facets_connectivities_begin[new_facet.element]);
    Vector<UInt> old_conn(facets_connectivities_begin[old_facet.element]);
    new_conn = old_conn;

    // this will fail if multiple facet types exists for a given element type
    // \TODO handle multiple sub-facet types
    auto nb_subfacet_per_facet = Mesh::getNbFacetsPerElement(old_facet.type);

    auto & subfacets_to_facets = mesh_accessor.getSubelementToElementNC(
        old_facet.type, old_facet.ghost_type);

    subfacets_to_facets.resize(
        nb_new_facets(old_facet.type, old_facet.ghost_type), ElementNull);

    auto subfacets_to_facets_begin =
        make_view(subfacets_to_facets, nb_subfacet_per_facet).begin();

    // copy the subfacet to facets information
    Vector<Element> old_subfacets_to_facet(
        subfacets_to_facets_begin[old_facet.element]);
    Vector<Element> new_subfacet_to_facet(
        subfacets_to_facets_begin[new_facet.element]);
    new_subfacet_to_facet = old_subfacets_to_facet;

    for (auto & subfacet : old_subfacets_to_facet) {
      if (subfacet == ElementNull)
        continue;

      /// update facet_to_subfacet array
      mesh_accessor.getElementToSubelement(subfacet).push_back(new_facet);
    }

    new_facets.getList().push_back(new_facet);
  }

  /// update facet_to_subfacet and _segment_3 facets if any
  if (not(dim == spatial_dimension - 1)) {
    updateSubelementToElement(dim - 1, true);
    updateElementToSubelement(dim - 1, true);
    updateQuadraticSegments(dim);
  } else {
    updateQuadraticSegments(dim);
  }

  mesh_facets.sendEvent(new_facets);
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <UInt dim> void CohesiveElementInserterHelper::findSubfacetToDouble() {
  AKANTU_DEBUG_IN();

  UInt spatial_dimension = mesh_facets.getSpatialDimension();

  MeshAccessor mesh_accessor(mesh_facets);

  auto & facets_to_double = *facets_to_double_by_dim[spatial_dimension - 1];
  auto & subfacets_to_facets = mesh_facets.getSubelementToElement();
  auto & elements_to_facets = mesh_accessor.getElementToSubelement();

  /// loop on every facet
  for (auto f : arange(2)) {
    for (auto && facet_to_double : make_view(facets_to_double, 2)) {
      auto old_facet = facet_to_double[f];
      auto new_facet = facet_to_double[1 - f];

      const auto & starting_element = elements_to_facets(new_facet, 0)[0];
      auto current_facet = old_facet;

      Vector<Element> subfacets_to_facet = subfacets_to_facets.get(old_facet);
      /// loop on every subfacet
      for (auto & subfacet : subfacets_to_facet) {
        if (subfacet == ElementNull)
          continue;

        if (dim == spatial_dimension - 2) {
          Vector<Element> subsubfacets_to_subfacet =
              subfacets_to_facets.get(subfacet);
          /// loop on every subsubfacet
          for (auto & subsubfacet : subsubfacets_to_subfacet) {
            if (subsubfacet == ElementNull)
              continue;

            Vector<UInt> subsubfacet_connectivity =
                mesh_facets.getConnectivity(subsubfacet);

            std::vector<Element> element_list;
            std::vector<Element> facet_list;
            std::vector<Element> subfacet_list;

            /// check if subsubfacet is to be doubled
            if (findElementsAroundSubfacet(
                    starting_element, current_facet, subsubfacet_connectivity,
                    element_list, facet_list, &subfacet_list) == false &&
                removeElementsInVector(
                    subfacet_list, elements_to_facets(subsubfacet)) == false) {
              Element new_subsubfacet{
                  subsubfacet.type,
                  nb_new_facets(subsubfacet.type, subsubfacet.ghost_type)++,
                  subsubfacet.ghost_type};
              facets_to_double_by_dim[dim - 1]->push_back(
                  Vector<Element>{subsubfacet, new_subsubfacet});
              elementsOfDimToElementsOfDim(dim, dim - 1)
                  .push_back(subfacet_list);
              elementsOfDimToElementsOfDim(dim + 1, dim - 1)
                  .push_back(facet_list);
              elementsOfDimToElementsOfDim(dim + 2, dim - 1)
                  .push_back(element_list);
            }
          }
        } else {
          std::vector<Element> element_list;
          std::vector<Element> facet_list;

          Vector<UInt> subfacet_connectivity =
              mesh_facets.getConnectivity(subfacet);

          /// check if subfacet is to be doubled
          if (findElementsAroundSubfacet(starting_element, current_facet,
                                         subfacet_connectivity, element_list,
                                         facet_list) == false and
              removeElementsInVector(facet_list,
                                     elements_to_facets(subfacet)) == false) {
            Element new_subfacet{
                subfacet.type,
                nb_new_facets(subfacet.type, subfacet.ghost_type)++,
                subfacet.ghost_type};
            facets_to_double_by_dim[dim - 1]->push_back(
                Vector<Element>{subfacet, new_subfacet});

            elementsOfDimToElementsOfDim(dim, dim - 1).push_back(facet_list);
            elementsOfDimToElementsOfDim(dim + 1, dim - 1)
                .push_back(element_list);
          }
        }
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
void CohesiveElementInserterHelper::doubleNodes(
    const std::vector<UInt> & old_nodes) {
  auto & position = mesh.getNodes();
  auto spatial_dimension = mesh.getSpatialDimension();

  auto old_nb_nodes = position.size();
  position.reserve(old_nb_nodes + old_nodes.size());
  doubled_nodes.reserve(doubled_nodes.size() + old_nodes.size());

  auto position_begin = position.begin(spatial_dimension);
  for (auto && data : enumerate(old_nodes)) {
    auto n = std::get<0>(data);
    auto old_node = std::get<1>(data);
    decltype(old_node) new_node = old_nb_nodes + n;

    /// store doubled nodes
    doubled_nodes.push_back(Vector<UInt>{old_node, new_node});

    /// update position
    Vector<Real> coords = Vector<Real>(position_begin[old_node]);
    position.push_back(coords);
  }
}

/* -------------------------------------------------------------------------- */
bool CohesiveElementInserterHelper::findElementsAroundSubfacet(
    const Element & starting_element, const Element & end_facet,
    const Vector<UInt> & subfacet_connectivity,
    std::vector<Element> & element_list, std::vector<Element> & facet_list,
    std::vector<Element> * subfacet_list) {
  bool facet_matched = false;

  element_list.push_back(starting_element);

  std::queue<Element> elements_to_check;
  elements_to_check.push(starting_element);

  /// keep going as long as there are elements to check
  while (not elements_to_check.empty()) {
    /// check current element
    Element & current_element = elements_to_check.front();

    const Vector<Element> facets_to_element =
        mesh_facets.getSubelementToElement(current_element);

    // for every facet of the element
    for (auto & current_facet : facets_to_element) {
      if (current_facet == ElementNull)
        continue;

      if (current_facet == end_facet)
        facet_matched = true;

      auto find_facet =
          std::find(facet_list.begin(), facet_list.end(), current_facet);
      // facet already listed
      if (find_facet != facet_list.end())
        continue;

      // subfacet_connectivity is not in the connectivity of current_facet;
      if ((find_facet != facet_list.end()) or
          not hasElement(mesh_facets.getConnectivity(current_facet),
                         subfacet_connectivity))
        continue;

      facet_list.push_back(current_facet);

      if (subfacet_list != nullptr) {
        const Vector<Element> subfacets_of_facet =
            mesh_facets.getSubelementToElement(current_facet);

        /// check subfacets
        for (const auto & current_subfacet : subfacets_of_facet) {
          if (current_subfacet == ElementNull)
            continue;

          if ((std::find(subfacet_list->begin(), subfacet_list->end(),
                         current_subfacet) == subfacet_list->end()) and
              hasElement(mesh_facets.getConnectivity(current_subfacet),
                         subfacet_connectivity))
            subfacet_list->push_back(current_subfacet);
        }
      }

      /// consider opposing element
      const auto & elements_to_facet =
          mesh_facets.getElementToSubelement(current_facet);
      UInt opposing = 0;
      if (elements_to_facet[0] == current_element)
        opposing = 1;

      auto & opposing_element = elements_to_facet[opposing];

      /// skip null elements since they are on a boundary
      if (opposing_element == ElementNull)
        continue;

      /// skip this element if already added
      if (std::find(element_list.begin(), element_list.end(),
                    opposing_element) != element_list.end())
        continue;

      /// only regular elements have to be checked
      if (opposing_element.kind() == _ek_regular)
        elements_to_check.push(opposing_element);

      element_list.push_back(opposing_element);

      AKANTU_DEBUG_ASSERT(hasElement(mesh.getConnectivity(opposing_element),
                                     subfacet_connectivity),
                          "Subfacet doesn't belong to this element");
    }

    /// erased checked element from the list
    elements_to_check.pop();
  }

  return facet_matched;
}

/* -------------------------------------------------------------------------- */
void CohesiveElementInserterHelper::updateCohesiveData() {

  UInt spatial_dimension = mesh.getSpatialDimension();
  bool third_dimension = spatial_dimension == 3;

  MeshAccessor mesh_accessor(mesh);
  MeshAccessor mesh_facets_accessor(mesh_facets);

  for (auto ghost_type : ghost_types) {
    for (auto cohesive_type : mesh.elementTypes(_element_kind = _ek_cohesive)) {
      auto nb_cohesive_elements = mesh.getNbElement(cohesive_type, ghost_type);
      nb_new_facets(cohesive_type, ghost_type) = nb_cohesive_elements;

      mesh_facets_accessor.getSubelementToElement(cohesive_type, ghost_type);
    }
  }

  auto & facets_to_double = *facets_to_double_by_dim[spatial_dimension - 1];

  new_elements.reserve(new_elements.size() + facets_to_double.size());

  std::array<Element, 2> facets;

  auto & element_to_facet = mesh_facets_accessor.getElementToSubelement();
  auto & facets_to_cohesive_elements =
      mesh_facets_accessor.getSubelementToElement();

  for (auto && facet_to_double : make_view(facets_to_double, 2)) {
    auto & old_facet = facet_to_double[0];

    /// (in 3D cohesive elements connectivity is inverted)
    facets[third_dimension ? 1 : 0] = old_facet;
    facets[third_dimension ? 0 : 1] = facet_to_double[1];

    auto type_cohesive = FEEngine::getCohesiveElementType(old_facet.type);

    auto & facet_connectivities =
        mesh_facets.getConnectivity(old_facet.type, old_facet.ghost_type);
    auto facet_connectivity_it =
        facet_connectivities.begin(facet_connectivities.getNbComponent());

    auto cohesive_element = Element{
        type_cohesive, nb_new_facets(type_cohesive, old_facet.ghost_type)++,
        old_facet.ghost_type};

    auto & cohesives_connectivities =
        mesh_accessor.getConnectivity(type_cohesive, old_facet.ghost_type);
    Matrix<UInt> connectivity(facet_connectivities.getNbComponent(), 2);
    Vector<Element> facets_to_cohesive_element(2);

    for (auto s : arange(2)) {
      /// store doubled facets
      facets_to_cohesive_element[s] = facets[s];

      // update connectivities
      connectivity(s) = Vector<UInt>(facet_connectivity_it[facets[s].element]);

      /// update element_to_facet vectors
      element_to_facet(facets[s], 0)[1] = cohesive_element;
    }

    cohesives_connectivities.push_back(connectivity);

    facets_to_cohesive_elements(type_cohesive, old_facet.ghost_type)
        .push_back(facets_to_cohesive_element);

    /// add cohesive element to the element event list
    new_elements.push_back(cohesive_element);
  }
}

/* -------------------------------------------------------------------------- */
void CohesiveElementInserterHelper::doublePointFacet() {
  UInt spatial_dimension = mesh.getSpatialDimension();
  if (spatial_dimension != 1)
    return;

  auto & facets_to_double = *facets_to_double_by_dim[spatial_dimension - 1];
  auto & element_to_facet = mesh_facets.getElementToSubelement();
  auto & position = mesh.getNodes();
  MeshAccessor mesh_accessor(mesh_facets);

  for (auto ghost_type : ghost_types) {
    for (auto facet_type : nb_new_facets.elementTypes(
             spatial_dimension - 1, ghost_type, _ek_regular)) {
      auto nb_new_element = nb_new_facets(facet_type, ghost_type);
      auto & connectivities =
          mesh_accessor.getConnectivity(facet_type, ghost_type);
      connectivities.resize(connectivities.size() + nb_new_element);
    }
  }

  for (auto facet_to_double : make_view(facets_to_double, 2)) {
    auto & old_facet = facet_to_double[0];
    auto & new_facet = facet_to_double[1];

    auto element = element_to_facet(new_facet)[0];
    auto & facet_connectivities =
        mesh_accessor.getConnectivity(old_facet.type, old_facet.ghost_type);
    auto old_node = facet_connectivities(old_facet.element);
    auto new_node = position.size();

    /// update position
    position.push_back(position(old_node));
    facet_connectivities(new_facet.element) = new_node;

    Vector<UInt> segment_conectivity = mesh.getConnectivity(element);

    /// update facet connectivity
    auto it = std::find(segment_conectivity.begin(), segment_conectivity.end(),
                        old_node);
    *it = new_node;

    doubled_nodes.push_back(Vector<UInt>{old_node, new_node});
  }
}

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
void CohesiveElementInserterHelper::doubleSubfacet() {
  if (spatial_dimension == 1)
    return;

  std::vector<UInt> nodes_to_double;
  MeshAccessor mesh_accessor(mesh_facets);
  auto & connectivities = mesh_accessor.getConnectivities();

  auto & facets_to_double = *facets_to_double_by_dim[0];
  ElementTypeMap<Int> nb_new_elements;
  for (auto && facet_to_double : make_view(facets_to_double, 2)) {
    auto && old_element = facet_to_double[0];
    nb_new_elements(old_element.type, old_element.ghost_type) = 0;
  }

  for (auto && facet_to_double : make_view(facets_to_double, 2)) {
    auto && old_element = facet_to_double[0];
    ++nb_new_elements(old_element.type, old_element.ghost_type);
  }

  for (auto ghost_type : ghost_types) {
    for (auto facet_type :
         nb_new_elements.elementTypes(0, ghost_type, _ek_regular)) {
      auto & connectivities =
          mesh_accessor.getConnectivity(facet_type, ghost_type);
      connectivities.resize(connectivities.size() +
                            nb_new_elements(facet_type, ghost_type));
    }
  }

  for (auto && facet_to_double : make_view(facets_to_double, 2)) {
    auto & old_facet = facet_to_double(0);
    // auto & new_facet = facet_to_double(1);

    AKANTU_DEBUG_ASSERT(
        old_facet.type == _point_1,
        "Only _point_1 subfacet doubling is supported at the moment");

    /// list nodes double
    nodes_to_double.push_back(connectivities(old_facet));
  }

  auto old_nb_doubled_nodes = doubled_nodes.size();
  doubleNodes(nodes_to_double);

  auto double_nodes_view = make_view(doubled_nodes, 2);

  for (auto && data :
       zip(make_view(facets_to_double, 2),
           range(double_nodes_view.begin() + old_nb_doubled_nodes,
                 double_nodes_view.end()),
           arange(facets_to_double.size()))) {
    // auto & old_facet = std::get<0>(data)[0];
    auto & new_facet = std::get<0>(data)[1];

    auto & nodes = std::get<1>(data);
    auto old_node = nodes(0);
    auto new_node = nodes(1);

    auto f = std::get<2>(data);

    /// add new nodes in connectivity
    connectivities(new_facet) = new_node;

    updateElementalConnectivity(mesh, old_node, new_node,
                                elementsOfDimToElementsOfDim(2, 0)(f),
                                &elementsOfDimToElementsOfDim(1, 0)(f));

    updateElementalConnectivity(mesh_facets, old_node, new_node,
                                elementsOfDimToElementsOfDim(1, 0)(f));

    if (spatial_dimension == 3)
      updateElementalConnectivity(mesh_facets, old_node, new_node,
                                  elementsOfDimToElementsOfDim(0, 0)(f));
  }

  updateSubelementToElement(0, spatial_dimension == 2);
  updateElementToSubelement(0, spatial_dimension == 2);
}

} // namespace akantu

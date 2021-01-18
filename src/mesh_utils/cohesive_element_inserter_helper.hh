/**
 * @file   cohesive_element_inserter_helper.hh
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
#include "aka_array.hh"
#include "mesh.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_COHESIVE_ELEMENT_INSERTER_HELPER_HH__
#define __AKANTU_COHESIVE_ELEMENT_INSERTER_HELPER_HH__

namespace akantu {

/* -------------------------------------------------------------------------- */
class CohesiveElementInserterHelper {
public:
  CohesiveElementInserterHelper(
      Mesh & mesh, const ElementTypeMapArray<bool> & facet_insertion);

  UInt insertCohesiveElement();
  UInt insertFacetsOnly();

private:
  template <UInt dim> UInt insertFacetsOnlyImpl();
  template <UInt dim> void doubleFacets();
  template <UInt dim> void findSubfacetToDouble();

  void doubleNodes(const std::vector<UInt> & old_nodes);

  bool findElementsAroundSubfacet(
      const Element & starting_element, const Element & end_facet,
      const Vector<UInt> & subfacet_connectivity,
      std::vector<Element> & element_list, std::vector<Element> & facet_list,
      std::vector<Element> * subfacet_list = nullptr);

  static inline bool hasElement(const Vector<UInt> & nodes_element,
                                const Vector<UInt> & nodes);
  static inline bool
  removeElementsInVector(const std::vector<Element> & elem_to_remove,
                         std::vector<Element> & elem_list);

  void updateElementalConnectivity(
      Mesh & mesh, UInt old_node, UInt new_node,
      const std::vector<Element> & element_list,
      const std::vector<Element> * facet_list = nullptr);

  // update functions
  void updateElementToSubelement(UInt dim, bool facet_mode);
  void updateSubelementToElement(UInt dim, bool facet_mode);
  void updateQuadraticSegments(UInt dim);

  void updateCohesiveData();
  void doublePointFacet();
  template <UInt spatial_dimension> void doubleSubfacet();

  decltype(auto) elementsOfDimToElementsOfDim(Int dim1, Int dim2) {
    AKANTU_DEBUG_ASSERT(dim1 >= 0 and dim1 <= 3,
                        "dimension of target element out of range");
    AKANTU_DEBUG_ASSERT(dim2 >= 0 and dim2 <= 3,
                        "dimension of source element out of range");

    auto & array = dimelements_to_dimelements[dim1][dim2];
    if (not array) {
      array = std::make_unique<Array<std::vector<Element>>>();
    }

    return (*array);
  }

public:
  decltype(auto) getNewElements() const { return (new_elements); }
  decltype(auto) getDoubledNodes() const { return (doubled_nodes); }

private:
  std::array<std::unique_ptr<Array<Element>>, 3> facets_to_double_by_dim;
  std::array<std::array<std::unique_ptr<Array<std::vector<Element>>>, 2>, 4>
      dimelements_to_dimelements;

  Array<UInt> doubled_nodes;
  Array<Element> new_elements;

  Mesh & mesh;
  Mesh & mesh_facets;

  ElementTypeMap<UInt> nb_new_facets;
};

} // namespace akantu

#endif /* __AKANTU_COHESIVE_ELEMENT_INSERTER_HELPER_HH__ */

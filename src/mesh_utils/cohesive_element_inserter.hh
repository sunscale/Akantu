/**
 * @file   cohesive_element_inserter.hh
 *
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Wed Dec 04 2013
 * @date last modification: Sun Feb 04 2018
 *
 * @brief  Cohesive element inserter
 *
 * @section LICENSE
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
#include "data_accessor.hh"
#include "mesh_utils.hh"
#include "parsable.hh"
/* -------------------------------------------------------------------------- */
#include <numeric>
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_COHESIVE_ELEMENT_INSERTER_HH__
#define __AKANTU_COHESIVE_ELEMENT_INSERTER_HH__

namespace akantu {
class GlobalIdsUpdater;
class FacetSynchronizer;
} // namespace akantu

namespace akantu {

class CohesiveElementInserter : public DataAccessor<Element>, public Parsable {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  CohesiveElementInserter(Mesh & mesh,
                          const ID & id = "cohesive_element_inserter");

  ~CohesiveElementInserter() override;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// set range limitation for intrinsic cohesive element insertion
  void setLimit(SpatialDirection axis, Real first_limit, Real second_limit);

  /// insert intrinsic cohesive elements in a predefined range
  UInt insertIntrinsicElements();

  /// insert extrinsic cohesive elements (returns the number of new
  /// cohesive elements)
  UInt insertElements(bool only_double_facets = false);

  /// limit check facets to match given insertion limits
  void limitCheckFacets();

protected:
  void parseSection(const ParserSection & section) override;

protected:
  /// internal version of limitCheckFacets
  void limitCheckFacets(ElementTypeMapArray<bool> & check_facets);

  /// update facet insertion arrays after facets doubling
  void updateInsertionFacets();

  /// functions for parallel communications
  inline UInt getNbData(const Array<Element> & elements,
                        const SynchronizationTag & tag) const override;

  inline void packData(CommunicationBuffer & buffer,
                       const Array<Element> & elements,
                       const SynchronizationTag & tag) const override;

  inline void unpackData(CommunicationBuffer & buffer,
                         const Array<Element> & elements,
                         const SynchronizationTag & tag) override;
  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  AKANTU_GET_MACRO_NOT_CONST(InsertionFacetsByElement, insertion_facets,
                             ElementTypeMapArray<bool> &);
  AKANTU_GET_MACRO(InsertionFacetsByElement, insertion_facets,
                   const ElementTypeMapArray<bool> &);
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE(InsertionFacets, insertion_facets, bool);

  AKANTU_GET_MACRO(CheckFacets, check_facets,
                   const ElementTypeMapArray<bool> &);
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE(CheckFacets, check_facets, bool);
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(CheckFacets, check_facets, bool);

  AKANTU_GET_MACRO(MeshFacets, mesh_facets, const Mesh &);
  AKANTU_GET_MACRO_NOT_CONST(MeshFacets, mesh_facets, Mesh &);

  AKANTU_SET_MACRO(IsExtrinsic, is_extrinsic, bool);

public:
  friend class SolidMechanicsModelCohesive;

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// object id
  ID id;

  /// main mesh where to insert cohesive elements
  Mesh & mesh;

  /// mesh containing facets
  Mesh & mesh_facets;

  /// list of facets where to insert elements
  ElementTypeMapArray<bool> insertion_facets;

  /// limits for element insertion
  Matrix<Real> insertion_limits;

  /// list of groups to consider for insertion, ignored if empty
  std::set<ID> physical_surfaces;

  /// list of groups in between which an inside which element are insterted
  std::set<ID> physical_zones;

  /// vector containing facets in which extrinsic cohesive elements can be
  /// inserted
  ElementTypeMapArray<bool> check_facets;

  /// global connectivity ids updater
  std::unique_ptr<GlobalIdsUpdater> global_ids_updater;

  /// is this inserter used in extrinsic
  bool is_extrinsic{false};
};

class CohesiveNewNodesEvent : public NewNodesEvent {
public:
  CohesiveNewNodesEvent(const std::string & origin) : NewNodesEvent(origin) {}
  ~CohesiveNewNodesEvent() override = default;

  AKANTU_GET_MACRO_NOT_CONST(OldNodesList, old_nodes, Array<UInt> &);
  AKANTU_GET_MACRO(OldNodesList, old_nodes, const Array<UInt> &);

private:
  Array<UInt> old_nodes;
};

} // namespace akantu

#include "cohesive_element_inserter_inline_impl.cc"

#endif /* __AKANTU_COHESIVE_ELEMENT_INSERTER_HH__ */

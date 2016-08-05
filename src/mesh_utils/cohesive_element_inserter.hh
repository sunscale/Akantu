/**
 * @file   cohesive_element_inserter.hh
 *
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Wed Dec 04 2013
 * @date last modification: Fri Oct 02 2015
 *
 * @brief  Cohesive element inserter
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
#ifndef __AKANTU_COHESIVE_ELEMENT_INSERTER_HH__
#define __AKANTU_COHESIVE_ELEMENT_INSERTER_HH__

/* -------------------------------------------------------------------------- */
#include <numeric>
#include "data_accessor.hh"
#include "mesh_utils.hh"

#if defined(AKANTU_PARALLEL_COHESIVE_ELEMENT)
#  include "global_ids_updater.hh"
#  include "facet_synchronizer.hh"
#endif

/* -------------------------------------------------------------------------- */
__BEGIN_AKANTU__

class CohesiveElementInserter : public DataAccessor {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  CohesiveElementInserter(Mesh & mesh,
			  bool is_extrinsic = false,
			  DistributedSynchronizer * synchronizer = NULL,
			  const ID & id = "cohesive_element_inserter");

  virtual ~CohesiveElementInserter();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  /// init function
  void init(bool is_extrinsic);

  /// set range limitation for intrinsic cohesive element insertion
  void setLimit(SpacialDirection axis, Real first_limit, Real second_limit);

  /// insert intrinsic cohesive elements in a predefined range
  UInt insertIntrinsicElements();

  /// preset insertion of intrinsic cohesive elements along 
  /// a predefined group of facet and assign them a defined material index.
  /// insertElement() method has to be called to finalize insertion. 
  void insertIntrinsicElements(std::string physname, UInt material_index);

  /// insert extrinsic cohesive elements (returns the number of new
  /// cohesive elements)
  UInt insertElements(bool only_double_facets = false);

  /// function to print the contain of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;

  /// limit check facets to match given insertion limits
  void limitCheckFacets();

#if defined(AKANTU_PARALLEL_COHESIVE_ELEMENT)
  /// init parallel variables
  void initParallel(FacetSynchronizer * facet_synchronizer,
		    DistributedSynchronizer * distributed_synchronizer);
#endif

protected:

  /// init facets check
  void initFacetsCheck();

  /// update facet insertion arrays after facets doubling
  void updateInsertionFacets();

#if defined(AKANTU_PARALLEL_COHESIVE_ELEMENT)
  /// update nodes type and global ids for parallel simulations
  /// (locally, within each processor)
  UInt updateGlobalIDs(NewNodesEvent & node_event);

  /// synchronize the global ids among the processors in parallel simulations
  void synchronizeGlobalIDs(NewNodesEvent & node_event);

  /// update nodes type
  void updateNodesType(Mesh & mesh, NewNodesEvent & node_event);

  /// functions for parallel communications
  inline UInt getNbDataForElements(const Array<Element> & elements,
				   SynchronizationTag tag) const;

  inline void packElementData(CommunicationBuffer & buffer,
			      const Array<Element> & elements,
			      SynchronizationTag tag) const;

  inline void unpackElementData(CommunicationBuffer & buffer,
				const Array<Element> & elements,
				SynchronizationTag tag);

  template<bool pack_mode>
  inline void packUnpackGlobalConnectivity(CommunicationBuffer & buffer,
					   const Array<Element> & elements) const;

  template<bool pack_mode>
  inline void packUnpackGroupedInsertionData(CommunicationBuffer & buffer,
					     const Array<Element> & elements) const;
#endif

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  AKANTU_GET_MACRO_NOT_CONST(InsertionFacetsByElement,
			     insertion_facets,
			     ElementTypeMapArray<bool> &);

  AKANTU_GET_MACRO(InsertionFacetsByElement,
		   insertion_facets,
		   const ElementTypeMapArray<bool> &);

  AKANTU_GET_MACRO_BY_ELEMENT_TYPE(InsertionFacets, insertion_facets, bool);

  AKANTU_GET_MACRO_BY_ELEMENT_TYPE(CheckFacets, check_facets, bool);
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(CheckFacets, check_facets, bool);
  AKANTU_GET_MACRO(MeshFacets, mesh_facets, const Mesh &);
  AKANTU_GET_MACRO_NOT_CONST(MeshFacets, mesh_facets, Mesh &);

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
  Array<Real> insertion_limits;

  /// vector containing facets in which extrinsic cohesive elements can be inserted
  ElementTypeMapArray<bool> check_facets;

#if defined(AKANTU_PARALLEL_COHESIVE_ELEMENT)
  /// facet synchronizer
  FacetSynchronizer * facet_synchronizer;

  /// global connectivity ids updater
  GlobalIdsUpdater * global_ids_updater;
#endif
};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#if defined(AKANTU_PARALLEL_COHESIVE_ELEMENT)
#  include "cohesive_element_inserter_inline_impl.cc"
#endif

/// standard output stream operator
inline std::ostream & operator <<(std::ostream & stream, const CohesiveElementInserter & _this)
{
  _this.printself(stream);
  return stream;
}


__END_AKANTU__


#endif /* __AKANTU_COHESIVE_ELEMENT_INSERTER_HH__ */

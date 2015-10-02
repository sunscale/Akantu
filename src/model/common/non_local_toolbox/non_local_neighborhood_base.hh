/**
 * @file   non_local_neighborhood_base.hh
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @date   Mon Sep 21 15:43:26 2015
 *
 * @brief  Non-local neighborhood base class
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
#ifndef __AKANTU_NON_LOCAL_NEIGHBORHOOD_BASE_HH__
#define __AKANTU_NON_LOCAL_NEIGHBORHOOD_BASE_HH__
/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "solid_mechanics_model.hh"
#include "aka_grid_dynamic.hh"
#include "grid_synchronizer.hh"
#include "aka_memory.hh"
/* -------------------------------------------------------------------------- */



__BEGIN_AKANTU__

class NonLocalNeighborhoodBase : public Memory {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  NonLocalNeighborhoodBase(const SolidMechanicsModel & model, Real radius,
			   const ElementTypeMapReal & quad_coordinates,
			   const ID & id = "neighborhood",
			   const MemoryID & memory_id = 0);
  virtual ~NonLocalNeighborhoodBase();

  typedef std::vector< std::pair<QuadraturePoint, QuadraturePoint> > PairList;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */

public:
  /// initialize the material computed parameter
  inline void insertQuad(const QuadraturePoint & quad, const Vector<Real> & coords);

  /// create the pairs of quadrature points
  void updatePairList();

  /// save the pairs of quadrature points in a file
  void savePairs(const std::string & filename) const;

  /// save the coordinates of all neighbors of a quad
  void saveNeighborCoords(const std::string & filename) const;

  /// create grid synchronizer and exchange ghost cells
  void createGridSynchronizer();

  /// compute weights, for instance needed for non-local damage computation
  virtual void computeWeights() {};

  /// compute the non-local counter part for a given element type map
  void weightedAvergageOnNeighbours(const ElementTypeMapReal & to_accumulate,
				    ElementTypeMapReal & accumulated,
				    UInt nb_degree_of_freedom,
				    const GhostType & ghost_type2) const {};

protected:

  /// intialize the neighborhood
  void initNeighborhood();

  /// create the grid
  void createGrid();

  void cleanupExtraGhostElement(const ElementTypeMap<UInt> & nb_ghost_protected) {};

  // virtual inline UInt getNbDataForElements(const Array<Element> & elements,
  // 					 SynchronizationTag tag) const;

  // virtual inline void packElementData(CommunicationBuffer & buffer,
  // 				      const Array<Element> & elements,
  // 				      SynchronizationTag tag) const;

  // virtual inline void unpackElementData(CommunicationBuffer & buffer,
  // 					const Array<Element> & elements,
  // 					SynchronizationTag tag);

  // virtual inline void onElementsAdded(const Array<Element> & element_list);
  // virtual inline void onElementsRemoved(const Array<Element> & element_list,
  // 					const ElementTypeMapArray<UInt> & new_numbering,
  // 					const RemovedElementsEvent & event) {};

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  
  /// the model to which the neighborhood belongs
  const SolidMechanicsModel & model;

  /// Radius of non-local neighborhood
  Real non_local_radius;

  /**
   * the pairs of quadrature points
   * 0: not ghost to not ghost
   * 1: not ghost to ghost
   */
  PairList pair_list[2];

  /// the regular grid to construct/update the pair lists
  SpatialGrid<QuadraturePoint> * spatial_grid;

  bool is_creating_grid;

  /// the grid synchronizer for parallel computations
  GridSynchronizer * grid_synchronizer;

  /// the quadrature point positions
  const ElementTypeMapReal & quad_coordinates;
};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "non_local_neighborhood_base_inline_impl.cc"

__END_AKANTU__
#endif /* __AKANTU_NON_LOCAL_NEIGHBORHOOD_BASE_HH__ */

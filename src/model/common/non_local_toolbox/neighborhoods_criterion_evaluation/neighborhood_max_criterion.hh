/**
 * @file   neighborhood_max_criterion.hh
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 *
 * @date creation: Sat Sep 26 2015
 * @date last modification: Wed Jan 31 2018
 *
 * @brief  Neighborhood to find a maximum value in a neighborhood
 *
 * @section LICENSE
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
#ifndef __AKANTU_NEIGHBORHOOD_MAX_CRITERION_BASE_HH__
#define __AKANTU_NEIGHBORHOOD_MAX_CRITERION_BASE_HH__
/* -------------------------------------------------------------------------- */
#include "neighborhood_base.hh"
#include "parsable.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

class NeighborhoodMaxCriterion : public NeighborhoodBase, public Parsable {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  NeighborhoodMaxCriterion(Model & model,
                           const ElementTypeMapReal & quad_coordinates,
                           const ID & criterion_id,
                           const ID & id = "neighborhood_max_criterion",
                           const MemoryID & memory_id = 0);
  ~NeighborhoodMaxCriterion() override;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */

public:
  /// initialize the neighborhood
  void initNeighborhood() override;

  /// create grid synchronizer and exchange ghost cells
  void createGridSynchronizer() override;

  /// find the quads which have the maximum criterion in their neighborhood
  void findMaxQuads(std::vector<IntegrationPoint> & max_quads);

protected:
  /// remove unneccessary ghost elements
  void
  cleanupExtraGhostElements(const ElementTypeMap<UInt> & nb_ghost_protected);

  /// insert the quadrature points in the grid
  void insertAllQuads(const GhostType & ghost_type);

  /// compare criterion with neighbors
  void checkNeighbors(const GhostType & ghost_type);

  /* --------------------------------------------------------------------------
   */
  /* DataAccessor inherited members */
  /* --------------------------------------------------------------------------
   */
public:
  virtual inline UInt getNbDataForElements(const Array<Element> & elements,
                                           SynchronizationTag tag) const;

  virtual inline void packElementData(CommunicationBuffer & buffer,
                                      const Array<Element> & elements,
                                      SynchronizationTag tag) const;

  virtual inline void unpackElementData(CommunicationBuffer & buffer,
                                        const Array<Element> & elements,
                                        SynchronizationTag tag);

  /* --------------------------------------------------------------------------
   */
  /* Accessors */
  /* --------------------------------------------------------------------------
   */
public:
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// a boolean to store the information if a quad has the max
  /// criterion in the neighborhood
  ElementTypeMapArray<bool> is_highest;

  /// an element type map to store the flattened internal of the criterion
  ElementTypeMapReal criterion;
};

} // namespace akantu

#include "neighborhood_max_criterion_inline_impl.cc"

#endif /* __AKANTU_NEIGHBORHOOD_MAX_CRITERION_BASE_HH__ */

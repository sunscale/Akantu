/**
 * @file   non_local_neighborhood_base.hh
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Sat Sep 26 2015
 * @date last modification: Tue Feb 20 2018
 *
 * @brief  Non-local neighborhood base class
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
#include "neighborhood_base.hh"
#include "parsable.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_NON_LOCAL_NEIGHBORHOOD_BASE_HH_
#define AKANTU_NON_LOCAL_NEIGHBORHOOD_BASE_HH_

namespace akantu {
class Model;
}
/* -------------------------------------------------------------------------- */
namespace akantu {

class NonLocalNeighborhoodBase : public NeighborhoodBase, public Parsable {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  NonLocalNeighborhoodBase(Model & model,
                           const ElementTypeMapReal & quad_coordinates,
                           const ID & id = "non_local_neighborhood");
  ~NonLocalNeighborhoodBase() override;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// create grid synchronizer and exchange ghost cells
  void createGridSynchronizer() override;

  void synchronize(DataAccessor<Element> & data_accessor,
                   const SynchronizationTag & tag) override;

  /// compute weights, for instance needed for non-local damage computation
  virtual void computeWeights(){};

  // compute the non-local counter part for a given element type map
  virtual void
  weightedAverageOnNeighbours(const ElementTypeMapReal & to_accumulate,
                              ElementTypeMapReal & accumulated,
                              UInt nb_degree_of_freedom,
                              GhostType ghost_type2) const = 0;

  /// update the weights for the non-local averaging
  virtual void updateWeights() = 0;

  /// update the weights for the non-local averaging
  virtual void saveWeights(const std::string & /*unused*/) const {
    AKANTU_TO_IMPLEMENT();
  }

  /// register a new non-local variable in the neighborhood
  virtual void registerNonLocalVariable(const ID & id);

  /// clean up the unneccessary ghosts
  void cleanupExtraGhostElements(std::set<Element> & relevant_ghost_elements);

  /// list releveant ghosts
  void getRelevantGhostElements(std::set<Element> & relevant_ghost_elements);

protected:
  /// create the grid
  void createGrid();

  /* --------------------------------------------------------------------------
   */
  /* DataAccessor inherited members */
  /* --------------------------------------------------------------------------
   */
public:
  inline UInt getNbData(const Array<Element> & /*elements*/,
                        const SynchronizationTag & /*tag*/) const override {
    return 0;
  }

  inline void packData(CommunicationBuffer & /*buffer*/,
                       const Array<Element> & /*element*/,
                       const SynchronizationTag & /*tag*/) const override {}

  inline void unpackData(CommunicationBuffer & /*buffer*/,
                         const Array<Element> & /*element*/,
                         const SynchronizationTag & /*tag*/) override {}

  /* --------------------------------------------------------------------------
   */
  /* Accessors */
  /* --------------------------------------------------------------------------
   */
public:
  AKANTU_GET_MACRO(NonLocalVariables, non_local_variables,
                   const std::set<ID> &);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// list of non-local variables associated to the neighborhood
  std::set<ID> non_local_variables;
};

} // namespace akantu

#endif /* AKANTU_NON_LOCAL_NEIGHBORHOOD_BASE_HH_ */

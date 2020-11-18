/**
 * @file   non_local_manager.hh
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Wed Nov 08 2017
 *
 * @brief  Classes that manages all the non-local neighborhoods
 *
 *
 * Copyright (©)  2010-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "communication_buffer.hh"
#include "data_accessor.hh"
#include "mesh_events.hh"
#include "non_local_manager_callback.hh"
#include "parsable.hh"
/* -------------------------------------------------------------------------- */
#include <map>
#include <set>
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
#ifndef AKANTU_NON_LOCAL_MANAGER_HH_
#define AKANTU_NON_LOCAL_MANAGER_HH_

namespace akantu {
class Model;
class NonLocalNeighborhoodBase;
class GridSynchronizer;
class SynchronizerRegistry;
class IntegrationPoint;
template <typename T> class SpatialGrid;
class FEEngine;
} // namespace akantu

namespace akantu {

class NonLocalManager : protected Memory,
                        public MeshEventHandler,
                        public Parsable {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  NonLocalManager(Model & model, NonLocalManagerCallback & callback,
                  const ID & id = "non_local_manager",
                  const MemoryID & memory_id = 0);
  ~NonLocalManager() override;
  using NeighborhoodMap =
      std::map<ID, std::unique_ptr<NonLocalNeighborhoodBase>>;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* -----------------------------------------------------------------------  */
public:
  /// register a new internal needed for the weight computations
  ElementTypeMapReal & registerWeightFunctionInternal(const ID & field_name);

  /// register a non-local variable
  void registerNonLocalVariable(const ID & variable_name,
                                const ID & nl_variable_name, UInt nb_component);

  /// register non-local neighborhood
  inline void registerNeighborhood(const ID & neighborhood,
                                   const ID & weight_func_id);

  // void registerNonLocalManagerCallback(NonLocalManagerCallback & callback);

  /// average the internals and compute the non-local stresses
  virtual void computeAllNonLocalStresses();

  /// initialize the non-local manager: compute pair lists and weights for all
  /// neighborhoods
  virtual void initialize();

  /// synchronize once on a given tag using the neighborhoods synchronizer
  void synchronize(DataAccessor<Element> & data_accessor,
                   const SynchronizationTag & /*tag*/);

protected:
  /// create the grid synchronizers for each neighborhood
  void createNeighborhoodSynchronizers();

  /// compute the weights in each neighborhood for non-local averaging
  void computeWeights();

  /// compute the weights in each neighborhood for non-local averaging
  void updatePairLists();

  /// average the non-local variables
  void averageInternals(GhostType ghost_type = _not_ghost);

  /// update the flattened version of the weight function internals
  void updateWeightFunctionInternals();

protected:
  /// create a new neighborhood for a given domain ID
  void createNeighborhood(const ID & weight_func, const ID & neighborhood);

  /// set the values of the jacobians
  void setJacobians(const FEEngine & fe_engine, ElementKind kind);

  /// allocation of eelment type maps
  // void initElementTypeMap(UInt nb_component,
  //                         ElementTypeMapReal & element_map,
  //                         const FEEngine & fe_engine,
  //                         const ElementKind el_kind = _ek_regular);

  /// resizing of element type maps
  void resizeElementTypeMap(UInt nb_component, ElementTypeMapReal & element_map,
                            const FEEngine & fee,
                             ElementKind el_kind = _ek_regular);

  /// remove integration points from element type maps
  static void removeIntegrationPointsFromMap(
      const ElementTypeMapArray<UInt> & new_numbering, UInt nb_component,
      ElementTypeMapReal & element_map, const FEEngine & fee,
      ElementKind el_kind = _ek_regular);

  /// allocate the non-local variables
  void initNonLocalVariables();

  /// cleanup unneccessary ghosts
  void
  cleanupExtraGhostElements(); // ElementTypeMap<UInt> & nb_ghost_protected);

  /* ------------------------------------------------------------------------ */
  /* DataAccessor kind of interface                                           */
  /* ------------------------------------------------------------------------ */
public:
  /// get Nb data for synchronization in parallel
  UInt getNbData(const Array<Element> & elements, const ID & id) const;

  /// pack data for synchronization in parallel
  void packData(CommunicationBuffer & buffer, const Array<Element> & elements,
                const ID & id) const;

  /// unpack data for synchronization in parallel
  void unpackData(CommunicationBuffer & buffer, const Array<Element> & elements,
                  const ID & id) const;

  /* ------------------------------------------------------------------------ */
  /* MeshEventHandler inherited members                                       */
  /* ------------------------------------------------------------------------ */
public:
  void onElementsRemoved(const Array<Element> & element_list,
                         const ElementTypeMapArray<UInt> & new_numbering,
                         const RemovedElementsEvent & event) override;
  void onElementsAdded(const Array<Element> & element_list,
                       const NewElementsEvent & event) override;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  AKANTU_GET_MACRO(SpatialDimension, spatial_dimension, UInt);
  AKANTU_GET_MACRO(Model, model, const Model &);
  AKANTU_GET_MACRO_NOT_CONST(Model, model, Model &);
  AKANTU_GET_MACRO_NOT_CONST(Volumes, volumes, ElementTypeMapReal &)
  AKANTU_GET_MACRO(NbStressCalls, compute_stress_calls, UInt);

  /// return the fem object associated with a provided name
  inline NonLocalNeighborhoodBase & getNeighborhood(const ID & name) const;

  inline const Array<Real> & getJacobians(ElementType type,
                                          GhostType ghost_type) {
    return *jacobians(type, ghost_type);
  }

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// the spatial dimension
  const UInt spatial_dimension;

protected:
  /// the non-local neighborhoods present
  NeighborhoodMap neighborhoods;

  /// list of all the non-local materials in the model
  //  std::vector<ID> non_local_materials;

  struct NonLocalVariable {
    NonLocalVariable(const ID & variable_name, const ID & nl_variable_name,
                     const ID & id, UInt nb_component)
        : local(variable_name, id, 0), non_local(nl_variable_name, id, 0),
          nb_component(nb_component) {}
    ElementTypeMapReal local;
    ElementTypeMapReal non_local;
    UInt nb_component;
  };

  /// the non-local variables associated to a certain neighborhood
  std::map<ID, std::unique_ptr<NonLocalVariable>> non_local_variables;

  /// reference to the model
  Model & model;

  /// jacobians for all the elements in the mesh
  ElementTypeMap<const Array<Real> *> jacobians;

  /// store the position of the quadrature points
  ElementTypeMapReal integration_points_positions;

  /// store the volume of each quadrature point for the non-local weight
  /// normalization
  ElementTypeMapReal volumes;

  /// counter for computeStress calls
  UInt compute_stress_calls;

  /// map to store weight function types from input file
  std::map<ID, ParserSection> weight_function_types;

  /// map to store the internals needed by the weight functions
  std::map<ID, std::unique_ptr<ElementTypeMapReal>> weight_function_internals;
  /* --------------------------------------------------------------------------
   */
  /// the following are members needed to make this processor participate in the
  /// grid creation of neighborhoods he doesn't own as a member. For details see
  /// createGridSynchronizers function

  /// synchronizer registry for dummy grid synchronizers
  std::unique_ptr<SynchronizerRegistry> dummy_registry;

  /// map of dummy synchronizers
  std::map<ID, std::unique_ptr<GridSynchronizer>> dummy_synchronizers;

  /// dummy spatial grid
  std::unique_ptr<SpatialGrid<IntegrationPoint>> dummy_grid;

  /// create a set of all neighborhoods present in the simulation
  std::set<ID> global_neighborhoods;

  class DummyDataAccessor : public DataAccessor<Element> {
  public:
    inline UInt getNbData(const Array<Element> & /*elements*/,
                          const SynchronizationTag & /*tag*/) const override {
      return 0;
    };

    inline void packData(CommunicationBuffer & /*buffer*/,
                         const Array<Element> & /*element*/,
                         const SynchronizationTag & /*tag*/) const override{};

    inline void unpackData(CommunicationBuffer & /*buffer*/,
                           const Array<Element> & /*element*/,
                           const SynchronizationTag & /*tag*/) override{};
  };

  DummyDataAccessor dummy_accessor;

  /* ------------------------------------------------------------------------ */
  NonLocalManagerCallback * callback;
};

} // namespace akantu

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "non_local_manager_inline_impl.hh"

#endif /* AKANTU_NON_LOCAL_MANAGER_HH_ */

/**
 * @file   non_local_manager.hh
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Tue Dec 08 2015
 *
 * @brief  Classes that manages all the non-local neighborhoods
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2012, 2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de
 * Lausanne)  Laboratory (LSMS  -  Laboratoire de  Simulation  en Mécanique  des
 * Solides)
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
#ifndef __AKANTU_NON_LOCAL_MANAGER_HH__
#define __AKANTU_NON_LOCAL_MANAGER_HH__
/* -------------------------------------------------------------------------- */
#include "aka_memory.hh"
#include "solid_mechanics_model.hh"
#include "non_local_neighborhood_base.hh"
#include "mesh_events.hh"
#include "parsable.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

class NonLocalManager : public Memory,
                        public Parsable,
                        public MeshEventHandler {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  NonLocalManager(SolidMechanicsModel & model,
                  const ID & id,
                  const MemoryID & memory_id);
  virtual ~NonLocalManager();
  typedef std::map<ID, NonLocalNeighborhoodBase *> NeighborhoodMap;
  typedef std::pair<ID, ID> KeyCOO;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* -----------------------------------------------------------------------  */
public:
  /// initialize the non-local manager: compute pair lists and weights for all
  /// neighborhoods
  virtual void init();

  /// insert new quadrature point in the grid
  inline void insertQuad(const IntegrationPoint & quad,
                         const Vector<Real> & coords, const ID & neighborhood);

  /// register non-local neighborhood
  inline void registerNeighborhood(const ID & neighborhood,
                                   const ID & weight_func_id);

  /// associate a non-local variable to a neighborhood
  void nonLocalVariableToNeighborhood(const ID & id, const ID & neighborhood);

  /// return the fem object associated with a provided name
  inline NonLocalNeighborhoodBase & getNeighborhood(const ID & name) const;

  /// create the grid synchronizers for each neighborhood
  void createNeighborhoodSynchronizers();

  /// compute the weights in each neighborhood for non-local averaging
  inline void computeWeights();

  /// compute the weights in each neighborhood for non-local averaging
  inline void updatePairLists();

  /// register a new non-local material
  inline void registerNonLocalMaterial(Material & new_mat);

  /// register a non-local variable
  inline void registerNonLocalVariable(const ID & variable_name,
                                       const ID & nl_variable_name,
                                       UInt nb_component);

  /// average the non-local variables
  void averageInternals(const GhostType & ghost_type = _not_ghost);

  /// average the internals and compute the non-local stresses
  virtual void computeAllNonLocalStresses();

  /// register a new internal needed for the weight computations
  inline ElementTypeMapReal &
  registerWeightFunctionInternal(const ID & field_name);

  /// update the flattened version of the weight function internals
  inline void updateWeightFunctionInternals();

  /// get Nb data for synchronization in parallel
  inline UInt getNbDataForElements(const Array<Element> & elements,
                                   const ID & id) const;

  /// pack data for synchronization in parallel
  inline void packElementData(CommunicationBuffer & buffer,
                              const Array<Element> & elements,
                              SynchronizationTag tag, const ID & id) const;

  /// unpack data for synchronization in parallel
  inline void unpackElementData(CommunicationBuffer & buffer,
                                const Array<Element> & elements,
                                SynchronizationTag tag, const ID & id) const;

protected:
  /// create a new neighborhood for a given domain ID
  void createNeighborhood(const ID & weight_func, const ID & neighborhood);

  /// flatten the material internal fields needed for the non-local computations
  void flattenInternal(ElementTypeMapReal & internal_flat,
                       const GhostType & ghost_type, const ElementKind & kind);

  /// set the values of the jacobians
  void setJacobians(const FEEngine & fe_engine, const ElementKind & kind);

  /// allocation of eelment type maps
  void initElementTypeMap(UInt nb_component, ElementTypeMapReal & element_map,
                          const FEEngine & fe_engine,
                          const ElementKind el_kind = _ek_regular);

  /// resizing of element type maps
  void resizeElementTypeMap(UInt nb_component, ElementTypeMapReal & element_map,
			    const FEEngine & fee,
                            const ElementKind el_kind = _ek_regular);

  /// remove integration points from element type maps
  void removeIntegrationPointsFromMap(
      const ElementTypeMapArray<UInt> & new_numbering, UInt nb_component,
      ElementTypeMapReal & element_map, const FEEngine & fee,
      const ElementKind el_kind = _ek_regular);

  /// allocate the non-local variables
  void initNonLocalVariables();

  /// copy the results of the averaging in the materials
  void distributeInternals(ElementKind kind);

  /// cleanup unneccessary ghosts
  void cleanupExtraGhostElements(ElementTypeMap<UInt> & nb_ghost_protected);

  /* ------------------------------------------------------------------------ */
  /* MeshEventHandler inherited members                                       */
  /* ------------------------------------------------------------------------ */
public:
  virtual void onElementsRemoved(const Array<Element> & element_list,
				 const ElementTypeMapArray<UInt> & new_numbering,
				 const RemovedElementsEvent & event);
  
  virtual void onElementsAdded(const Array<Element> & element_list,
                               const NewElementsEvent & event);

  virtual void onElementsChanged(__attribute__((unused)) const Array<Element> & old_elements_list,
				 __attribute__((unused)) const Array<Element> & new_elements_list,
				 __attribute__((unused)) const ElementTypeMapArray<UInt> & new_numbering,
				 __attribute__((unused)) const ChangedElementsEvent & event) {};

  virtual void onNodesAdded(__attribute__((unused)) const Array<UInt> & nodes_list,
                            __attribute__((unused)) const NewNodesEvent & event) {};
  virtual void onNodesRemoved(__attribute__((unused)) const Array<UInt> & nodes_list,
                              __attribute__((unused)) const Array<UInt> & new_numbering,
                              __attribute__((unused)) const RemovedNodesEvent & event) {};

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  AKANTU_GET_MACRO(SpatialDimension, spatial_dimension, UInt);
  AKANTU_GET_MACRO(Model, model, const SolidMechanicsModel &);
  AKANTU_GET_MACRO_NOT_CONST(Volumes, volumes, ElementTypeMapReal &)
  AKANTU_GET_MACRO(NbStressCalls, compute_stress_calls, UInt);

  inline const Array<Real> & getJacobians(const ElementType & type,
                                          const GhostType & ghost_type) {
    return *jacobians(type, ghost_type);
  }

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// the non-local neighborhoods present
  NeighborhoodMap neighborhoods;

  /// list of all the non-local materials in the model
  std::vector<Material *> non_local_materials;

  struct NonLocalVariable {
    NonLocalVariable(const ID & variable_name, const ID & nl_variable_name,
                     const ID & id, UInt nb_component)
      : local(variable_name, id, 0),
	non_local(nl_variable_name, id, 0),
	nb_component(nb_component) {}
    ElementTypeMapReal local;
    ElementTypeMapReal non_local;
    UInt nb_component;
  };

  /// the non-local variables associated to a certain neighborhood
  std::map<ID, NonLocalVariable *> non_local_variables;

  /// reference to the model
  SolidMechanicsModel & model;

  /// jacobians for all the elements in the mesh
  ElementTypeMap<const Array<Real> *> jacobians;

  /// store the position of the quadrature points
  ElementTypeMapReal quad_positions;

  /// store the volume of each quadrature point for the non-local weight
  /// normalization
  ElementTypeMapReal volumes;

  /// the spatial dimension
  const UInt spatial_dimension;

  /// counter for computeStress calls
  UInt compute_stress_calls;

  /// map to store weight function types from input file
  std::map<ID, ParserSection> weight_function_types;

  /// map to store the internals needed by the weight functions
  std::map<ID, ElementTypeMapReal *> weight_function_internals;
  /* --------------------------------------------------------------------------
   */
  /// the following are members needed to make this processor participate in the
  /// grid creation of neighborhoods he doesn't own as a member. For details see
  /// createGridSynchronizers function

  /// synchronizer registry for dummy grid synchronizers
  SynchronizerRegistry * dummy_registry;

  /// map of dummy synchronizers
  std::map<ID, GridSynchronizer *> dummy_synchronizers;

  /// dummy spatial grid
  SpatialGrid<IntegrationPoint> * dummy_grid;

  /// create a set of all neighborhoods present in the simulation
  std::set<ID> global_neighborhoods;

  class DummyDataAccessor : public DataAccessor {
  public:
    virtual inline UInt getNbDataForElements(__attribute__((unused)) const Array<Element> & elements,
					     __attribute__((unused)) SynchronizationTag tag) const { return 0; };
    
    virtual inline void packElementData(__attribute__((unused)) CommunicationBuffer & buffer,
					__attribute__((unused)) const Array<Element> & element,
					__attribute__((unused)) SynchronizationTag tag) const {};
    
    virtual inline void unpackElementData(__attribute__((unused)) CommunicationBuffer & buffer,
					  __attribute__((unused)) const Array<Element> & element,
					  __attribute__((unused)) SynchronizationTag tag)  {};
  };

  DummyDataAccessor dummy_accessor;
};

} // akantu

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "non_local_manager_inline_impl.cc"

#endif /* __AKANTU_NON_LOCAL_MANAGER_HH__ */

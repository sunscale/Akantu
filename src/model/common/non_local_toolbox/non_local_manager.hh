/**
 * @file   non_local_manager.hh
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @date   Mon Sep 21 14:21:33 2015
 *
 * @brief  Classes that manages all the non-local neighborhoods
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
#ifndef __AKANTU_NON_LOCAL_MANAGER_HH__
#define __AKANTU_NON_LOCAL_MANAGER_HH__
/* -------------------------------------------------------------------------- */
#include "aka_memory.hh"
#include "solid_mechanics_model.hh"
#include "non_local_neighborhood_base.hh"
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

class NonLocalManager : public Memory {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  NonLocalManager(SolidMechanicsModel & model, 
		  const ID & id = "non_local_manager",
		  const MemoryID & memory_id = 0);
  virtual ~NonLocalManager();
  typedef std::map<ID, NonLocalNeighborhoodBase *> NeighborhoodMap;
  typedef std::pair<ID, ID> KeyCOO;


/* -------------------------------------------------------------------------- */
/* Methods                                                                    */
/* -------------------------------------------------------------------------- */
public:
  
  /// initialize the non-local manager: compute pair lists and weights for all neighborhoods
  void init();

  /// insert new quadrature point in the grid
  inline void insertQuad(const QuadraturePoint & quad, const Vector<Real> & coords, Real radius, const ID & type, ID name = "");

  /// return the fem object associated with a provided name
  inline NonLocalNeighborhoodBase & getNeighborhood(const ID & name = "") const;

  /// create the grid synchronizers for each neighborhood
  void createNeighborhoodSynchronizers();

  /// compute the weights in each neighborhood for non-local averaging
  inline void computeWeights();

  /// compute the weights in each neighborhood for non-local averaging
  inline void updatePairLists();

  /// register a new non-local material
  inline void registerNonLocalMaterial(Material & new_mat);

  /// register a non-local variable
  inline void registerNonLocalVariable(const ID & variable_name, const ID & nl_variable_name, UInt nb_component);

  void averageInternals(const GhostType & ghost_type = _not_ghost);

private:

  /// create a new neighborhood for a given domain ID
  void createNeighborhood(const ID & type, Real radius, const ID & name);

  /// flatten the material internal fields needed for the non-local computations
  void flattenInternal(ElementTypeMapReal & internal_flat,
		       const GhostType & ghost_type,
		       const ElementKind & kind);

  /// set the values of the jacobians
  void setJacobians(const FEEngine & fe_engine, const ElementKind & kind);


  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  AKANTU_GET_MACRO(Model, model, const SolidMechanicsModel &);
  AKANTU_GET_MACRO_NOT_CONST(Volumes, volumes, ElementTypeMapReal &)

  inline const Array<Real> & getJacobians(const ElementType & type, const GhostType & ghost_type) {
    return *jacobians(type, ghost_type);
  }

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:

  /// the non-local neighborhoods present
  NeighborhoodMap neighborhoods;

  /// list of all the non-local materials in the model
  std::vector<Material * > non_local_materials;

  struct NonLocalVariable {
    NonLocalVariable(const ID & variable_name, const ID & nl_variable_name, const ID & id, UInt nb_component) :
      local(variable_name, id),
      non_local(nl_variable_name, id),
      nb_component(nb_component){
    } 
    ElementTypeMapReal local;
    ElementTypeMapReal non_local;
    UInt nb_component;
  };

  /// the non-local variables associated to a certain neighborhood
  std::map<ID, NonLocalVariable *> non_local_variables;

  /// reference to the model
  SolidMechanicsModel & model;

  /// jacobians for all the elements in the mesh
  ElementTypeMap<const Array<Real> * > jacobians;

  /// default neighborhood object
  std::string default_neighborhood;

  /// store the position of the quadrature points
  ElementTypeMapReal quad_positions;

  /// store the volume of each quadrature point for the non-local weight normalization
  ElementTypeMapReal volumes; 

};

__END_AKANTU__

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "non_local_manager_inline_impl.cc"


#endif /* __AKANTU_NON_LOCAL_MANAGER_HH__ */

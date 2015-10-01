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
  /// typedef std::map<ID, NonLocalVariable *> NonLocalVariableMap;

/* -------------------------------------------------------------------------- */
/* Methods                                                                    */
/* -------------------------------------------------------------------------- */
public:
  /// insert new quadrature point in the grid
  inline void insertQuad(const QuadraturePoint & quad, const ID & id = "");

  /// return the fem object associated with a provided name
  inline NonLocalNeighborhoodBase & getNeighborhood(const ID & name = "") const;

  /// create a new neighborhood for a given domain ID
  void createNeighborhood(const ID & type, Real radius, const ID & name = "");

  /// set the values of the jacobians
  void setJacobians(const FEEngine & fe_engine, const ElementKind & kind);

  /// create the grid synchronizers for each neighborhood
  void createNeighborhoodSynchronizers();

  /// compute the weights in each neighborhood for non-local averaging
  inline void computeWeights();

  /// compute the weights in each neighborhood for non-local averaging
  inline void updatePairLists();

private:

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

  struct NonLocalVariable {
    ElementTypeMap<Real> * local;
    ElementTypeMap<Real> * non_local;
    UInt nb_component;
  };

  /// the non-local variables associated to a certain neighborhood
  std::map<ID, NonLocalVariable> non_local_variables;

  /// reference to the model
  SolidMechanicsModel & model;

  /// jacobians for all the elements in the mesh
  ElementTypeMap<const Array<Real> * > jacobians;

  /// default neighborhood object
  std::string default_neighborhood;

  /// store the volume of each quadrature point for the non-local weight normalization
  ElementTypeMapReal volumes; 
};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "non_local_manager_inline_impl.cc"

__END_AKANTU__

#endif /* __AKANTU_NON_LOCAL_MANAGER_HH__ */

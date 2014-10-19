/**
 * @file   boundary_condition.hh
 *
 * @author Dana Christen <dana.christen@gmail.com>
 *
 * @date creation: Fri May 03 2013
 * @date last modification: Thu Jun 05 2014
 *
 * @brief  XXX
 *
 * @section LICENSE
 *
 * Copyright (©) 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

#ifndef __AKANTU_BOUNDARY_CONDITION_HH__
#define __AKANTU_BOUNDARY_CONDITION_HH__

#include "aka_common.hh"
#include "boundary_condition_functor.hh"
#include "mesh.hh"
#include "fe_engine.hh"
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

template <class ModelType>
class BoundaryCondition {

  /* ------------------------------------------------------------------------ */
  /* Typedefs                                                                 */
  /* ------------------------------------------------------------------------ */
private:

  /* ------------------------------------------------------------------------ */
  /* Constructors / Destructors / Initializers                                */
  /* ------------------------------------------------------------------------ */
public:
  BoundaryCondition() : model(NULL), primal(NULL), dual(NULL), primal_increment(NULL) {}

  void initBC(ModelType & ptr, Array<Real> & primal, Array<Real> & dual);
  void initBC(ModelType & ptr, Array<Real> & primal,
	      Array<Real> & primal_increment, Array<Real> & dual);
  /* ------------------------------------------------------------------------ */
  /* Methods and accessors                                                    */
  /* ------------------------------------------------------------------------ */
public:

  //inline void initBoundaryCondition();
  template<typename FunctorType>
  inline void applyBC(const FunctorType & func);

  template<class FunctorType>
  inline void applyBC(const FunctorType & func, const std::string & group_name);

  template<class FunctorType>
  inline void applyBC(const FunctorType & func, const ElementGroup & element_group);

  AKANTU_GET_MACRO_NOT_CONST(Model, *model, ModelType &);
  AKANTU_GET_MACRO_NOT_CONST(Primal,*primal, Array<Real> &);
  AKANTU_GET_MACRO_NOT_CONST(Dual,  *dual,   Array<Real> &);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */

public:

  template<class FunctorType, BC::Functor::Type type = FunctorType::type>
  struct TemplateFunctionWrapper;

private:

  ModelType * model;
  Array<Real> * primal;
  Array<Real> * dual;
  Array<Real> * primal_increment;
};

__END_AKANTU__

#include "boundary_condition_tmpl.hh"

#endif /* __AKANTU_BOUNDARY_CONDITION_HH__ */


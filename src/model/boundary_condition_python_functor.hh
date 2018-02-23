/**
 * @file   boundary_condition_python_functor.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Thu Feb 21 2013
 * @date last modification: Wed Jan 31 2018
 *
 * @brief  interface for BC::Functor writen in python
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
#include "aka_common.hh"
#include "boundary_condition_functor.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_BOUNDARY_CONDITION_PYTHON_FUNCTOR_HH__
#define __AKANTU_BOUNDARY_CONDITION_PYTHON_FUNCTOR_HH__
/* -------------------------------------------------------------------------- */
#include "boundary_condition_functor.hh"
#include "python_functor.hh"
/* -------------------------------------------------------------------------- */
namespace akantu {

namespace BC {

  class PythonFunctorDirichlet : public PythonFunctor, public Functor {

    /* ------------------------------------------------------------------------
     */
    /* Constructors/Destructors */
    /* ------------------------------------------------------------------------
     */

  public:
    PythonFunctorDirichlet(PyObject * obj) : PythonFunctor(obj) {}

    /* ------------------------------------------------------------------------
     */
    /* Methods */
    /* ------------------------------------------------------------------------
     */

  public:
    void operator()(UInt node, Vector<bool> & flags, Vector<Real> & primal,
                    const Vector<Real> & coord) const;

    /* ------------------------------------------------------------------------
     */
    /* Class Members */
    /* ------------------------------------------------------------------------
     */

  public:
    static const Type type = _dirichlet;
  };

  /* --------------------------------------------------------------------------
   */

  class PythonFunctorNeumann : public PythonFunctor, public Functor {

    /* ------------------------------------------------------------------------
     */
    /* Constructors/Destructors */
    /* ------------------------------------------------------------------------
     */

  public:
    PythonFunctorNeumann(PyObject * obj) : PythonFunctor(obj) {}

    /* ------------------------------------------------------------------------
     */
    /* Methods */
    /* ------------------------------------------------------------------------
     */

  public:
    void operator()(const IntegrationPoint & quad_point, Vector<Real> & dual,
                    const Vector<Real> & coord,
                    const Vector<Real> & normals) const;

    /* ------------------------------------------------------------------------
     */
    /* Class Members */
    /* ------------------------------------------------------------------------
     */

  public:
    static const Type type = _neumann;
  };

} // end namespace BC

} // akantu

#endif /* __AKANTU_BOUNDARY_CONDITION_PYTHON_FUNCTOR_HH__ */

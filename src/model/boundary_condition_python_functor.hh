/**
 * @file   boundary_condition_python_functor.hh
 *

 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 *
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
#include "aka_common.hh"
#include "boundary_condition_functor.hh"
/* -------------------------------------------------------------------------- */
#include <Python.h>
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_BOUNDARY_CONDITION_PYTHON_FUNCTOR_HH__
#define __AKANTU_BOUNDARY_CONDITION_PYTHON_FUNCTOR_HH__

__BEGIN_AKANTU__


namespace BC {

class PythonFunctor : public Functor {

public:
  PythonFunctor(PyObject * obj) : python_obj(obj) {}

private:
  PyObject * python_obj;

};

class PythonFunctorDirichlet : public PythonFunctor {

public:
  PythonFunctorDirichlet(PyObject * obj) : PythonFunctor(obj) {}

public:
  void operator()(UInt node,
		  Vector<bool> & flags,
		  Vector<Real> & primal,
		  const Vector<Real> & coord) const;

public:
  static const Type type = _dirichlet;

  
};


    
}//end namespace BC


__END_AKANTU__

#endif /* __AKANTU_BOUNDARY_CONDITION_PYTHON_FUNCTOR_HH__ */

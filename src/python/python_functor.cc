/**
 * @file   python_functor.cc
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
#include "python_functor.hh"
#include "aka_common.hh"
/* -------------------------------------------------------------------------- */
__BEGIN_AKANTU__
/* -------------------------------------------------------------------------- */  

PythonFunctor::PythonFunctor(PyObject * obj):python_obj(obj){
}

/* -------------------------------------------------------------------------- */
  

PyObject * PythonFunctor::callFunctor(PyObject * functor,
				      PyObject * args,
				      PyObject * kwargs) const{
  
  
  if (!PyCallable_Check(functor))
    AKANTU_EXCEPTION("Provided functor is not a function");
  
  PyObject * pValue = PyObject_Call(functor, args,kwargs);
  
  PyObject* exception_type = PyErr_Occurred();
  if (exception_type){
    PyObject * exception;
    PyObject * traceback;
    PyErr_Fetch(&exception_type, &exception, &traceback);
    
    PyObject_Print(exception_type, stdout, Py_PRINT_RAW);
    PyObject_Print(exception, stdout, Py_PRINT_RAW);
    std::stringstream sstr;
    sstr << "Exception occured while calling the functor: ";
    
    PyObject * exception_mesg = PyObject_GetAttrString(exception,"message");
    if (exception_mesg && PyString_Check(exception_mesg))
      sstr << PyString_AsString(exception_mesg);
    else
      sstr << PyString_AsString(exception);
    
    AKANTU_EXCEPTION(sstr.str());
  }

  return pValue;
}

/* -------------------------------------------------------------------------- */

__END_AKANTU__
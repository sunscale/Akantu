/**
 * @file   python_functor.hh
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


#ifndef __AKANTU_PYTHON_FUNCTOR_HH__
#define __AKANTU_PYTHON_FUNCTOR_HH__
/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include <Python.h>
/* -------------------------------------------------------------------------- */
__BEGIN_AKANTU__
/* -------------------------------------------------------------------------- */


class PythonFunctor {


  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */

public:
  PythonFunctor(PyObject * obj);


  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */

protected:

  /// call the python functor
  PyObject * callFunctor(PyObject * functor,
			 PyObject * args,
			 PyObject * kwargs) const;

  /// call the python functor from variadic types
  template <typename return_type, typename... Params>
  return_type callFunctor(const std::string & functor_name,
			  Params... parameters) const;

  /// empty function to cose the recursive template loop
  inline void packArguments(std::vector<PyObject*> & pArgs) const;
  
  /// get the python function object
  inline PyObject * getPythonFunction(const std::string & functor_name) const;


  ///
  template<typename T, typename... Args>
  inline void packArguments(std::vector<PyObject*> & pArgs,
			    T & p,
			    Args&...params)const;

  /// convert an akantu object to python 
  template <typename T>
  inline PyObject *convertToNumpy(const T & akantu_obj) const;

  /// convert an akantu vector to python 
  template <typename T>
  inline PyObject *convertToNumpy(const Vector<T> & akantu_obj) const;

  /// returns the numpy data type code
  template <typename T>
  inline int getPythonDataTypeCode() const;

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */

  
  PyObject * python_obj;

};
/* -------------------------------------------------------------------------- */
__END_AKANTU__
/* -------------------------------------------------------------------------- */
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "python_functor_inline_impl.cc"
/* -------------------------------------------------------------------------- */


#endif /* __AKANTU_PYTHON_FUNCTOR_HH__ */

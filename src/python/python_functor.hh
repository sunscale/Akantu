/**
 * @file   python_functor.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Sun Nov 15 2015
 *
 * @brief  Python Functor interface
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
#include "aka_common.hh"
#include "internal_field.hh"
/* -------------------------------------------------------------------------- */
#include <Python.h>
#include <map>
#include <vector>
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_PYTHON_FUNCTOR_HH__
#define __AKANTU_PYTHON_FUNCTOR_HH__

namespace akantu {

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
  PyObject * callFunctor(PyObject * functor, PyObject * args,
                         PyObject * kwargs) const;

  /// call the python functor from variadic types
  template <typename return_type, typename... Params>
  return_type callFunctor(const std::string & functor_name,
                          Params &... parameters) const;

  /// empty function to cose the recursive template loop
  inline void packArguments(std::vector<PyObject *> & pArgs) const;

  /// get the python function object
  inline PyObject * getPythonFunction(const std::string & functor_name) const;

  /// variadic template for unknown number of arguments to unpack
  template <typename T, typename... Args>
  inline void packArguments(std::vector<PyObject *> & pArgs, T & p,
                            Args &... params) const;

  /// convert an akantu object to python
  template <typename T>
  inline PyObject * convertToPython(const T & akantu_obj) const;

  /// convert a stl vector to python
  template <typename T>
  inline PyObject * convertToPython(const std::vector<T> & akantu_obj) const;

  /// convert a stl vector to python
  template <typename T>
  inline PyObject * convertToPython(const std::vector<Array<T> *> & akantu_obj) const;

  /// convert a stl vector to python
  template <typename T>
  inline PyObject * convertToPython(const std::unique_ptr<T> & akantu_obj) const;
  
  /// convert a stl vector to python
  template <typename T1, typename T2>
  inline PyObject * convertToPython(const std::map<T1, T2> & akantu_obj) const;
  
  /// convert an akantu vector to python
  template <typename T>
  inline PyObject * convertToPython(const Vector<T> & akantu_obj) const;

  /// convert an akantu internal to python
  template <typename T>
  inline PyObject * convertToPython(const InternalField<T> & akantu_obj) const;
  
  /// convert an akantu vector to python
  template <typename T>
  inline PyObject * convertToPython(const Array<T> & akantu_obj) const;

  /// convert an akantu vector to python
  template <typename T>
  inline PyObject * convertToPython(Array<T> * akantu_obj) const;

  /// convert a akantu matrix to python
  template <typename T>
  inline PyObject * convertToPython(const Matrix<T> & akantu_obj) const;

  /// convert a python object to an akantu object
  template <typename return_type>
  inline return_type convertToAkantu(PyObject * python_obj) const;

  /// convert a python object to an akantu object
  template <typename T>
  inline std::vector<T> convertListToAkantu(PyObject * python_obj) const;

  /// returns the numpy data type code
  template <typename T> inline int getPythonDataTypeCode() const;

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  PyObject * python_obj;
};

} // akantu

#endif /* __AKANTU_PYTHON_FUNCTOR_HH__ */

#include "python_functor_inline_impl.cc"

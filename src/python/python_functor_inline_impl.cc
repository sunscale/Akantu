/**
 * @file   python_functor_inline_impl.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 *
 * @date creation: Fri Nov 13 2015
 * @date last modification: Tue Feb 20 2018
 *
 * @brief  Python functor interface
 *
 * @section LICENSE
 *
 * Copyright (©) 2015-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "integration_point.hh"
#include "internal_field.hh"
/* -------------------------------------------------------------------------- */
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
#include <numpy/arrayobject.h>
#include <typeinfo>
#if PY_MAJOR_VERSION >= 3
#include <codecvt>
#endif
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_PYTHON_FUNCTOR_INLINE_IMPL_CC__
#define __AKANTU_PYTHON_FUNCTOR_INLINE_IMPL_CC__

namespace akantu {

/* -------------------------------------------------------------------------- */
template <typename T> inline int PythonFunctor::getPythonDataTypeCode() {
  AKANTU_EXCEPTION("undefined type: " << debug::demangle(typeid(T).name()));
}

/* -------------------------------------------------------------------------- */
template <> inline int PythonFunctor::getPythonDataTypeCode<bool>() {
  int data_typecode = NPY_NOTYPE;
  size_t s = sizeof(bool);
  switch (s) {
  case 1:
    data_typecode = NPY_BOOL;
    break;
  case 2:
    data_typecode = NPY_UINT16;
    break;
  case 4:
    data_typecode = NPY_UINT32;
    break;
  case 8:
    data_typecode = NPY_UINT64;
    break;
  }
  return data_typecode;
}

/* -------------------------------------------------------------------------- */
template <> inline int PythonFunctor::getPythonDataTypeCode<double>() {
  return NPY_DOUBLE;
}

/* -------------------------------------------------------------------------- */
template <> inline int PythonFunctor::getPythonDataTypeCode<UInt>() {
  return NPY_UINT;
}

/* -------------------------------------------------------------------------- */
template <>
inline PyObject *
PythonFunctor::convertToPython<double>(const double & akantu_object) {
  return PyFloat_FromDouble(akantu_object);
}

/* -------------------------------------------------------------------------- */
template <>
inline PyObject *
PythonFunctor::convertToPython<UInt>(const UInt & akantu_object) {
#if PY_MAJOR_VERSION >= 3
  return PyLong_FromLong(akantu_object);
#else
  return PyInt_FromLong(akantu_object);
#endif
}

/* -------------------------------------------------------------------------- */

template <>
inline PyObject *
PythonFunctor::convertToPython<int>(const int & akantu_object) {
#if PY_MAJOR_VERSION >= 3
  return PyLong_FromLong(akantu_object);
#else
  return PyInt_FromLong(akantu_object);
#endif
}

/* -------------------------------------------------------------------------- */

template <>
inline PyObject *
PythonFunctor::convertToPython<bool>(const bool & akantu_object) {
  return PyBool_FromLong(long(akantu_object));
}
/* -------------------------------------------------------------------------- */

template <>
inline PyObject *
PythonFunctor::convertToPython<std::string>(const std::string & str) {
#if PY_MAJOR_VERSION >= 3
  return PyUnicode_FromString(str.c_str());
#else
  return PyString_FromString(str.c_str());
#endif
}

/* --------------------------------------------------------------------------
 */

template <>
inline PyObject *
PythonFunctor::convertToPython<NodeGroup>(const NodeGroup & group) {
  return PythonFunctor::convertToPython(group.getNodes());
}

/* --------------------------------------------------------------------------
 */

template <>
inline PyObject *
PythonFunctor::convertToPython<ElementGroup>(const ElementGroup & group) {

  PyObject * res = PyDict_New();

  PyDict_SetItem(res,
                 PythonFunctor::convertToPython(std::string("element_group")),
                 PythonFunctor::convertToPython(group.getElements()));

  PyDict_SetItem(res, PythonFunctor::convertToPython(std::string("node_group")),
                 PythonFunctor::convertToPython(group.getNodeGroup()));

  return res;
}

/* --------------------------------------------------------------------------
 */

template <>
inline PyObject *
PythonFunctor::convertToPython<ElementGroup *>(ElementGroup * const & group) {
  return PythonFunctor::convertToPython(*group);
}

/* --------------------------------------------------------------------------
 */

template <>
inline PyObject *
PythonFunctor::convertToPython<ElementType>(const ElementType & type) {
  // std::stringstream sstr;
  // sstr << type;
  // return PythonFunctor::convertToPython(sstr.str());
  return PythonFunctor::convertToPython(int(type));
}

/* --------------------------------------------------------------------------
 */

template <typename T>
inline PyObject *
PythonFunctor::convertToPython(const ElementTypeMapArray<T> & map) {

  std::map<std::string, Array<T> *> res;
  for (const auto & type : map.elementTypes()) {
    std::stringstream sstr;
    sstr << type;
    res[sstr.str()] = const_cast<Array<T> *>(&map(type));
  }
  return PythonFunctor::convertToPython(res);
}

/* --------------------------------------------------------------------------
 */

template <typename T> PyObject * PythonFunctor::convertToPython(const T &) {
  AKANTU_EXCEPTION(__PRETTY_FUNCTION__ << " : not implemented yet !"
                                       << std::endl
                                       << debug::demangle(typeid(T).name()));
}

/* -------------------------------------------------------------------------- */
template <typename T>
inline PyObject * PythonFunctor::convertToPython(const std::vector<T> & array) {
  int data_typecode = getPythonDataTypeCode<T>();
  npy_intp dims[1] = {int(array.size())};
  PyObject * obj = PyArray_SimpleNewFromData(1, dims, data_typecode,
                                             const_cast<T *>(&array[0]));
  auto * res = (PyArrayObject *)obj;
  return (PyObject *)res;
}
/* -------------------------------------------------------------------------- */
template <typename T>
inline PyObject *
PythonFunctor::convertToPython(const std::vector<Array<T> *> & array) {

  PyObject * res = PyDict_New();

  for (auto a : array) {
    PyObject * obj = PythonFunctor::convertToPython(*a);
    PyObject * name = PythonFunctor::convertToPython(a->getID());
    PyDict_SetItem(res, name, obj);
  }
  return (PyObject *)res;
}
/* -------------------------------------------------------------------------- */

template <typename T1, typename T2>
inline PyObject * PythonFunctor::convertToPython(const std::map<T1, T2> & map) {

  PyObject * res = PyDict_New();

  for (auto && a : map) {
    PyObject * key = PythonFunctor::convertToPython(a.first);
    PyObject * value = PythonFunctor::convertToPython(a.second);
    PyDict_SetItem(res, key, value);
  }
  return (PyObject *)res;
}

/* -------------------------------------------------------------------------- */
template <typename T>
PyObject * PythonFunctor::convertToPython(const Vector<T> & array) {
  int data_typecode = getPythonDataTypeCode<T>();
  npy_intp dims[1] = {array.size()};
  PyObject * obj =
      PyArray_SimpleNewFromData(1, dims, data_typecode, array.storage());
  auto * res = (PyArrayObject *)obj;
  return (PyObject *)res;
}

/* -------------------------------------------------------------------------- */
template <typename T>
inline PyObject *
PythonFunctor::convertToPython(const InternalField<T> & internals) {
  return convertToPython(
      static_cast<const ElementTypeMapArray<T> &>(internals));
}

/* --------------------------------------------------------------------- */

template <typename T>
inline PyObject * PythonFunctor::convertToPython(const std::unique_ptr<T> & u) {
  return convertToPython(*u);
}

/* --------------------------------------------------------------------- */

template <typename T>
PyObject * PythonFunctor::convertToPython(const Array<T> & array) {
  int data_typecode = getPythonDataTypeCode<T>();
  npy_intp dims[2] = {array.size(), array.getNbComponent()};
  PyObject * obj =
      PyArray_SimpleNewFromData(2, dims, data_typecode, array.storage());
  auto * res = (PyArrayObject *)obj;
  return (PyObject *)res;
}
/* ---------------------------------------------------------------------- */
template <typename T>
PyObject * PythonFunctor::convertToPython(Array<T> * array) {
  return PythonFunctor::convertToPython(*array);
}

/* ---------------------------------------------------------------------- */
template <typename T>
PyObject * PythonFunctor::convertToPython(const Matrix<T> & mat) {
  int data_typecode = getPythonDataTypeCode<T>();
  npy_intp dims[2] = {mat.size(0), mat.size(1)};
  PyObject * obj =
      PyArray_SimpleNewFromData(2, dims, data_typecode, mat.storage());
  auto * res = (PyArrayObject *)obj;
  return (PyObject *)res;
}

/* ---------------------------------------------------------------------- */

template <>
inline PyObject *
PythonFunctor::convertToPython<IntegrationPoint>(const IntegrationPoint & qp) {

  PyObject * input = PyDict_New();
  PyObject * num_point = PythonFunctor::convertToPython(qp.num_point);
  PyObject * global_num = PythonFunctor::convertToPython(qp.global_num);
  PyObject * material_id = PythonFunctor::convertToPython(qp.material_id);
  PyObject * position = PythonFunctor::convertToPython(qp.getPosition());
  PyDict_SetItemString(input, "num_point", num_point);
  PyDict_SetItemString(input, "global_num", global_num);
  PyDict_SetItemString(input, "material_id", material_id);
  PyDict_SetItemString(input, "position", position);
  return input;
}

/* --------------------------------------------------------------------------
 */
inline PyObject *
PythonFunctor::getPythonFunction(const std::string & functor_name) const {
#if PY_MAJOR_VERSION < 3
  if (!PyInstance_Check(this->python_obj))
    AKANTU_EXCEPTION("Python object is not an instance");
#else
// does not make sense to check everything is an instance of object in python 3
#endif

  if (not PyObject_HasAttrString(this->python_obj, functor_name.c_str()))
    AKANTU_EXCEPTION("Python dictionary has no " << functor_name << " entry");

  PyObject * pFunctor =
      PyObject_GetAttrString(this->python_obj, functor_name.c_str());

  return pFunctor;
}

/* --------------------------------------------------------------------------
 */
inline void PythonFunctor::packArguments(__attribute__((unused))
                                         std::vector<PyObject *> & p_args) {}

/* --------------------------------------------------------------------------
 */
template <typename T, typename... Args>
inline void PythonFunctor::packArguments(std::vector<PyObject *> & p_args,
                                         T & p, Args &... params) {
  p_args.push_back(PythonFunctor::convertToPython(p));
  if (sizeof...(params) != 0)
    PythonFunctor::packArguments(p_args, params...);
}

/* --------------------------------------------------------------------------
 */
template <typename return_type, typename... Params>
return_type PythonFunctor::callFunctor(const std::string & functor_name,
                                       Params &... parameters) const {
  _import_array();

  std::vector<PyObject *> arg_vector;
  this->packArguments(arg_vector, parameters...);

  PyObject * pArgs = PyTuple_New(arg_vector.size());
  for (UInt i = 0; i < arg_vector.size(); ++i) {
    PyTuple_SetItem(pArgs, i, arg_vector[i]);
  }

  PyObject * kwargs = PyDict_New();

  PyObject * pFunctor = this->getPythonFunction(functor_name);
  PyObject * res = this->callFunctor(pFunctor, pArgs, kwargs);

  Py_XDECREF(pArgs);
  Py_XDECREF(kwargs);

  return this->convertToAkantu<return_type>(res);
}

/* --------------------------------------------------------------------------
 */
template <typename return_type>
inline return_type PythonFunctor::convertToAkantu(PyObject * python_obj) {

  if (PyList_Check(python_obj)) {
    return PythonFunctor::convertListToAkantu<typename return_type::value_type>(
        python_obj);
  }
  AKANTU_TO_IMPLEMENT();
}

/* --------------------------------------------------------------------------
 */
template <>
inline void PythonFunctor::convertToAkantu<void>(PyObject * python_obj) {
  if (python_obj != Py_None)
    AKANTU_DEBUG_WARNING(
        "functor return a value while none was expected: ignored");
}

/* --------------------------------------------------------------------------
 */
template <>
inline std::string
PythonFunctor::convertToAkantu<std::string>(PyObject * python_obj) {
#if PY_MAJOR_VERSION >= 3
  if (!PyUnicode_Check(python_obj))
    AKANTU_EXCEPTION("cannot convert object to string");

  std::wstring unicode_str(PyUnicode_AsWideCharString(python_obj, NULL));
  std::wstring_convert<std::codecvt_utf8<wchar_t>, wchar_t> converter;
  return converter.to_bytes(unicode_str);
#else
  if (!PyString_Check(python_obj))
    AKANTU_EXCEPTION("cannot convert object to string");

  return PyString_AsString(python_obj);
#endif
}

/* --------------------------------------------------------------------------
 */
template <>
inline Real PythonFunctor::convertToAkantu<Real>(PyObject * python_obj) {
  if (!PyFloat_Check(python_obj))
    AKANTU_EXCEPTION("cannot convert object to float");
  return PyFloat_AsDouble(python_obj);
}
/* --------------------------------------------------------------------------
 */
template <>
inline UInt PythonFunctor::convertToAkantu<UInt>(PyObject * python_obj) {
#if PY_MAJOR_VERSION >= 3
  if (!PyLong_Check(python_obj))
    AKANTU_EXCEPTION("cannot convert object to integer");
  return PyLong_AsLong(python_obj);
#else
  if (!PyInt_Check(python_obj))
    AKANTU_EXCEPTION("cannot convert object to integer");
  return PyInt_AsLong(python_obj);
#endif
}

/* --------------------------------------------------------------------------
 */
template <typename T>
inline std::vector<T>
PythonFunctor::convertListToAkantu(PyObject * python_obj) {
  std::vector<T> res;
  UInt size = PyList_Size(python_obj);
  for (UInt i = 0; i < size; ++i) {
    PyObject * item = PyList_GET_ITEM(python_obj, i);
    res.push_back(PythonFunctor::convertToAkantu<T>(item));
  }
  return res;
}

/* --------------------------------------------------------------------------
 */

} // akantu

#endif /* __AKANTU_PYTHON_FUNCTOR_INLINE_IMPL_CC__ */

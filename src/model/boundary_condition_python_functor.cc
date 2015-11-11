/**
 * @file   boundary_condition_python_functor.cc
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
#include "boundary_condition_python_functor.hh"
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__


namespace BC {

  PythonFunctor::PythonFunctor(PyObject * obj) : python_obj(obj) {
    import_array();
  }

  
  template <typename T> int getPythonDataTypeCode(){ AKANTU_EXCEPTION("undefined type");}
  template <> int getPythonDataTypeCode<bool>(){
    int data_typecode = NPY_NOTYPE;
    size_t s = sizeof(bool);
    switch(s) {
    case 1: data_typecode = NPY_BOOL; break;
    case 2: data_typecode = NPY_UINT16; break;
    case 4: data_typecode = NPY_UINT32; break;
    case 8: data_typecode = NPY_UINT64; break;
    }
    return data_typecode;
  }

  template <> int getPythonDataTypeCode<double>(){return NPY_DOUBLE;}

  
  template <typename T>
  PyObject *PythonFunctor::convertArrayToNumpy(const Vector<T> & array) const{

    int data_typecode = getPythonDataTypeCode< T >();
    npy_intp dims[1] = {array.size()};
    PyObject* obj = PyArray_SimpleNewFromData(1, dims, data_typecode, array.storage());
    PyArrayObject* res = (PyArrayObject*) obj;
    return (PyObject*)res;
  }

  
  
  void PythonFunctor::callFunctor(PyObject * functor,
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
    
  }


  void PythonFunctorDirichlet::operator ()(UInt node,
					   Vector<bool> & flags,
					   Vector<Real> & primal,
					   const Vector<Real> & coord) const{

    if (!PyInstance_Check(this->python_obj))
      AKANTU_EXCEPTION("Python object is not an instance");
    
    // if (!PyDict_Check(this->python_obj))
    //   AKANTU_EXCEPTION("Python object is not a dictionary");

    //PyObject * keys = PyDict_Keys(dict);
    //PyObject_Print(keys, stdout, Py_PRINT_RAW);    

    PyObject * pFunctor = PyObject_GetAttrString(this->python_obj, "operator");
    if (!pFunctor) AKANTU_EXCEPTION("Python dictionary has no 'functor' entry");

    PyObject * pArgs = PyTuple_New(4);
    PyTuple_SetItem(pArgs,0,PyInt_FromLong(node));
    PyTuple_SetItem(pArgs,1,this->convertArrayToNumpy(flags));
    PyTuple_SetItem(pArgs,2,this->convertArrayToNumpy(primal));
    PyTuple_SetItem(pArgs,3,this->convertArrayToNumpy(coord));
    
    PyObject * kwargs = PyDict_New();

    this->callFunctor(pFunctor,pArgs,kwargs);

    //Py_XDECREF(pArgs);
    //Py_XDECREF(kwargs);

  }


  
  void PythonFunctorNeumann::operator()(const IntegrationPoint & quad_point,
					Vector<Real> & dual,
					const Vector<Real> & coord,
					const Vector<Real> & normals) const{
    throw;
  }

  
}//end namespace BC


__END_AKANTU__


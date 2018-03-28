/**
 * @file   aka_common.i
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Dec 12 2014
 * @date last modification: Wed Jan 13 2016
 *
 * @brief  wrapper to aka_common.hh
 *
 * @section LICENSE
 *
 * Copyright (©) 2015 EPFL (Ecole Polytechnique Fédérale de Lausanne) Laboratory
 * (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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

%{
//#include "aka_config.hh"
#include "aka_common.hh"
#include "aka_csr.hh"
#include "element.hh"
#include "python_functor.hh"
#include "parser.hh"
%}

namespace akantu {
  %ignore getUserParser;
  %ignore ghost_types;
  %ignore initialize(int & argc, char ** & argv);
  %ignore initialize(const std::string & input_file, int & argc, char ** & argv);
  %ignore finalize;
  extern const Array<UInt> empty_filter;
}

%include "stl.i"
%include "parser.i"


%typemap(in) (int argc, char *argv[]) {
  int i = 0;
  if (!PyList_Check($input)) {
    PyErr_SetString(PyExc_ValueError, "Expecting a list");
    return NULL;
  }

  $1 = PyList_Size($input);
  $2 = new char *[$1+1];

  for (i = 0; i < $1; i++) {
    PyObject *s = PyList_GetItem($input,i);
    if (!PyString_Check(s)) {
        free($2);
        PyErr_SetString(PyExc_ValueError, "List items must be strings");
        return NULL;
    }
    $2[i] = PyString_AsString(s);
  }
  $2[i] = 0;
}

%typemap(freearg) (int argc, char *argv[]) {
%#if defined(__INTEL_COMPILER)
//#pragma warning ( disable : 383 )
%#elif defined (__clang__) // test clang to be sure that when we test for gnu it is only gnu
%#elif (defined(__GNUC__) || defined(__GNUG__))
%#  if __cplusplus > 199711L
%#    pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
%#  endif
%#endif

  delete [] $2;

%#if defined(__INTEL_COMPILER)
//#pragma warning ( disable : 383 )
%#elif defined (__clang__) // test clang to be sure that when we test for gnu it is only gnu
%#elif (defined(__GNUC__) || defined(__GNUG__))
%#  if __cplusplus > 199711L
%#    pragma GCC diagnostic pop
%#  endif
%#endif
}

%inline %{
  namespace akantu {
#if defined(AKANTU_USE_MPI)
    const int MPI=1;
#else
    const int MPI=0;
#endif
    void _initializeWithoutArgv(const std::string & input_file) {
      int nb_args = 0;
      char ** null = nullptr;
      initialize(input_file, nb_args, null);
    }

  void __finalize() {
    finalize();
  }
  
#define AKANTU_PP_STR_TO_TYPE2(s, data, elem)                    \
    ({ BOOST_PP_STRINGIZE(elem) , elem})

    PyObject * getElementTypes(){

    std::map<std::string, akantu::ElementType> element_types{
      BOOST_PP_SEQ_FOR_EACH_I(
          AKANTU_PP_ENUM, BOOST_PP_SEQ_SIZE(AKANTU_ek_regular_ELEMENT_TYPE),
          BOOST_PP_SEQ_TRANSFORM(AKANTU_PP_STR_TO_TYPE2, akantu, AKANTU_ek_regular_ELEMENT_TYPE))};

    return akantu::PythonFunctor::convertToPython(element_types);
    }
  }
%}

%pythoncode %{
  import sys as _aka_sys
  def __initialize(input_file="", argv=_aka_sys.argv):
      if _aka_sys.modules[__name__].MPI == 1:
         try:
           from mpi4py import MPI
         except ImportError:
           pass

      _initializeWithoutArgv(input_file)

  def initialize(*args, **kwargs):
      raise RuntimeError("No need to call initialize,"
                      " use parseInput to read an input file")

  def finalize(*args, **kwargs):
      raise RuntimeError("No need to call finalize")

  def parseInput(input_file):
      getStaticParser().parse(input_file)
%}

%include "aka_config.hh"
%include "aka_common.hh"
%include "aka_element_classes_info.hh"
%include "element.hh"
%include "aka_error.hh"

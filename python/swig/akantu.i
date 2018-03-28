/**
 * @file   akantu.i
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Dec 12 2014
 * @date last modification: Mon Nov 23 2015
 *
 * @brief  Main swig file for akantu' python interface
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

%module akantu

%exception {
  try {
    $action
  } catch (akantu::debug::Exception e) {
    PyErr_SetString(PyExc_IndexError,e.what());
    return NULL;
  }
 }

#define __attribute__(x)

%ignore akantu::operator <<;

%include "aka_common.i"
%include "aka_csr.i"
%include "aka_array.i"

%define print_self(MY_CLASS)
  %extend akantu::MY_CLASS {
    std::string __str__() {
      std::stringstream sstr;
      sstr << *($self);
      return sstr.str();
    }
 }
%enddef

%include "mesh.i"
%include "mesh_utils.i"
%include "model.i"
%include "communicator.i"

%include "solid_mechanics_model.i"
#if defined(AKANTU_COHESIVE_ELEMENT)
%include "solid_mechanics_model_cohesive.i"
#endif

#if defined(AKANTU_HEAT_TRANSFER)
%include "heat_transfer_model.i"
#endif


#if defined(AKANTU_STRUCTURAL_MECHANICS)
%include "load_functions.i"
%include "structural_mechanics_model.i"
#endif

%pythoncode %{
  __initialize()
%}

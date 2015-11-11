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

__BEGIN_AKANTU__


namespace BC {

  void PythonFunctorDirichlet::operator ()(UInt node,
					   Vector<bool> & flags,
					   Vector<Real> & primal,
					   const Vector<Real> & coord) const{
    throw;
  }
    
}//end namespace BC


__END_AKANTU__


/**
 * @file   aka_optimize.cc
 *
 * @author Alejandro M. Aragón <alejandro.aragon@epfl.ch>
 *
 * @date creation: Fri Jan 04 2013
 * @date last modification: Tue May 13 2014
 *
 * @brief  Objects that can be used to carry out optimization
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

#include "aka_common.hh"
#include "aka_optimize.hh"

__BEGIN_AKANTU__

std::ostream& operator<<(std::ostream& os, nlopt::result r) {
  
  switch (r) {
      
    case 1:
      os<<"NLOPT_SUCCESS: Generic success return value."<<endl;
      break;
    case 2:
      os<<"NLOPT_STOPVAL_REACHED: Optimization stopped because stopval (above) was reached."<<endl;
      break;
    case 3:
      os<<"NLOPT_FTOL_REACHED: Optimization stopped because ftol_rel or ftol_abs (above) was reached."<<endl;
      break;
    case 4:
      os<<"NLOPT_XTOL_REACHED: Optimization stopped because xtol_rel or xtol_abs (above) was reached."<<endl;
      break;
    case 5:
      os<<"NLOPT_MAXEVAL_REACHED: Optimization stopped because maxeval (above) was reached."<<endl;
      break;
    case 6:
      os<<"NLOPT_MAXTIME_REACHED: Optimization stopped because maxtime (above) was reached."<<endl;
      break;
    case -1:
      os<<"NLOPT_FAILURE: Generic failure code."<<endl;
      break;
    case -2:
      os<<"NLOPT_INVALID_ARGS: Invalid arguments (e.g. lower bounds are bigger than upper bounds, an unknown algorithm was specified, etc.)"<<endl;
      break;
    case -3:
      os<<"NLOPT_OUT_OF_MEMORY: Ran out of memory."<<endl;
      break;
    case -4:
      os<<"NLOPT_ROUNDOFF_LIMITED: Halted because roundoff errors limited progress. (In this case, the optimization still typically returns a useful result."<<endl;
      break;
    case -5:
      os<<"NLOPT_FORCED_STOP: Halted because of a forced termination: the user called nlopt_force_stop(opt) on the optimization’s nlopt_opt object opt from the user’s objective function or constraints."<<endl;
      break;
    default:
      os<<"*** ERROR *** Unrecognized optimization code."<<endl;
      exit(1);
  }
  return os;
}


__END_AKANTU__

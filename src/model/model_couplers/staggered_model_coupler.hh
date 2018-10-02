/**
 * @file   staggered_model_coupler.hh
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 * 
 * @date creation: Fri Sep 28 2018
 * @date last modification: Fri Sep 28 2018
 *
 * @brief  class for staggered coupling of models 
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "model_coupler.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_STAGGERED_MODEL_COUPLER_HH__
#define __AKANTU_STAGGERED_MODEL_COUPLER_HH__

namespace akantu {

template <typename... Models>  
class StaggeredModelCoupler : public ModelCoupler<Models...>  {

  /* ------------------------------------------------------------------------ */
  /* Constructor/Destructors                                                  */
  /* ------------------------------------------------------------------------ */
public:

  StaggeredModelCoupler();

  virtual ~StaggeredModelCoupler();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  virtual void setCouplingDOFs();
  
  virtual void solve();

  /* ------------------------------------------------------------------------ */
  /* Members                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  
};


} // akantu


#endif

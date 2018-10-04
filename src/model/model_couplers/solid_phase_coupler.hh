/**
 * @file   solid_phase_coupler.hh
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 * 
 * @date creation: Fri Sep 28 2018
 * @date last modification: Fri Sep 28 2018
 *
 * @brief  class for coupling of solid mechancis and phasefield model 
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
#include "aka_common.hh"
#include "solid_mechanics_model.hh"
#include "phase_field_model.hh"
#include "material.hh"
#include "material_phasefield.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_SOLID_PHASE_COUPLER_HH__
#define __AKANTU_SOLID_PHASE_COUPLER_HH__

namespace akantu {

template <typename SolidType, typename PhaseType> 
class SolidPhaseCoupler {

  /* ------------------------------------------------------------------------ */
  /*  Constructor/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
private:
    
public:
  
  SolidPhaseCoupler(SolidType &, PhaseType &);

  ~SolidPhaseCoupler();
  

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// computes damage on quad points for solid mechanics model from
  /// damage array from phasefield model
  void computeDamageOnQuadPoints(const GhostType &);

  /// computes strain on quadrature points for phasefield model from
  /// displacement gradient from solid mechanics model
  void computeStrainOnQuadPoints(const GhostType & ghost_type);
  
  /// solve the coupled model
  void solve();

private:
  /// computes small strain from displacement gradient
  void gradUToEpsilon(const Matrix<Real> & grad_u,
		      Matrix<Real> & epsilon);

  /* ------------------------------------------------------------------------ */
  /* Members                                                                  */
  /* ------------------------------------------------------------------------ */
private:
  ///
  SolidType & solid;
  ///
  PhaseType & phase;
  
  /// Spatial dimension of models
  UInt spatial_dimension;
  
};
  
} // akantu 



#endif /* __AKANTU_SOLID_PHASE_COUPLER_HH__ */

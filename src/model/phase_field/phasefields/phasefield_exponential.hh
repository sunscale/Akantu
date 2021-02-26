/**
 * @file   phasefield_exponential.hh
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Fri Jun 18 2020
 * @date last modification: Mon Jan 29 2020
 *
 * @brief  Phasefield exponential
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
#include "phasefield.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_PHASEFIELD_EXPONENTIAL_HH__
#define __AKANTU_PHASEFIELD_EXPONENTIAL_HH__

namespace akantu {
class PhaseFieldExponential : public PhaseField {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  PhaseFieldExponential(PhaseFieldModel & model, const ID & id = "" );

  ~PhaseFieldExponential() override = default;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
protected:
  void computePhiOnQuad(const Matrix<Real> &, Real &, Real &);

  void computeDrivingForce(const ElementType & , GhostType ) override;
  
  inline void computeDrivingForceOnQuad(const Real &, Real &);

  inline void computeDamageEnergyDensityOnQuad(const Real &, Real &);


public:

  void updateInternalParameters() override;
};


/* -------------------------------------------------------------------------- */
inline void PhaseFieldExponential:: computeDrivingForceOnQuad(const Real & phi_quad,
							     Real & driving_force_quad){
  driving_force_quad = 2.0 * phi_quad;
}

/* -------------------------------------------------------------------------- */
inline void PhaseFieldExponential::computeDamageEnergyDensityOnQuad(const Real & phi_quad,
								    Real & dam_energy_quad) {
  dam_energy_quad = 2.0 * phi_quad + this->g_c/this->l0;
}
  
/* -------------------------------------------------------------------------- */
inline void PhaseFieldExponential::computePhiOnQuad(const Matrix<Real> & strain_quad,
						    Real & phi_quad, Real & phi_hist_quad)  {

  Matrix<Real> strain_plus(spatial_dimension, spatial_dimension);
  Matrix<Real> strain_minus(spatial_dimension, spatial_dimension);
  Matrix<Real> strain_dir(spatial_dimension, spatial_dimension);
  Matrix<Real> strain_diag_plus(spatial_dimension, spatial_dimension);
  Matrix<Real> strain_diag_minus(spatial_dimension, spatial_dimension);

  Vector<Real> strain_values(spatial_dimension);

  Real trace_plus, trace_minus;

  strain_plus.clear();
  strain_minus.clear();
  strain_dir.clear();
  strain_values.clear();
  strain_diag_plus.clear();
  strain_diag_minus.clear();
  
  strain_quad.eig(strain_values, strain_dir);

  for (UInt i = 0; i < spatial_dimension; i++) {
    strain_diag_plus(i, i) = std::max(Real(0.), strain_values(i));
    strain_diag_minus(i, i) = std::min(Real(0.), strain_values(i));
  }

  Matrix<Real> mat_tmp(spatial_dimension, spatial_dimension);
  Matrix<Real> sigma_plus(spatial_dimension, spatial_dimension);
  Matrix<Real> sigma_minus(spatial_dimension, spatial_dimension);
  
  mat_tmp.mul<false, true>(strain_diag_plus, strain_dir);
  strain_plus.mul<false, false>(strain_dir, mat_tmp);
  mat_tmp.mul<false, true>(strain_diag_minus, strain_dir);
  strain_minus.mul<false, true>(strain_dir, mat_tmp);
  
  trace_plus = std::max(Real(0.), strain_quad.trace());
  trace_minus = std::min(Real(0.), strain_quad.trace());
  
  for (UInt i = 0; i < spatial_dimension; i++) {
    for (UInt j = 0; j < spatial_dimension; j++) {
      sigma_plus(i, j) =
	(i == j) * lambda * trace_plus + 2 * mu * strain_plus(i, j);
      sigma_minus(i, j) =
	(i == j) * lambda * trace_minus + 2 * mu * strain_minus(i, j);
    }
  }
  
  phi_quad = 0.5 * sigma_plus.doubleDot(strain_quad);
  if (phi_quad < phi_hist_quad) 
    phi_quad = phi_hist_quad;
}
  
}

#endif

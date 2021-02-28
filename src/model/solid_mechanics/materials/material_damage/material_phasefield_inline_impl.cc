/**
 * @file   material_phasefield_inline_impl.cc
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Wed Oct 02 2018
 * @date last modification: Wed Oct 02 2018
 *
 * @brief  Implementation of the inline functions of the material phasefield
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

template<UInt spatial_dimension>
inline void MaterialPhaseField<spatial_dimension>::computeStressOnQuad(
    Matrix<Real> & grad_u, Matrix<Real> & sigma, Real & dam) {

  MaterialElastic<spatial_dimension>::computeStressOnQuad(grad_u, sigma);

  Matrix<Real> strain(           spatial_dimension, spatial_dimension);     
  Matrix<Real> strain_plus(      spatial_dimension, spatial_dimension);
  Matrix<Real> strain_minus(     spatial_dimension, spatial_dimension);
  Matrix<Real> strain_dir(       spatial_dimension, spatial_dimension);
  Matrix<Real> strain_diag_plus( spatial_dimension, spatial_dimension);
  Matrix<Real> strain_diag_minus(spatial_dimension, spatial_dimension);

  Vector<Real> strain_values(spatial_dimension);
  
  Real trace_plus, trace_minus;
  
  this->template gradUToEpsilon<spatial_dimension>(grad_u, strain);
  
  strain.eig(strain_values, strain_dir);

  for (UInt i=0; i < spatial_dimension; i++) {
    strain_diag_plus(i, i)  = std::max(Real(0.), strain_values(i));
    strain_diag_minus(i, i) = std::min(Real(0.), strain_values(i));
  }

  Matrix<Real> mat_tmp(    spatial_dimension, spatial_dimension);
  Matrix<Real> sigma_plus( spatial_dimension, spatial_dimension);
  Matrix<Real> sigma_minus(spatial_dimension, spatial_dimension);
  
  mat_tmp.mul<false,true>(strain_diag_plus, strain_dir);
  strain_plus.mul<false, false>(strain_dir, mat_tmp);
  mat_tmp.mul<false, true>(strain_diag_minus, strain_dir);
  strain_minus.mul<false, true>(strain_dir, mat_tmp);
  
  trace_plus  = std::max(Real(0.), strain.trace());
  trace_minus = std::min(Real(0.), strain.trace());

  Real lambda = MaterialElastic<spatial_dimension>::getLambda();
  Real mu     = MaterialElastic<spatial_dimension>::getMu();
  
  for (UInt i=0; i < spatial_dimension; i++) {
    for (UInt j=0; j < spatial_dimension; j++) {
      sigma_plus(i, j)  = (i==j) * lambda * trace_plus
	+ 2 * mu * strain_plus(i, j);
      sigma_minus(i, j) = (i==j) * lambda * trace_minus
	+ 2 * mu * strain_minus(i, j);
    }
  }     

  
  //sigma = (1 - dam) * sigma_plus + sigma_minus;
  sigma *= (1- dam)*(1-dam);
}

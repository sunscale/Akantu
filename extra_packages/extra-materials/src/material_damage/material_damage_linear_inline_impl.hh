/**
 * @file   material_damage_linear_inline_impl.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Marion Estelle Chambart <marion.chambart@epfl.ch>
 *
 *
 * @brief  Implementation of the inline functions of the material damage linear
 *
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
inline void MaterialDamageLinear<spatial_dimension>::computeStressOnQuad(
    Matrix<Real> & grad_u, Matrix<Real> & sigma, Real & dam, Real & K) {
  Real Fdiag[3];
  Real Fdiagp[3];

  Math::matrix33_eigenvalues(grad_u.storage(), Fdiag);

  Fdiagp[0] = std::max(0., Fdiag[0]);
  Fdiagp[1] = std::max(0., Fdiag[1]);
  Fdiagp[2] = std::max(0., Fdiag[2]);

  Real Ehat = sqrt(Fdiagp[0] * Fdiagp[0] + Fdiagp[1] * Fdiagp[1] +
                   Fdiagp[2] * Fdiagp[2]);

  MaterialElastic<spatial_dimension>::computeStressOnQuad(grad_u, sigma);

  Real Fd = Ehat - K;

  if (Fd > 0) {
    dam = (Ehat - Epsmin) / (Epsmax - Epsmin) * (Ehat / Epsmax);
    dam = std::min(dam, 1.);
    K = Ehat;
  }

  sigma *= 1 - dam;
}

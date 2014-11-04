/**
 * @file   material_damage_linear_inline_impl.cc
 *
 * @author Marion Estelle Chambart <marion.chambart@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date   Mon Oct 03 11:53:09 2011
 *
 * @brief  Implementation of the inline functions of the material damage linear
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
template<UInt spatial_dimension>
inline void
MaterialDamageLinear<spatial_dimension>::computeStressOnQuad(Matrix<Real> & grad_u,
							     Matrix<Real> & sigma,
							     Real & dam,
							     Real & K) {
  Real Fdiag[3];
  Real Fdiagp[3];

  Math::matrix33_eigenvalues(grad_u.storage(), Fdiag);

  Fdiagp[0] = std::max(0., Fdiag[0]);
  Fdiagp[1] = std::max(0., Fdiag[1]);
  Fdiagp[2] = std::max(0., Fdiag[2]);

  Real Ehat=sqrt(Fdiagp[0]*Fdiagp[0] + Fdiagp[1]*Fdiagp[1] + Fdiagp[2]*Fdiagp[2]);

  MaterialElastic<spatial_dimension>::computeStressOnQuad(grad_u, sigma);

  Real Fd = Ehat-K;

  if (Fd > 0) {
    dam = (Ehat - Epsmin) / (Epsmax-Epsmin)*(Ehat/Epsmax);
    dam = std::min(dam, 1.);
    K=Ehat;
  }

  sigma *= 1-dam;
}


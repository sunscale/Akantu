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
  sigma *= 1- dam;
}

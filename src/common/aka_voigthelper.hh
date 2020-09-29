/**
 * @file   aka_voigthelper.hh
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 * @author Till Junge <till.junge@epfl.ch>
 * @author Daniel Pino Muñoz <daniel.pinomunoz@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Dec 20 2013
 * @date last modification: Mon Jan 29 2018
 *
 * @brief  Helper file for Voigt notation
 * Wikipedia convention: @f[2*\epsilon_{ij} (i!=j) = voigt_\epsilon_{I}@f]
 * http://en.wikipedia.org/wiki/Voigt_notation
 *
 *
 * Copyright (©) 2014-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "aka_types.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKA_VOIGTHELPER_HH_
#define AKA_VOIGTHELPER_HH_

namespace akantu {

/* -------------------------------------------------------------------------- */
template <UInt dim> class VoigtHelper {
  static_assert(dim > 0U, "Cannot be < 1D");
  static_assert(dim < 4U, "Cannot be > 3D");

public:
  /* ------------------------------------------------------------------------ */
  template <class M, class V>
  static inline void matrixToVoigt(M && matrix, V && vector);

  template <class M> static inline decltype(auto) matrixToVoigt(M && matrix);

  template <class M, class V>
  static inline void matrixToVoigtWithFactors(M && matrix, V && vector);

  template <class M>
  static inline decltype(auto) matrixToVoigtWithFactors(M && matrix);

  template <class M, class V>
  static inline void voigtToMatrix(V && vector, M && matrix);

  template <class V> static inline decltype(auto) voigtToMatrix(V && vector);
  /* ------------------------------------------------------------------------ */

  /// transfer the B matrix to a Voigt notation B matrix
  inline static void transferBMatrixToSymVoigtBMatrix(
      const Matrix<Real> & B, Matrix<Real> & Bvoigt, UInt nb_nodes_per_element);

  /// transfer the BNL matrix to a Voigt notation B matrix (See Bathe et al.
  /// IJNME vol 9, 1975)
  inline static void transferBMatrixToBNL(const Matrix<Real> & B,
                                          Matrix<Real> & Bvoigt,
                                          UInt nb_nodes_per_element);

  /// transfer the BL2 matrix to a Voigt notation B matrix (See Bathe et al.
  /// IJNME vol 9, 1975)
  inline static void transferBMatrixToBL2(const Matrix<Real> & B,
                                          const Matrix<Real> & grad_u,
                                          Matrix<Real> & Bvoigt,
                                          UInt nb_nodes_per_element);

public:
  static constexpr UInt size{(dim * (dim - 1)) / 2 + dim};
  // matrix of vector index I as function of tensor indices i,j
  static const UInt mat[dim][dim];
  // array of matrix indices ij as function of vector index I
  static const UInt vec[dim * dim][2];
  // factors to multiply the strain by for voigt notation
  static const Real factors[size];
};

} // namespace akantu

#include "aka_voigthelper_tmpl.hh"

#endif

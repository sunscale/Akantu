/**
 * @file   aka_voigthelper_tmpl.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Dec 20 2013
 * @date last modification: Wed Dec 06 2017
 *
 * @brief  implementation of the voight helper
 *
 * @section LICENSE
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

/**
 * @file   aka_voigthelper_tmpl.hh
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @date   Wed Nov 16 12:22:58 2016
 */

/* -------------------------------------------------------------------------- */
#include "aka_voigthelper.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_AKA_VOIGTHELPER_TMPL_HH__
#define __AKANTU_AKA_VOIGTHELPER_TMPL_HH__

namespace akantu {

template <UInt dim> constexpr UInt VoigtHelper<dim>::size;

/* -------------------------------------------------------------------------- */
template <UInt dim>
inline void VoigtHelper<dim>::transferBMatrixToSymVoigtBMatrix(
    const Matrix<Real> & B, Matrix<Real> & Bvoigt, UInt nb_nodes_per_element) {
  Bvoigt.clear();

  for (UInt i = 0; i < dim; ++i)
    for (UInt n = 0; n < nb_nodes_per_element; ++n)
      Bvoigt(i, i + n * dim) = B(i, n);

  if (dim == 2) {
    /// in 2D, fill the @f$ [\frac{\partial N_i}{\partial x}, \frac{\partial
    /// N_i}{\partial y}]@f$ row
    for (UInt n = 0; n < nb_nodes_per_element; ++n) {
      Bvoigt(2, 1 + n * 2) = B(0, n);
      Bvoigt(2, 0 + n * 2) = B(1, n);
    }
  }

  if (dim == 3) {
    for (UInt n = 0; n < nb_nodes_per_element; ++n) {
      Real dndx = B(0, n);
      Real dndy = B(1, n);
      Real dndz = B(2, n);

      /// in 3D, fill the @f$ [0, \frac{\partial N_i}{\partial y},
      /// \frac{N_i}{\partial z}]@f$ row
      Bvoigt(3, 1 + n * 3) = dndz;
      Bvoigt(3, 2 + n * 3) = dndy;

      /// in 3D, fill the @f$ [\frac{\partial N_i}{\partial x}, 0,
      /// \frac{N_i}{\partial z}]@f$ row
      Bvoigt(4, 0 + n * 3) = dndz;
      Bvoigt(4, 2 + n * 3) = dndx;

      /// in 3D, fill the @f$ [\frac{\partial N_i}{\partial x},
      /// \frac{N_i}{\partial y}, 0]@f$ row
      Bvoigt(5, 0 + n * 3) = dndy;
      Bvoigt(5, 1 + n * 3) = dndx;
    }
  }
}

/* -------------------------------------------------------------------------- */
template <UInt dim>
inline void VoigtHelper<dim>::transferBMatrixToBNL(const Matrix<Real> & B,
                                                   Matrix<Real> & Bvoigt,
                                                   UInt nb_nodes_per_element) {
  Bvoigt.clear();

  // see Finite element formulations for large deformation dynamic analysis,
  // Bathe et al. IJNME vol 9, 1975, page 364 B_{NL}
  for (UInt i = 0; i < dim; ++i) {
    for (UInt m = 0; m < nb_nodes_per_element; ++m) {
      for (UInt n = 0; n < dim; ++n) {
        // std::cout << B(n, m) << std::endl;
        Bvoigt(i * dim + n, m * dim + i) = B(n, m);
      }
    }
  }
  // TODO: Verify the 2D and 1D case
}

/* -------------------------------------------------------------------------- */
template <>
inline void VoigtHelper<1>::transferBMatrixToBL2(const Matrix<Real> & B,
                                                 const Matrix<Real> & grad_u,
                                                 Matrix<Real> & Bvoigt,
                                                 UInt nb_nodes_per_element) {
  Bvoigt.clear();
  for (UInt j = 0; j < nb_nodes_per_element; ++j)
    Bvoigt(0, j) = grad_u(0, 0) * B(0, j);
}

/* -------------------------------------------------------------------------- */
template <>
inline void VoigtHelper<3>::transferBMatrixToBL2(const Matrix<Real> & B,
                                                 const Matrix<Real> & grad_u,
                                                 Matrix<Real> & Bvoigt,
                                                 UInt nb_nodes_per_element) {
  Bvoigt.clear();

  for (UInt i = 0; i < 3; ++i)
    for (UInt j = 0; j < nb_nodes_per_element; ++j)
      for (UInt k = 0; k < 3; ++k)
        Bvoigt(i, j * 3 + k) = grad_u(k, i) * B(i, j);

  for (UInt i = 3; i < 6; ++i) {
    for (UInt j = 0; j < nb_nodes_per_element; ++j) {
      for (UInt k = 0; k < 3; ++k) {
        UInt aux = i - 3;
        for (UInt m = 0; m < 3; ++m) {
          if (m != aux) {
            UInt index1 = m;
            UInt index2 = 3 - m - aux;
            Bvoigt(i, j * 3 + k) += grad_u(k, index1) * B(index2, j);
          }
        }
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
template <>
inline void VoigtHelper<2>::transferBMatrixToBL2(const Matrix<Real> & B,
                                                 const Matrix<Real> & grad_u,
                                                 Matrix<Real> & Bvoigt,
                                                 UInt nb_nodes_per_element) {

  Bvoigt.clear();

  for (UInt i = 0; i < 2; ++i)
    for (UInt j = 0; j < nb_nodes_per_element; ++j)
      for (UInt k = 0; k < 2; ++k)
        Bvoigt(i, j * 2 + k) = grad_u(k, i) * B(i, j);

  for (UInt j = 0; j < nb_nodes_per_element; ++j) {
    for (UInt k = 0; k < 2; ++k) {
      for (UInt m = 0; m < 2; ++m) {
        UInt index1 = m;
        UInt index2 = (2 - 1) - m;
        Bvoigt(2, j * 2 + k) += grad_u(k, index1) * B(index2, j);
      }
    }
  }
}

} // namespace akantu

#endif /* __AKANTU_AKA_VOIGTHELPER_TMPL_HH__ */

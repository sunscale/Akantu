/**
 * @file   aka_math.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Marion Estelle Chambart <marion.chambart@epfl.ch>
 * @author David Simon Kammer <david.kammer@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Leonardo Snozzi <leonardo.snozzi@epfl.ch>
 * @author Peter Spijker <peter.spijker@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Wed Aug 04 2010
 * @date last modification: Sun Aug 13 2017
 *
 * @brief  Implementation of the math toolbox
 *
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
#include "aka_math.hh"
#include "aka_array.hh"
#include "aka_iterators.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {
namespace Math {
  /* --------------------------------------------------------------------------
   */
  void matrix_vector(UInt m, UInt n, const Array<Real> & A,
                     const Array<Real> & x, Array<Real> & y, Real alpha) {
    AKANTU_DEBUG_IN();

    AKANTU_DEBUG_ASSERT(A.size() == x.size(),
                        "The vector A(" << A.getID() << ") and the vector x("
                                        << x.getID()
                                        << ") must have the same size");

    AKANTU_DEBUG_ASSERT(
        A.getNbComponent() == m * n,
        "The vector A(" << A.getID() << ") has the good number of component.");

    AKANTU_DEBUG_ASSERT(x.getNbComponent() == n,
                        "The vector x("
                            << x.getID()
                            << ") do not the good number of component.");

    AKANTU_DEBUG_ASSERT(y.getNbComponent() == n,
                        "The vector y("
                            << y.getID()
                            << ") do not the good number of component.");

    UInt nb_element = A.size();
    UInt offset_A = A.getNbComponent();
    UInt offset_x = x.getNbComponent();

    y.resize(nb_element);

    Real * A_val = A.storage();
    Real * x_val = x.storage();
    Real * y_val = y.storage();

    for (UInt el = 0; el < nb_element; ++el) {
      matrix_vector(m, n, A_val, x_val, y_val, alpha);

      A_val += offset_A;
      x_val += offset_x;
      y_val += offset_x;
    }

    AKANTU_DEBUG_OUT();
  }

  /* --------------------------------------------------------------------------
   */
  void matrix_matrix(UInt m, UInt n, UInt k, const Array<Real> & A,
                     const Array<Real> & B, Array<Real> & C, Real alpha) {
    AKANTU_DEBUG_IN();

    AKANTU_DEBUG_ASSERT(A.size() == B.size(),
                        "The vector A(" << A.getID() << ") and the vector B("
                                        << B.getID()
                                        << ") must have the same size");

    AKANTU_DEBUG_ASSERT(
        A.getNbComponent() == m * k,
        "The vector A(" << A.getID() << ") has the good number of component.");

    AKANTU_DEBUG_ASSERT(B.getNbComponent() == k * n,
                        "The vector B("
                            << B.getID()
                            << ") do not the good number of component.");

    AKANTU_DEBUG_ASSERT(C.getNbComponent() == m * n,
                        "The vector C("
                            << C.getID()
                            << ") do not the good number of component.");

    UInt nb_element = A.size();
    UInt offset_A = A.getNbComponent();
    UInt offset_B = B.getNbComponent();
    UInt offset_C = C.getNbComponent();

    C.resize(nb_element);

    Real * A_val = A.storage();
    Real * B_val = B.storage();
    Real * C_val = C.storage();

    for (UInt el = 0; el < nb_element; ++el) {
      matrix_matrix(m, n, k, A_val, B_val, C_val, alpha);

      A_val += offset_A;
      B_val += offset_B;
      C_val += offset_C;
    }

    AKANTU_DEBUG_OUT();
  }

  /* --------------------------------------------------------------------------
   */
  void matrix_matrixt(UInt m, UInt n, UInt k, const Array<Real> & A,
                      const Array<Real> & B, Array<Real> & C, Real alpha) {
    AKANTU_DEBUG_IN();

    AKANTU_DEBUG_ASSERT(A.size() == B.size(),
                        "The vector A(" << A.getID() << ") and the vector B("
                                        << B.getID()
                                        << ") must have the same size");

    AKANTU_DEBUG_ASSERT(
        A.getNbComponent() == m * k,
        "The vector A(" << A.getID() << ") has the good number of component.");

    AKANTU_DEBUG_ASSERT(B.getNbComponent() == k * n,
                        "The vector B("
                            << B.getID()
                            << ") do not the good number of component.");

    AKANTU_DEBUG_ASSERT(C.getNbComponent() == m * n,
                        "The vector C("
                            << C.getID()
                            << ") do not the good number of component.");

    UInt nb_element = A.size();
    UInt offset_A = A.getNbComponent();
    UInt offset_B = B.getNbComponent();
    UInt offset_C = C.getNbComponent();

    C.resize(nb_element);

    Real * A_val = A.storage();
    Real * B_val = B.storage();
    Real * C_val = C.storage();

    for (UInt el = 0; el < nb_element; ++el) {
      matrix_matrixt(m, n, k, A_val, B_val, C_val, alpha);

      A_val += offset_A;
      B_val += offset_B;
      C_val += offset_C;
    }

    AKANTU_DEBUG_OUT();
  }

  /* --------------------------------------------------------------------------
   */
  void compute_tangents(const Array<Real> & normals, Array<Real> & tangents) {
    AKANTU_DEBUG_IN();

    if (normals.empty()) {
      return;
    }

    auto spatial_dimension = normals.getNbComponent();
    auto tangent_components = spatial_dimension * (spatial_dimension - 1);

    if (tangent_components == 0) {
      return;
    }

    AKANTU_DEBUG_ASSERT(
        tangent_components == tangents.getNbComponent(),
        "Cannot compute the tangents, the storage array for tangents"
            << " does not have the good amount of components.");

    auto nb_normals = normals.size();
    tangents.resize(nb_normals);
    tangents.zero();

    /// compute first tangent
    for (auto && data : zip(make_view(normals, spatial_dimension),
                            make_view(tangents, tangent_components))) {
      const auto & normal = std::get<0>(data);
      auto & tangent = std::get<1>(data);

      if (are_float_equal(norm2(normal.storage()), 0.)) {
        tangent(0) = 1.;
      } else {
        normal2(normal.storage(), tangent.storage());
      }
    }

    /// compute second tangent (3D case)
    if (spatial_dimension == 3) {
      for (auto && data : zip(make_view(normals, spatial_dimension),
                              make_view(tangents, tangent_components))) {
        const auto & normal = std::get<0>(data);
        auto & tangent = std::get<1>(data);

        normal3(normal.storage(), tangent.storage(),
                tangent.storage() + spatial_dimension);
      }
    }

    AKANTU_DEBUG_OUT();
  } // namespace akantu

  /* --------------------------------------------------------------------------
   */
  Real reduce(Array<Real> & array) {
    UInt nb_values = array.size();
    if (nb_values == 0) {
      return 0.;
    }

    UInt nb_values_to_sum = nb_values >> 1;

    std::sort(array.begin(), array.end());

    // as long as the half is not empty
    while (nb_values_to_sum != 0U) {
      UInt remaining = (nb_values - 2 * nb_values_to_sum);
      if (remaining != 0U) {
        array(nb_values - 2) += array(nb_values - 1);
      }

      // sum to consecutive values and store the sum in the first half
      for (UInt i = 0; i < nb_values_to_sum; ++i) {
        array(i) = array(2 * i) + array(2 * i + 1);
      }

      nb_values = nb_values_to_sum;
      nb_values_to_sum >>= 1;
    }

    return array(0);
  }

} // namespace Math
} // namespace akantu

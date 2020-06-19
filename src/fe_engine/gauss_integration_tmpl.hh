/**
 * @file   gauss_integration_tmpl.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Tue May 10 2016
 * @date last modification: Wed Nov 29 2017
 *
 * @brief  implementation of the gauss integration helpers
 *
 *
 * Copyright (©) 2016-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

#ifndef __AKANTU_GAUSS_INTEGRATION_TMPL_HH__
#define __AKANTU_GAUSS_INTEGRATION_TMPL_HH__

namespace akantu {

/* -------------------------------------------------------------------------- */
/* GaussIntegrationElement                                                    */
/* -------------------------------------------------------------------------- */
namespace _aka_gauss_helpers {
  template <GaussIntegrationType type, UInt n>
  struct GaussIntegrationNbPoints {
    static const UInt nb_points = 0;
  };

#if !defined(DOXYGEN)
  template <UInt n> struct GaussIntegrationNbPoints<_git_not_defined, n> {
    static const UInt nb_points = 0;
  };

  template <UInt n> struct GaussIntegrationNbPoints<_git_point, n> {
    static const UInt nb_points = 1;
  };

  template <UInt n> struct GaussIntegrationNbPoints<_git_segment, n> {
    static const UInt nb_points = (n + 1) / 2 + ((n + 1) % 2 ? 1 : 0);
  };

#define DECLARE_GAUSS_NB_POINTS(type, order, points)                           \
  template <> struct GaussIntegrationNbPoints<type, order> {                   \
    static const UInt nb_points = points;                                      \
  }

#define DECLARE_GAUSS_NB_POINTS_PENT(type, order, xo, yo)                      \
  template <> struct GaussIntegrationNbPoints<type, order> {                   \
    static const UInt x_order = xo;                                            \
    static const UInt yz_order = yo;                                           \
    static const UInt nb_points = 1;                                           \
  }

  DECLARE_GAUSS_NB_POINTS(_git_triangle, 1, 1);
  DECLARE_GAUSS_NB_POINTS(_git_triangle, 2, 3);
  DECLARE_GAUSS_NB_POINTS(_git_triangle, 3, 4);
  DECLARE_GAUSS_NB_POINTS(_git_triangle, 4, 6);
  DECLARE_GAUSS_NB_POINTS(_git_triangle, 5, 7);
  DECLARE_GAUSS_NB_POINTS(_git_tetrahedron, 1, 1);
  DECLARE_GAUSS_NB_POINTS(_git_tetrahedron, 2, 4);
  DECLARE_GAUSS_NB_POINTS(_git_tetrahedron, 3, 5);
  DECLARE_GAUSS_NB_POINTS(_git_tetrahedron, 4, 15);
  DECLARE_GAUSS_NB_POINTS(_git_tetrahedron, 5, 15);
  DECLARE_GAUSS_NB_POINTS_PENT(_git_pentahedron, 1, 3,
                               2); // order 3 in x, order 2 in y and z
  DECLARE_GAUSS_NB_POINTS_PENT(_git_pentahedron, 2, 3,
                               2); // order 3 in x, order 2 in y and z
  DECLARE_GAUSS_NB_POINTS_PENT(_git_pentahedron, 3, 3,
                               3); // order 3 in x, order 3 in y and z
  DECLARE_GAUSS_NB_POINTS_PENT(_git_pentahedron, 4, 5,
                               5); // order 5 in x, order 5 in y and z
  DECLARE_GAUSS_NB_POINTS_PENT(_git_pentahedron, 5, 5,
                               5); // order 5 in x, order 5 in y and z

  template <GaussIntegrationType type, UInt n, UInt on = n,
            bool end_recurse = false>
  struct GaussIntegrationNbPointsHelper {
    static const UInt pnp = GaussIntegrationNbPoints<type, n>::nb_points;
    static const UInt order = n;
    static const UInt nb_points = pnp;
  };

  template <GaussIntegrationType type, UInt n, UInt on>
  struct GaussIntegrationNbPointsHelper<type, n, on, true> {
    static const UInt nb_points = 0;
  };
#endif

  /* ------------------------------------------------------------------------ */
  /* Generic helper                                                           */
  /* ------------------------------------------------------------------------ */
  template <GaussIntegrationType type, UInt dimension, UInt n>
  struct GaussIntegrationTypeDataHelper {
    using git_np = GaussIntegrationNbPoints<type, n>;
    using git_data = GaussIntegrationTypeData<type, git_np::nb_points>;

    static UInt getNbQuadraturePoints() { return git_np::nb_points; }

    static const Matrix<Real> getQuadraturePoints() {
      return Matrix<Real>(git_data::quad_positions, dimension,
                          git_np::nb_points);
    }

    static const Vector<Real> getWeights() {
      return Vector<Real>(git_data::quad_weights, git_np::nb_points);
    }
  };

#if !defined(DOXYGEN)
  /* ------------------------------------------------------------------------ */
  /* helper for _segment _quadrangle _hexahedron                              */
  /* ------------------------------------------------------------------------ */
  template <UInt dimension, UInt dp>
  struct GaussIntegrationTypeDataHelper<_git_segment, dimension, dp> {
    using git_np = GaussIntegrationNbPoints<_git_segment, dp>;
    using git_data = GaussIntegrationTypeData<_git_segment, git_np::nb_points>;

    static UInt getNbQuadraturePoints() {
      return Math::pow<dimension>(git_np::nb_points);
    }

    static const Matrix<Real> getQuadraturePoints() {
      UInt tot_nquad = getNbQuadraturePoints();
      UInt nquad = git_np::nb_points;

      Matrix<Real> quads(dimension, tot_nquad);
      Vector<Real> pos(git_data::quad_positions, nquad);

      UInt offset = 1;
      for (UInt d = 0; d < dimension; ++d) {
        for (UInt n = 0, q = 0; n < tot_nquad; ++n, q += offset) {
          UInt rq = q % tot_nquad + q / tot_nquad;
          quads(d, rq) = pos(n % nquad);
        }
        offset *= nquad;
      }
      return quads;
    }

    static const Vector<Real> getWeights() {
      UInt tot_nquad = getNbQuadraturePoints();
      UInt nquad = git_np::nb_points;

      Vector<Real> quads_weights(tot_nquad, 1.);
      Vector<Real> weights(git_data::quad_weights, nquad);

      UInt offset = 1;
      for (UInt d = 0; d < dimension; ++d) {
        for (UInt n = 0, q = 0; n < tot_nquad; ++n, q += offset) {
          UInt rq = q % tot_nquad + q / tot_nquad;
          quads_weights(rq) *= weights(n % nquad);
        }
        offset *= nquad;
      }
      return quads_weights;
    }
  };

  /* ------------------------------------------------------------------------ */
  /* helper for _pentahedron                                                  */
  /* ------------------------------------------------------------------------ */
  template <UInt dimension, UInt dp>
  struct GaussIntegrationTypeDataHelper<_git_pentahedron, dimension, dp> {
    using git_info = GaussIntegrationNbPoints<_git_pentahedron, dp>;
    using git_np_seg =
        GaussIntegrationNbPoints<_git_segment, git_info::x_order>;
    using git_np_tri =
        GaussIntegrationNbPoints<_git_triangle, git_info::yz_order>;
    using git_data_seg =
        GaussIntegrationTypeData<_git_segment, git_np_seg::nb_points>;
    using git_data_tri =
        GaussIntegrationTypeData<_git_triangle, git_np_tri::nb_points>;

    static UInt getNbQuadraturePoints() {
      return git_np_seg::nb_points * git_np_tri::nb_points;
    }

    static const Matrix<Real> getQuadraturePoints() {
      UInt tot_nquad = getNbQuadraturePoints();
      UInt nquad_seg = git_np_seg::nb_points;
      UInt nquad_tri = git_np_tri::nb_points;

      Matrix<Real> quads(dimension, tot_nquad);
      Matrix<Real> pos_seg_w(git_data_seg::quad_positions, 1, nquad_seg);
      Matrix<Real> pos_tri_w(git_data_tri::quad_positions, 2, nquad_tri);

      for (UInt ns = 0, q = 0; ns < nquad_seg; ++ns) {
        Vector<Real> pos_seg = pos_seg_w(ns);
        for (UInt nt = 0; nt < nquad_tri; ++nt, ++q) {
          Vector<Real> pos_tri = pos_tri_w(nt);
          Vector<Real> quad = quads(q);
          quad(_x) = pos_seg(_x);
          quad(_y) = pos_tri(_x);
          quad(_z) = pos_tri(_y);
        }
      }
      return quads;
    }

    static const Vector<Real> getWeights() {
      UInt tot_nquad = getNbQuadraturePoints();
      UInt nquad_seg = git_np_seg::nb_points;
      UInt nquad_tri = git_np_tri::nb_points;

      Vector<Real> quads_weights(tot_nquad);

      Vector<Real> weight_seg(git_data_seg::quad_weights, nquad_seg);
      Vector<Real> weight_tri(git_data_tri::quad_weights, nquad_tri);

      for (UInt ns = 0, q = 0; ns < nquad_seg; ++ns) {
        for (UInt nt = 0; nt < nquad_tri; ++nt, ++q) {
          quads_weights(q) = weight_seg(ns) * weight_tri(nt);
        }
      }
      return quads_weights;
    }
  };
#endif
} // namespace _aka_gauss_helpers

template <ElementType element_type, UInt n>
const Matrix<Real>
GaussIntegrationElement<element_type, n>::getQuadraturePoints() {
  const InterpolationType itp_type =
      ElementClassProperty<element_type>::interpolation_type;
  using interpolation_property = InterpolationProperty<itp_type>;
  using data_helper = _aka_gauss_helpers::GaussIntegrationTypeDataHelper<
      ElementClassProperty<element_type>::gauss_integration_type,
      interpolation_property::natural_space_dimension, n>;
  Matrix<Real> tmp(data_helper::getQuadraturePoints());
  return tmp;
}

/* -------------------------------------------------------------------------- */
template <ElementType element_type, UInt n>
const Vector<Real> GaussIntegrationElement<element_type, n>::getWeights() {
  const InterpolationType itp_type =
      ElementClassProperty<element_type>::interpolation_type;
  using interpolation_property = InterpolationProperty<itp_type>;
  using data_helper = _aka_gauss_helpers::GaussIntegrationTypeDataHelper<
      ElementClassProperty<element_type>::gauss_integration_type,
      interpolation_property::natural_space_dimension, n>;
  Vector<Real> tmp(data_helper::getWeights());
  return tmp;
}

/* -------------------------------------------------------------------------- */
template <ElementType element_type, UInt n>
UInt GaussIntegrationElement<element_type, n>::getNbQuadraturePoints() {
  const InterpolationType itp_type =
      ElementClassProperty<element_type>::interpolation_type;
  using interpolation_property = InterpolationProperty<itp_type>;
  using data_helper = _aka_gauss_helpers::GaussIntegrationTypeDataHelper<
      ElementClassProperty<element_type>::gauss_integration_type,
      interpolation_property::natural_space_dimension, n>;
  return data_helper::getNbQuadraturePoints();
}

} // namespace akantu

#endif /* __AKANTU_GAUSS_INTEGRATION_TMPL_HH__ */

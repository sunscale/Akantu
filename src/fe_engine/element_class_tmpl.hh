/**
 * @file   element_class_tmpl.hh
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Thomas Menouillard <tmenouillard@stucky.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Thu Feb 21 2013
 * @date last modification: Wed Nov 29 2017
 *
 * @brief  Implementation of the inline templated function of the element class
 * descriptions
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
#include "element_class.hh"
#include "gauss_integration_tmpl.hh"
/* -------------------------------------------------------------------------- */
#include <type_traits>
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_ELEMENT_CLASS_TMPL_HH_
#define AKANTU_ELEMENT_CLASS_TMPL_HH_

namespace akantu {

template <ElementType element_type, ElementKind element_kind>
inline constexpr auto
ElementClass<element_type, element_kind>::getFacetTypes() {
  return VectorProxy<const ElementType>(
      element_class_extra_geom_property::facet_type.data(),
      geometrical_element::getNbFacetTypes());
}

/* -------------------------------------------------------------------------- */
/* GeometricalElement                                                         */
/* -------------------------------------------------------------------------- */
template <GeometricalType geometrical_type, GeometricalShapeType shape>
inline constexpr auto
GeometricalElement<geometrical_type,
                   shape>::getFacetLocalConnectivityPerElement(UInt t) {
  int pos = 0;
  for (UInt i = 0; i < t; ++i) {
    pos += geometrical_property::nb_facets[i] *
           geometrical_property::nb_nodes_per_facet[i];
  }

  return MatrixProxy<const UInt>(
      geometrical_property::facet_connectivity_vect.data() + pos,
      geometrical_property::nb_facets[t],
      geometrical_property::nb_nodes_per_facet[t]);
}

/* -------------------------------------------------------------------------- */
template <GeometricalType geometrical_type, GeometricalShapeType shape>
inline UInt
GeometricalElement<geometrical_type, shape>::getNbFacetsPerElement() {
  UInt total_nb_facets = 0;
  for (UInt n = 0; n < geometrical_property::nb_facet_types; ++n) {
    total_nb_facets += geometrical_property::nb_facets[n];
  }

  return total_nb_facets;
}

/* -------------------------------------------------------------------------- */
template <GeometricalType geometrical_type, GeometricalShapeType shape>
inline UInt
GeometricalElement<geometrical_type, shape>::getNbFacetsPerElement(UInt t) {
  return geometrical_property::nb_facets[t];
}

/* -------------------------------------------------------------------------- */
template <GeometricalType geometrical_type, GeometricalShapeType shape>
template <class vector_type>
inline bool GeometricalElement<geometrical_type, shape>::contains(
    const vector_type & coords) {
  return GeometricalShapeContains<shape>::contains(coords);
}

/* -------------------------------------------------------------------------- */
template <>
template <class vector_type>
inline bool
GeometricalShapeContains<_gst_point>::contains(const vector_type & coords) {
  return (coords(0) < std::numeric_limits<Real>::epsilon());
}

/* -------------------------------------------------------------------------- */
template <>
template <class vector_type>
inline bool
GeometricalShapeContains<_gst_square>::contains(const vector_type & coords) {
  bool in = true;
  for (UInt i = 0; i < coords.size() && in; ++i) {
    in &= ((coords(i) >= -(1. + std::numeric_limits<Real>::epsilon())) &&
           (coords(i) <= (1. + std::numeric_limits<Real>::epsilon())));
  }
  return in;
}

/* -------------------------------------------------------------------------- */
template <>
template <class vector_type>
inline bool
GeometricalShapeContains<_gst_triangle>::contains(const vector_type & coords) {
  bool in = true;
  Real sum = 0;
  for (UInt i = 0; (i < coords.size()) && in; ++i) {
    in &= ((coords(i) >= -(Math::getTolerance())) &&
           (coords(i) <= (1. + Math::getTolerance())));
    sum += coords(i);
  }
  if (in) {
    return (in && (sum <= (1. + Math::getTolerance())));
  }
  return in;
}

/* -------------------------------------------------------------------------- */
template <>
template <class vector_type>
inline bool
GeometricalShapeContains<_gst_prism>::contains(const vector_type & coords) {
  bool in = ((coords(0) >= -1.) && (coords(0) <= 1.)); // x in segment [-1, 1]

  // y and z in triangle
  in &= ((coords(1) >= 0) && (coords(1) <= 1.));
  in &= ((coords(2) >= 0) && (coords(2) <= 1.));
  Real sum = coords(1) + coords(2);

  return (in && (sum <= 1));
}

/* -------------------------------------------------------------------------- */
/* InterpolationElement                                                       */
/* -------------------------------------------------------------------------- */
template <InterpolationType interpolation_type, InterpolationKind kind>
inline void InterpolationElement<interpolation_type, kind>::computeShapes(
    const Matrix<Real> & natural_coord, Matrix<Real> & N) {
  UInt nb_points = natural_coord.cols();
  for (UInt p = 0; p < nb_points; ++p) {
    Vector<Real> Np(N(p));
    Vector<Real> ncoord_p(natural_coord(p));
    computeShapes(ncoord_p, Np);
  }
}

/* -------------------------------------------------------------------------- */
template <InterpolationType interpolation_type, InterpolationKind kind>
inline void InterpolationElement<interpolation_type, kind>::computeDNDS(
    const Matrix<Real> & natural_coord, Tensor3<Real> & dnds) {
  UInt nb_points = natural_coord.cols();
  for (UInt p = 0; p < nb_points; ++p) {
    Matrix<Real> dnds_p(dnds(p));
    Vector<Real> ncoord_p(natural_coord(p));
    computeDNDS(ncoord_p, dnds_p);
  }
}

/* -------------------------------------------------------------------------- */
/**
 * interpolate on a point a field for which values are given on the
 * node of the element using the shape functions at this interpolation point
 *
 * @param nodal_values values of the function per node @f$ f_{ij} = f_{n_i j}
 *@f$ so it should be a matrix of size nb_nodes_per_element @f$\times@f$
 *nb_degree_of_freedom
 * @param shapes value of shape functions at the interpolation point
 * @param interpolated interpolated value of f @f$ f_j(\xi) = \sum_i f_{n_i j}
 *N_i @f$
 */
template <InterpolationType interpolation_type, InterpolationKind kind>
inline void InterpolationElement<interpolation_type, kind>::interpolate(
    const Matrix<Real> & nodal_values, const Vector<Real> & shapes,
    Vector<Real> & interpolated) {
  Matrix<Real> interpm(interpolated.storage(), nodal_values.rows(), 1);
  Matrix<Real> shapesm(
      shapes.storage(),
      InterpolationProperty<interpolation_type>::nb_nodes_per_element, 1);
  interpm.mul<false, false>(nodal_values, shapesm);
}

/* -------------------------------------------------------------------------- */
/**
 * interpolate on several points a field  for which values are given on the
 * node of the element using the shape functions at the interpolation point
 *
 * @param nodal_values values of the function per node @f$ f_{ij} = f_{n_i j}
 *@f$ so it should be a matrix of size nb_nodes_per_element @f$\times@f$
 *nb_degree_of_freedom
 * @param shapes value of shape functions at the interpolation point
 * @param interpolated interpolated values of f @f$ f_j(\xi) = \sum_i f_{n_i j}
 *N_i @f$
 */
template <InterpolationType interpolation_type, InterpolationKind kind>
inline void InterpolationElement<interpolation_type, kind>::interpolate(
    const Matrix<Real> & nodal_values, const Matrix<Real> & shapes,
    Matrix<Real> & interpolated) {
  UInt nb_points = shapes.cols();
  for (UInt p = 0; p < nb_points; ++p) {
    Vector<Real> Np(shapes(p));
    Vector<Real> interpolated_p(interpolated(p));
    interpolate(nodal_values, Np, interpolated_p);
  }
}

/* -------------------------------------------------------------------------- */
/**
 * interpolate the field on a point given in natural coordinates the field which
 * values are given on the node of the element
 *
 * @param natural_coords natural coordinates of point where to interpolate \xi
 * @param nodal_values values of the function per node @f$ f_{ij} = f_{n_i j}
 *@f$ so it should be a matrix of size nb_nodes_per_element @f$\times@f$
 *nb_degree_of_freedom
 * @param interpolated interpolated value of f @f$ f_j(\xi) = \sum_i f_{n_i j}
 *N_i @f$
 */
template <InterpolationType interpolation_type, InterpolationKind kind>
inline void
InterpolationElement<interpolation_type, kind>::interpolateOnNaturalCoordinates(
    const Vector<Real> & natural_coords, const Matrix<Real> & nodal_values,
    Vector<Real> & interpolated) {
  Vector<Real> shapes(
      InterpolationProperty<interpolation_type>::nb_nodes_per_element);
  computeShapes(natural_coords, shapes);

  interpolate(nodal_values, shapes, interpolated);
}

/* -------------------------------------------------------------------------- */
/// @f$ gradient_{ij} = \frac{\partial f_j}{\partial s_i} = \sum_k
/// \frac{\partial N_k}{\partial s_i}f_{j n_k} @f$
template <InterpolationType interpolation_type, InterpolationKind kind>
inline void
InterpolationElement<interpolation_type, kind>::gradientOnNaturalCoordinates(
    const Vector<Real> & natural_coords, const Matrix<Real> & f,
    Matrix<Real> & gradient) {
  Matrix<Real> dnds(
      InterpolationProperty<interpolation_type>::natural_space_dimension,
      InterpolationProperty<interpolation_type>::nb_nodes_per_element);
  computeDNDS(natural_coords, dnds);
  gradient.mul<false, true>(f, dnds);
}

/* -------------------------------------------------------------------------- */
/* ElementClass                                                               */
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
template <ElementType type, ElementKind kind>
inline void
ElementClass<type, kind>::computeJMat(const Tensor3<Real> & dnds,
                                      const Matrix<Real> & node_coords,
                                      Tensor3<Real> & J) {
  UInt nb_points = dnds.size(2);
  for (UInt p = 0; p < nb_points; ++p) {
    Matrix<Real> J_p(J(p));
    Matrix<Real> dnds_p(dnds(p));
    computeJMat(dnds_p, node_coords, J_p);
  }
}

/* -------------------------------------------------------------------------- */
template <ElementType type, ElementKind kind>
inline void
ElementClass<type, kind>::computeJMat(const Matrix<Real> & dnds,
                                      const Matrix<Real> & node_coords,
                                      Matrix<Real> & J) {
  /// @f$ J = dxds = dnds * x @f$
  J.mul<false, true>(dnds, node_coords);
}

/* -------------------------------------------------------------------------- */
template <ElementType type, ElementKind kind>
inline void
ElementClass<type, kind>::computeJacobian(const Matrix<Real> & natural_coords,
                                          const Matrix<Real> & node_coords,
                                          Vector<Real> & jacobians) {
  UInt nb_points = natural_coords.cols();
  Matrix<Real> dnds(interpolation_property::natural_space_dimension,
                    interpolation_property::nb_nodes_per_element);
  Matrix<Real> J(natural_coords.rows(), node_coords.rows());

  for (UInt p = 0; p < nb_points; ++p) {
    Vector<Real> ncoord_p(natural_coords(p));
    interpolation_element::computeDNDS(ncoord_p, dnds);
    computeJMat(dnds, node_coords, J);
    computeJacobian(J, jacobians(p));
  }
}

/* -------------------------------------------------------------------------- */
template <ElementType type, ElementKind kind>
inline void
ElementClass<type, kind>::computeJacobian(const Tensor3<Real> & J,
                                          Vector<Real> & jacobians) {
  UInt nb_points = J.size(2);
  for (UInt p = 0; p < nb_points; ++p) {
    computeJacobian(J(p), jacobians(p));
  }
}

/* -------------------------------------------------------------------------- */
template <ElementType type, ElementKind kind>
inline void ElementClass<type, kind>::computeJacobian(const Matrix<Real> & J,
                                                      Real & jacobians) {
  if (J.rows() == J.cols()) {
    jacobians = Math::det<element_property::spatial_dimension>(J.storage());
  } else {
    interpolation_element::computeSpecialJacobian(J, jacobians);
  }
}

/* -------------------------------------------------------------------------- */
template <ElementType type, ElementKind kind>
inline void
ElementClass<type, kind>::computeShapeDerivatives(const Tensor3<Real> & J,
                                                  const Tensor3<Real> & dnds,
                                                  Tensor3<Real> & shape_deriv) {
  UInt nb_points = J.size(2);
  for (UInt p = 0; p < nb_points; ++p) {
    Matrix<Real> shape_deriv_p(shape_deriv(p));
    computeShapeDerivatives(J(p), dnds(p), shape_deriv_p);
  }
}

/* -------------------------------------------------------------------------- */
template <ElementType type, ElementKind kind>
inline void
ElementClass<type, kind>::computeShapeDerivatives(const Matrix<Real> & J,
                                                  const Matrix<Real> & dnds,
                                                  Matrix<Real> & shape_deriv) {
  Matrix<Real> inv_J(J.rows(), J.cols());
  Math::inv<element_property::spatial_dimension>(J.storage(), inv_J.storage());

  shape_deriv.mul<false, false>(inv_J, dnds);
}

/* -------------------------------------------------------------------------- */
template <ElementType type, ElementKind kind>
inline void ElementClass<type, kind>::computeNormalsOnNaturalCoordinates(
    const Matrix<Real> & coord, Matrix<Real> & f, Matrix<Real> & normals) {
  UInt dimension = normals.rows();
  UInt nb_points = coord.cols();

  AKANTU_DEBUG_ASSERT((dimension - 1) ==
                          interpolation_property::natural_space_dimension,
                      "cannot extract a normal because of dimension mismatch "
                          << dimension - 1 << " "
                          << interpolation_property::natural_space_dimension);

  Matrix<Real> J(dimension, interpolation_property::natural_space_dimension);
  for (UInt p = 0; p < nb_points; ++p) {
    interpolation_element::gradientOnNaturalCoordinates(coord(p), f, J);
    if (dimension == 2) {
      Math::normal2(J.storage(), normals(p).storage());
    }
    if (dimension == 3) {
      Math::normal3(J(0).storage(), J(1).storage(), normals(p).storage());
    }
  }
}

/* ------------------------------------------------------------------------- */
/**
 * In the non linear cases we need to iterate to find the natural coordinates
 *@f$\xi@f$
 * provided real coordinates @f$x@f$.
 *
 * We want to solve: @f$ x- \phi(\xi) = 0@f$ with @f$\phi(\xi) = \sum_I N_I(\xi)
 *x_I@f$
 * the mapping function which uses the nodal coordinates @f$x_I@f$.
 *
 * To that end we use the Newton method and the following series:
 *
 * @f$ \frac{\partial \phi(x_k)}{\partial \xi} \left( \xi_{k+1} - \xi_k \right)
 *= x - \phi(x_k)@f$
 *
 * When we consider elements embedded in a dimension higher than them (2D
 *triangle in a 3D space for example)
 * @f$ J = \frac{\partial \phi(\xi_k)}{\partial \xi}@f$ is of dimension
 *@f$dim_{space} \times dim_{elem}@f$ which
 * is not invertible in most cases. Rather we can solve the problem:
 *
 * @f$ J^T J \left( \xi_{k+1} - \xi_k \right) = J^T \left( x - \phi(\xi_k)
 *\right) @f$
 *
 * So that
 *
 * @f$ d\xi = \xi_{k+1} - \xi_k = (J^T J)^{-1} J^T \left( x - \phi(\xi_k)
 *\right) @f$
 *
 * So that if the series converges we have:
 *
 * @f$ 0 = J^T \left( \phi(\xi_\infty) - x \right) @f$
 *
 * And we see that this is ill-posed only if @f$ J^T x = 0@f$ which means that
 *the vector provided
 * is normal to any tangent which means it is outside of the element itself.
 *
 * @param real_coords: the real coordinates the natural coordinates are sought
 *for
 * @param node_coords: the coordinates of the nodes forming the element
 * @param natural_coords: output->the sought natural coordinates
 * @param spatial_dimension: spatial dimension of the problem
 *
 **/
template <ElementType type, ElementKind kind>
inline void ElementClass<type, kind>::inverseMap(
    const Vector<Real> & real_coords, const Matrix<Real> & node_coords,
    Vector<Real> & natural_coords, Real tolerance) {
  UInt spatial_dimension = real_coords.size();
  UInt dimension = natural_coords.size();

  // matrix copy of the real_coords
  Matrix<Real> mreal_coords(real_coords.storage(), spatial_dimension, 1);

  // initial guess
  //  Matrix<Real> natural_guess(natural_coords.storage(), dimension, 1);
  natural_coords.zero();

  // real space coordinates provided by initial guess
  Matrix<Real> physical_guess(dimension, 1);

  // objective function f = real_coords - physical_guess
  Matrix<Real> f(dimension, 1);

  // dnds computed on the natural_guess
  //  Matrix<Real> dnds(interpolation_element::nb_nodes_per_element,
  //  spatial_dimension);

  // J Jacobian matrix computed on the natural_guess
  Matrix<Real> J(spatial_dimension, dimension);

  // G = J^t * J
  Matrix<Real> G(spatial_dimension, spatial_dimension);

  // Ginv = G^{-1}
  Matrix<Real> Ginv(spatial_dimension, spatial_dimension);

  // J = Ginv * J^t
  Matrix<Real> F(spatial_dimension, dimension);

  // dxi = \xi_{k+1} - \xi in the iterative process
  Matrix<Real> dxi(spatial_dimension, 1);

  /* --------------------------- */
  /* init before iteration loop  */
  /* --------------------------- */
  // do interpolation
  auto update_f = [&f, &physical_guess, &natural_coords, &node_coords,
                   &mreal_coords, dimension]() {
    Vector<Real> physical_guess_v(physical_guess.storage(), dimension);
    interpolation_element::interpolateOnNaturalCoordinates(
        natural_coords, node_coords, physical_guess_v);

    // compute initial objective function value f = real_coords - physical_guess
    f = mreal_coords;
    f -= physical_guess;

    // compute initial error
    auto error = f.norm<L_2>();
    return error;
  };

  auto inverse_map_error = update_f();

  /* --------------------------- */
  /* iteration loop              */
  /* --------------------------- */
  while (tolerance < inverse_map_error) {
    // compute J^t
    interpolation_element::gradientOnNaturalCoordinates(natural_coords,
                                                        node_coords, J);

    // compute G
    G.mul<true, false>(J, J);

    // inverse G
    Ginv.inverse(G);

    // compute F
    F.mul<false, true>(Ginv, J);

    // compute increment
    dxi.mul<false, false>(F, f);

    // update our guess
    natural_coords += Vector<Real>(dxi(0));

    inverse_map_error = update_f();
  }
  //  memcpy(natural_coords.storage(), natural_guess.storage(), sizeof(Real) *
  //  natural_coords.size());
}

/* -------------------------------------------------------------------------- */
template <ElementType type, ElementKind kind>
inline void ElementClass<type, kind>::inverseMap(
    const Matrix<Real> & real_coords, const Matrix<Real> & node_coords,
    Matrix<Real> & natural_coords, Real tolerance) {
  UInt nb_points = real_coords.cols();
  for (UInt p = 0; p < nb_points; ++p) {
    Vector<Real> X(real_coords(p));
    Vector<Real> ncoord_p(natural_coords(p));
    inverseMap(X, node_coords, ncoord_p, tolerance);
  }
}

} // namespace akantu

#endif /* AKANTU_ELEMENT_CLASS_TMPL_HH_ */

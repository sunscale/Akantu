/**
 * @file   integrator_gauss_inline_impl.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Tue Feb 15 2011
 * @date last modification: Tue Feb 20 2018
 *
 * @brief  inline function of gauss integrator
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
#include "fe_engine.hh"
#include "mesh_iterators.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {
namespace debug {
  struct IntegratorGaussException : public Exception {};
} // namespace debug
/* -------------------------------------------------------------------------- */
template <ElementKind kind, class IntegrationOrderFunctor>
template <ElementType type>
inline void IntegratorGauss<kind, IntegrationOrderFunctor>::integrateOnElement(
    const Array<Real> & f, Real * intf, UInt nb_degree_of_freedom,
    const UInt elem, GhostType ghost_type) const {
  Array<Real> & jac_loc = jacobians(type, ghost_type);

  UInt nb_quadrature_points = ElementClass<type>::getNbQuadraturePoints();
  AKANTU_DEBUG_ASSERT(f.getNbComponent() == nb_degree_of_freedom,
                      "The vector f do not have the good number of component.");

  Real * f_val = f.storage() + elem * f.getNbComponent();
  Real * jac_val = jac_loc.storage() + elem * nb_quadrature_points;

  integrate(f_val, jac_val, intf, nb_degree_of_freedom, nb_quadrature_points);
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind, class IntegrationOrderFunctor>
template <ElementType type>
inline Real IntegratorGauss<kind, IntegrationOrderFunctor>::integrate(
    const Vector<Real> & in_f, UInt index, GhostType ghost_type) const {
  const Array<Real> & jac_loc = jacobians(type, ghost_type);

  UInt nb_quadrature_points =
      GaussIntegrationElement<type>::getNbQuadraturePoints();
  AKANTU_DEBUG_ASSERT(in_f.size() == nb_quadrature_points,
                      "The vector f do not have nb_quadrature_points entries.");

  Real * jac_val = jac_loc.storage() + index * nb_quadrature_points;
  Real intf;

  integrate(in_f.storage(), jac_val, &intf, 1, nb_quadrature_points);

  return intf;
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind, class IntegrationOrderFunctor>
inline void IntegratorGauss<kind, IntegrationOrderFunctor>::integrate(
    Real * f, Real * jac, Real * inte, UInt nb_degree_of_freedom,
    UInt nb_quadrature_points) const {
  std::fill_n(inte, nb_degree_of_freedom, 0.);

  Real * cjac = jac;
  for (UInt q = 0; q < nb_quadrature_points; ++q) {
    for (UInt dof = 0; dof < nb_degree_of_freedom; ++dof) {
      inte[dof] += *f * *cjac;
      ++f;
    }
    ++cjac;
  }
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind, class IntegrationOrderFunctor>
template <ElementType type>
inline const Matrix<Real> &
IntegratorGauss<kind, IntegrationOrderFunctor>::getIntegrationPoints(
    GhostType ghost_type) const {
  AKANTU_DEBUG_ASSERT(
      quadrature_points.exists(type, ghost_type),
      "Quadrature points for type "
          << quadrature_points.printType(type, ghost_type)
          << " have not been initialized."
          << " Did you use 'computeQuadraturePoints' function ?");
  return quadrature_points(type, ghost_type);
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind, class IntegrationOrderFunctor>
template <ElementType type>
inline UInt
IntegratorGauss<kind, IntegrationOrderFunctor>::getNbIntegrationPoints(
    GhostType ghost_type) const {
  AKANTU_DEBUG_ASSERT(
      quadrature_points.exists(type, ghost_type),
      "Quadrature points for type "
          << quadrature_points.printType(type, ghost_type)
          << " have not been initialized."
          << " Did you use 'computeQuadraturePoints' function ?");
  return quadrature_points(type, ghost_type).cols();
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind, class IntegrationOrderFunctor>
template <ElementType type, UInt polynomial_degree>
inline Matrix<Real>
IntegratorGauss<kind, IntegrationOrderFunctor>::getIntegrationPoints() const {
  return GaussIntegrationElement<type,
                                 polynomial_degree>::getQuadraturePoints();
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind, class IntegrationOrderFunctor>
template <ElementType type, UInt polynomial_degree>
inline Vector<Real>
IntegratorGauss<kind, IntegrationOrderFunctor>::getIntegrationWeights() const {
  return GaussIntegrationElement<type, polynomial_degree>::getWeights();
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind, class IntegrationOrderFunctor>
template <ElementType type>
inline void
IntegratorGauss<kind, IntegrationOrderFunctor>::computeQuadraturePoints(
    GhostType ghost_type) {
  Matrix<Real> & quads = quadrature_points(type, ghost_type);
  const UInt polynomial_degree =
      IntegrationOrderFunctor::template getOrder<type>();
  quads =
      GaussIntegrationElement<type, polynomial_degree>::getQuadraturePoints();
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind, class IntegrationOrderFunctor>
template <ElementType type>
inline void IntegratorGauss<kind, IntegrationOrderFunctor>::
    computeJacobianOnQuadPointsByElement(const Matrix<Real> & node_coords,
                                         const Matrix<Real> & quad,
                                         Vector<Real> & jacobians) const {
  // jacobian
  ElementClass<type>::computeJacobian(quad, node_coords, jacobians);
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind, class IntegrationOrderFunctor>
IntegratorGauss<kind, IntegrationOrderFunctor>::IntegratorGauss(
    const Mesh & mesh, UInt spatial_dimension, const ID & id)
    : Integrator(mesh, spatial_dimension, id) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind, class IntegrationOrderFunctor>
template <ElementType type>
void IntegratorGauss<kind, IntegrationOrderFunctor>::checkJacobians(
    GhostType ghost_type) const {
  AKANTU_DEBUG_IN();

  UInt nb_quadrature_points = this->quadrature_points(type, ghost_type).cols();

  UInt nb_element = mesh.getConnectivity(type, ghost_type).size();

  Real * jacobians_val = jacobians(type, ghost_type).storage();

  for (UInt i = 0; i < nb_element * nb_quadrature_points;
       ++i, ++jacobians_val) {
    if (*jacobians_val < 0) {
      AKANTU_CUSTOM_EXCEPTION_INFO(debug::IntegratorGaussException{},
                                   "Negative jacobian computed,"
                                       << " possible problem in the element "
                                          "node ordering (Quadrature Point "
                                       << i % nb_quadrature_points << ":"
                                       << i / nb_quadrature_points << ":"
                                       << type << ":" << ghost_type << ")");
    }
  }
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind, class IntegrationOrderFunctor>
template <ElementType type>
void IntegratorGauss<kind, IntegrationOrderFunctor>::
    computeJacobiansOnIntegrationPoints(
        const Array<Real> & nodes, const Matrix<Real> & quad_points,
        Array<Real> & jacobians, GhostType ghost_type,
        const Array<UInt> & filter_elements) const {
  AKANTU_DEBUG_IN();

  UInt spatial_dimension = mesh.getSpatialDimension();
  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
  UInt nb_quadrature_points = quad_points.cols();

  UInt nb_element = mesh.getNbElement(type, ghost_type);
  jacobians.resize(nb_element * nb_quadrature_points);

  auto jacobians_it =
      jacobians.begin_reinterpret(nb_quadrature_points, nb_element);
  auto jacobians_begin = jacobians_it;

  Array<Real> x_el(0, spatial_dimension * nb_nodes_per_element);
  FEEngine::extractNodalToElementField(mesh, nodes, x_el, type, ghost_type,
                                       filter_elements);

  auto x_it = x_el.begin(spatial_dimension, nb_nodes_per_element);

  nb_element = x_el.size();

  //  Matrix<Real> local_coord(spatial_dimension, nb_nodes_per_element);
  for (UInt elem = 0; elem < nb_element; ++elem, ++x_it) {
    const Matrix<Real> & x = *x_it;
    if (filter_elements != empty_filter) {
      jacobians_it = jacobians_begin + filter_elements(elem);
    }

    Vector<Real> & J = *jacobians_it;
    computeJacobianOnQuadPointsByElement<type>(x, quad_points, J);

    if (filter_elements == empty_filter) {
      ++jacobians_it;
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
#if defined(AKANTU_STRUCTURAL_MECHANICS)
template <>
template <ElementType type>
void IntegratorGauss<_ek_structural, DefaultIntegrationOrderFunctor>::
    computeJacobiansOnIntegrationPoints(
        const Array<Real> & nodes, const Matrix<Real> & quad_points,
        Array<Real> & jacobians, GhostType ghost_type,
        const Array<UInt> & filter_elements) const {
  AKANTU_DEBUG_IN();

  const UInt spatial_dimension = mesh.getSpatialDimension();
  const UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
  const UInt nb_quadrature_points = quad_points.cols();
  const UInt nb_dofs = ElementClass<type>::getNbDegreeOfFreedom();

  UInt nb_element = mesh.getNbElement(type, ghost_type);
  jacobians.resize(nb_element * nb_quadrature_points);

  auto jacobians_it =
      jacobians.begin_reinterpret(nb_quadrature_points, nb_element);
  auto jacobians_begin = jacobians_it;

  Array<Real> x_el(0, spatial_dimension * nb_nodes_per_element);
  FEEngine::extractNodalToElementField(mesh, nodes, x_el, type, ghost_type,
                                       filter_elements);

  auto x_it = x_el.begin(spatial_dimension, nb_nodes_per_element);

  nb_element = x_el.size();

  const bool has_extra_normal =
      mesh.hasData<Real>("extra_normal", type, ghost_type);
  Array<Real>::const_vector_iterator extra_normal;
  Array<Real>::const_vector_iterator extra_normal_begin;
  if (has_extra_normal) {
    extra_normal = mesh.getData<Real>("extra_normal", type, ghost_type)
                       .begin(spatial_dimension);
    extra_normal_begin = extra_normal;
  }

  //  Matrix<Real> local_coord(spatial_dimension, nb_nodes_per_element);
  for (UInt elem = 0; elem < nb_element; ++elem, ++x_it) {
    if (filter_elements != empty_filter) {
      jacobians_it = jacobians_begin + filter_elements(elem);
      extra_normal = extra_normal_begin + filter_elements(elem);
    }

    const Matrix<Real> & X = *x_it;
    Vector<Real> & J = *jacobians_it;
    Matrix<Real> R(nb_dofs, nb_dofs);

    if (has_extra_normal) {
      ElementClass<type>::computeRotationMatrix(R, X, *extra_normal);
    } else {
      ElementClass<type>::computeRotationMatrix(R, X, Vector<Real>(X.rows()));
    }
    // Extracting relevant lines
    auto x = (R.block(0, 0, spatial_dimension, spatial_dimension) * X)
                 .block(0, 0, ElementClass<type>::getNaturalSpaceDimension(),
                        ElementClass<type>::getNbNodesPerElement());

    computeJacobianOnQuadPointsByElement<type>(x, quad_points, J);

    if (filter_elements == empty_filter) {
      ++jacobians_it;
      ++extra_normal;
    }
  }

  AKANTU_DEBUG_OUT();
}
#endif

/* -------------------------------------------------------------------------- */
#if defined(AKANTU_COHESIVE_ELEMENT)
template <>
template <ElementType type>
void IntegratorGauss<_ek_cohesive, DefaultIntegrationOrderFunctor>::
    computeJacobiansOnIntegrationPoints(
        const Array<Real> & nodes, const Matrix<Real> & quad_points,
        Array<Real> & jacobians, GhostType ghost_type,
        const Array<UInt> & filter_elements) const {
  AKANTU_DEBUG_IN();

  UInt spatial_dimension = mesh.getSpatialDimension();
  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
  UInt nb_quadrature_points = quad_points.cols();

  UInt nb_element = mesh.getNbElement(type, ghost_type);
  jacobians.resize(nb_element * nb_quadrature_points);

  auto jacobians_begin =
      jacobians.begin_reinterpret(nb_quadrature_points, nb_element);

  Array<Real> x_el(0, spatial_dimension * nb_nodes_per_element);
  FEEngine::extractNodalToElementField(mesh, nodes, x_el, type, ghost_type,
                                       filter_elements);

  auto x_it = x_el.begin(spatial_dimension, nb_nodes_per_element);

  UInt nb_nodes_per_subelement = nb_nodes_per_element / 2;
  Matrix<Real> x(spatial_dimension, nb_nodes_per_subelement);

  nb_element = x_el.size();
  UInt l_el = 0;
  auto compute = [&](const auto & el) {
    Vector<Real> J(jacobians_begin[el]);
    Matrix<Real> X(x_it[l_el]);
    ++l_el;

    for (UInt n = 0; n < nb_nodes_per_subelement; ++n) {
      Vector<Real>(x(n)) =
          (Vector<Real>(X(n)) + Vector<Real>(X(n + nb_nodes_per_subelement))) /
          2.;
    }

    if (type == _cohesive_1d_2) {
      J(0) = 1;
    } else {
      this->computeJacobianOnQuadPointsByElement<type>(x, quad_points, J);
    }
  };

  for_each_element(nb_element, filter_elements, compute);

  AKANTU_DEBUG_OUT();
}
#endif
/* -------------------------------------------------------------------------- */
template <ElementKind kind, class IntegrationOrderFunctor>
template <ElementType type>
void IntegratorGauss<kind, IntegrationOrderFunctor>::
    precomputeJacobiansOnQuadraturePoints(const Array<Real> & nodes,
                                          GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Array<Real> & jacobians_tmp = jacobians.alloc(0, 1, type, ghost_type);

  this->computeJacobiansOnIntegrationPoints<type>(
      nodes, quadrature_points(type, ghost_type), jacobians_tmp, ghost_type);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind, class IntegrationOrderFunctor>
template <ElementType type, UInt polynomial_degree>
void IntegratorGauss<kind, IntegrationOrderFunctor>::multiplyJacobiansByWeights(
    Array<Real> & jacobians, const Array<UInt> & filter_elements) const {
  AKANTU_DEBUG_IN();

  UInt nb_quadrature_points =
      GaussIntegrationElement<type, polynomial_degree>::getNbQuadraturePoints();

  Vector<Real> weights =
      GaussIntegrationElement<type, polynomial_degree>::getWeights();

  auto && view = make_view(jacobians, nb_quadrature_points);

  if (filter_elements != empty_filter) {
    auto J_it = view.begin();
    for (auto el : filter_elements) {
      Vector<Real> J(J_it[el]);
      J *= weights;
    }
  } else {
    for (auto & J : make_view(jacobians, nb_quadrature_points)) {
      J *= weights;
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind, class IntegrationOrderFunctor>
void IntegratorGauss<kind, IntegrationOrderFunctor>::integrate(
    const Array<Real> & in_f, Array<Real> & intf, UInt nb_degree_of_freedom,
    const Array<Real> & jacobians, UInt nb_element) const {
  AKANTU_DEBUG_IN();

  intf.resize(nb_element);
  if (nb_element == 0) {
    return;
  }

  UInt nb_points = jacobians.size() / nb_element;

  Array<Real>::const_matrix_iterator J_it;
  Array<Real>::matrix_iterator inte_it;
  Array<Real>::const_matrix_iterator f_it;

  f_it = in_f.begin_reinterpret(nb_degree_of_freedom, nb_points, nb_element);
  inte_it = intf.begin_reinterpret(nb_degree_of_freedom, 1, nb_element);
  J_it = jacobians.begin_reinterpret(nb_points, 1, nb_element);

  for (UInt el = 0; el < nb_element; ++el, ++J_it, ++f_it, ++inte_it) {
    const Matrix<Real> & f = *f_it;
    const Matrix<Real> & J = *J_it;
    Matrix<Real> & inte_f = *inte_it;

    inte_f.mul<false, false>(f, J);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind, class IntegrationOrderFunctor>
template <ElementType type>
void IntegratorGauss<kind, IntegrationOrderFunctor>::integrate(
    const Array<Real> & in_f, Array<Real> & intf, UInt nb_degree_of_freedom,
    GhostType ghost_type, const Array<UInt> & filter_elements) const {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_ASSERT(jacobians.exists(type, ghost_type),
                      "No jacobians for the type "
                          << jacobians.printType(type, ghost_type));

  const Array<Real> & jac_loc = jacobians(type, ghost_type);
  if (filter_elements != empty_filter) {
    UInt nb_element = filter_elements.size();
    auto * filtered_J = new Array<Real>(0, jac_loc.getNbComponent());
    FEEngine::filterElementalData(mesh, jac_loc, *filtered_J, type, ghost_type,
                                  filter_elements);
    this->integrate(in_f, intf, nb_degree_of_freedom, *filtered_J, nb_element);
    delete filtered_J;
  } else {
    UInt nb_element = mesh.getNbElement(type, ghost_type);
    this->integrate(in_f, intf, nb_degree_of_freedom, jac_loc, nb_element);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind, class IntegrationOrderFunctor>
template <ElementType type, UInt polynomial_degree>
void IntegratorGauss<kind, IntegrationOrderFunctor>::integrate(
    const Array<Real> & in_f, Array<Real> & intf, UInt nb_degree_of_freedom,
    GhostType ghost_type) const {
  AKANTU_DEBUG_IN();

  Matrix<Real> quads = this->getIntegrationPoints<type, polynomial_degree>();

  Array<Real> jacobians;
  this->computeJacobiansOnIntegrationPoints<type>(mesh.getNodes(), quads,
                                                  jacobians, ghost_type);
  this->multiplyJacobiansByWeights<type, polynomial_degree>(jacobians);

  this->integrate(in_f, intf, nb_degree_of_freedom, jacobians,
                  mesh.getNbElement(type, ghost_type));

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind, class IntegrationOrderFunctor>
template <ElementType type, UInt polynomial_degree>
Real IntegratorGauss<kind, IntegrationOrderFunctor>::integrate(
    const Array<Real> & in_f, GhostType ghost_type) const {
  AKANTU_DEBUG_IN();

  Array<Real> intfv(0, 1);
  integrate<type, polynomial_degree>(in_f, intfv, 1, ghost_type);

  Real res = Math::reduce(intfv);

  AKANTU_DEBUG_OUT();
  return res;
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind, class IntegrationOrderFunctor>
template <ElementType type>
Real IntegratorGauss<kind, IntegrationOrderFunctor>::integrate(
    const Array<Real> & in_f, GhostType ghost_type,
    const Array<UInt> & filter_elements) const {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_ASSERT(jacobians.exists(type, ghost_type),
                      "No jacobians for the type "
                          << jacobians.printType(type, ghost_type));

  Array<Real> intfv(0, 1);
  integrate<type>(in_f, intfv, 1, ghost_type, filter_elements);

  Real res = Math::reduce(intfv);

  AKANTU_DEBUG_OUT();
  return res;
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind, class IntegrationOrderFunctor>
void IntegratorGauss<kind, IntegrationOrderFunctor>::
    integrateOnIntegrationPoints(const Array<Real> & in_f, Array<Real> & intf,
                                 UInt nb_degree_of_freedom,
                                 const Array<Real> & jacobians,
                                 UInt nb_element) const {
  AKANTU_DEBUG_IN();

  UInt nb_points = jacobians.size() / nb_element;

  intf.resize(nb_element * nb_points);

  auto J_it = jacobians.begin();
  auto f_it = in_f.begin(nb_degree_of_freedom);
  auto inte_it = intf.begin(nb_degree_of_freedom);

  for (UInt el = 0; el < nb_element; ++el, ++J_it, ++f_it, ++inte_it) {
    const Real & J = *J_it;
    const Vector<Real> & f = *f_it;
    Vector<Real> & inte_f = *inte_it;

    inte_f = f;
    inte_f *= J;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind, class IntegrationOrderFunctor>
template <ElementType type>
void IntegratorGauss<kind, IntegrationOrderFunctor>::
    integrateOnIntegrationPoints(const Array<Real> & in_f, Array<Real> & intf,
                                 UInt nb_degree_of_freedom,
                                 GhostType ghost_type,
                                 const Array<UInt> & filter_elements) const {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_ASSERT(jacobians.exists(type, ghost_type),
                      "No jacobians for the type "
                          << jacobians.printType(type, ghost_type));

  const Array<Real> & jac_loc = this->jacobians(type, ghost_type);

  if (filter_elements != empty_filter) {

    UInt nb_element = filter_elements.size();
    auto * filtered_J = new Array<Real>(0, jac_loc.getNbComponent());
    FEEngine::filterElementalData(mesh, jac_loc, *filtered_J, type, ghost_type,
                                  filter_elements);

    this->integrateOnIntegrationPoints(in_f, intf, nb_degree_of_freedom,
                                       *filtered_J, nb_element);
  } else {
    UInt nb_element = mesh.getNbElement(type, ghost_type);
    this->integrateOnIntegrationPoints(in_f, intf, nb_degree_of_freedom,
                                       jac_loc, nb_element);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind, class IntegrationOrderFunctor>
template <ElementType type>
inline void
IntegratorGauss<kind, IntegrationOrderFunctor>::onElementsAddedByType(
    const Array<UInt> & elements, GhostType ghost_type) {
  const auto & nodes = mesh.getNodes();

  if (not quadrature_points.exists(type, ghost_type)) {
    computeQuadraturePoints<type>(ghost_type);
  }

  if (not jacobians.exists(type, ghost_type)) {
    jacobians.alloc(0, 1, type, ghost_type);
  }

  this->computeJacobiansOnIntegrationPoints(
      nodes, quadrature_points(type, ghost_type), jacobians(type, ghost_type),
      type, ghost_type, elements);

  constexpr UInt polynomial_degree =
      IntegrationOrderFunctor::template getOrder<type>();

  multiplyJacobiansByWeights<type, polynomial_degree>(
      this->jacobians(type, ghost_type), elements);
}

/* -------------------------------------------------------------------------- */
namespace integrator {
  namespace details {
    template <ElementKind kind> struct IntegratorOnElementsAddedHelper {};

#define ON_ELEMENT_ADDED(type)                                                 \
  integrator.template onElementsAddedByType<type>(elements, ghost_type);

#define AKANTU_SPECIALIZE_ON_ELEMENT_ADDED_HELPER(kind)                        \
  template <> struct IntegratorOnElementsAddedHelper<kind> {                   \
    template <class I>                                                         \
    static void call(I & integrator, const Array<UInt> & elements,             \
                     ElementType type, GhostType ghost_type) {                 \
      AKANTU_BOOST_KIND_ELEMENT_SWITCH(ON_ELEMENT_ADDED, kind);                \
    }                                                                          \
  };

    AKANTU_BOOST_ALL_KIND(AKANTU_SPECIALIZE_ON_ELEMENT_ADDED_HELPER)

#undef AKANTU_SPECIALIZE_ON_ELEMENT_ADDED_HELPER
#undef ON_ELEMENT_ADDED
  } // namespace details
} // namespace integrator

/* -------------------------------------------------------------------------- */
template <ElementKind kind, class IntegrationOrderFunctor>
void IntegratorGauss<kind, IntegrationOrderFunctor>::onElementsAdded(
    const Array<Element> & new_elements) {
  for (auto elements_range : MeshElementsByTypes(new_elements)) {
    auto type = elements_range.getType();
    auto ghost_type = elements_range.getGhostType();

    if (mesh.getSpatialDimension(type) != _spatial_dimension) {
      continue;
    }

    if (mesh.getKind(type) != kind) {
      continue;
    }

    integrator::details::IntegratorOnElementsAddedHelper<kind>::call(
        *this, elements_range.getElements(), type, ghost_type);
  }
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind, class IntegrationOrderFunctor>
template <ElementType type>
inline void IntegratorGauss<kind, IntegrationOrderFunctor>::initIntegrator(
    const Array<Real> & nodes, GhostType ghost_type) {
  computeQuadraturePoints<type>(ghost_type);
  precomputeJacobiansOnQuadraturePoints<type>(nodes, ghost_type);
  checkJacobians<type>(ghost_type);
  constexpr UInt polynomial_degree =
      IntegrationOrderFunctor::template getOrder<type>();
  multiplyJacobiansByWeights<type, polynomial_degree>(
      this->jacobians(type, ghost_type));
}

namespace integrator {
  namespace details {
    template <ElementKind kind> struct GaussIntegratorInitHelper {};

#define INIT_INTEGRATOR(type)                                                  \
  _int.template initIntegrator<type>(nodes, ghost_type)

#define AKANTU_GAUSS_INTERGRATOR_INIT_HELPER(kind)                             \
  template <> struct GaussIntegratorInitHelper<kind> {                         \
    template <ElementKind k, class IOF>                                        \
    static void call(IntegratorGauss<k, IOF> & _int,                           \
                     const Array<Real> & nodes, ElementType type,              \
                     GhostType ghost_type) {                                   \
      AKANTU_BOOST_KIND_ELEMENT_SWITCH(INIT_INTEGRATOR, kind);                 \
    }                                                                          \
  };

    AKANTU_BOOST_ALL_KIND(AKANTU_GAUSS_INTERGRATOR_INIT_HELPER)

#undef AKANTU_GAUSS_INTERGRATOR_INIT_HELPER
#undef INIT_INTEGRATOR
  } // namespace details
} // namespace integrator

template <ElementKind kind, class IntegrationOrderFunctor>
inline void IntegratorGauss<kind, IntegrationOrderFunctor>::initIntegrator(
    const Array<Real> & nodes, ElementType type, GhostType ghost_type) {
  integrator::details::GaussIntegratorInitHelper<kind>::call(*this, nodes, type,
                                                             ghost_type);
}

namespace integrator {
  namespace details {
    template <ElementKind kind> struct GaussIntegratorComputeJacobiansHelper {};

#define AKANTU_COMPUTE_JACOBIANS(type)                                         \
  _int.template computeJacobiansOnIntegrationPoints<type>(                     \
      nodes, quad_points, jacobians, ghost_type, filter_elements);

#define AKANTU_GAUSS_INTERGRATOR_COMPUTE_JACOBIANS(kind)                       \
  template <> struct GaussIntegratorComputeJacobiansHelper<kind> {             \
    template <ElementKind k, class IOF>                                        \
    static void                                                                \
    call(const IntegratorGauss<k, IOF> & _int, const Array<Real> & nodes,      \
         const Matrix<Real> & quad_points, Array<Real> & jacobians,            \
         ElementType type, GhostType ghost_type,                               \
         const Array<UInt> & filter_elements) {                                \
      AKANTU_BOOST_KIND_ELEMENT_SWITCH(AKANTU_COMPUTE_JACOBIANS, kind);        \
    }                                                                          \
  };

    AKANTU_BOOST_ALL_KIND(AKANTU_GAUSS_INTERGRATOR_COMPUTE_JACOBIANS)

#undef AKANTU_GAUSS_INTERGRATOR_COMPUTE_JACOBIANS
#undef AKANTU_COMPUTE_JACOBIANS
  } // namespace details
} // namespace integrator

template <ElementKind kind, class IntegrationOrderFunctor>
void IntegratorGauss<kind, IntegrationOrderFunctor>::
    computeJacobiansOnIntegrationPoints(
        const Array<Real> & nodes, const Matrix<Real> & quad_points,
        Array<Real> & jacobians, ElementType type, GhostType ghost_type,
        const Array<UInt> & filter_elements) const {
  integrator::details::GaussIntegratorComputeJacobiansHelper<kind>::call(
      *this, nodes, quad_points, jacobians, type, ghost_type, filter_elements);
}

} // namespace akantu

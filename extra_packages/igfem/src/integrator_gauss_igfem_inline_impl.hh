/**
 * @file   integrator_gauss_igfem.hh
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 *
 *
 * @brief  Inline functions of gauss integrator for the case of IGFEM
 *
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

/* -------------------------------------------------------------------------- */

} // namespace akantu

#include "fe_engine.hh"
#if defined(AKANTU_DEBUG_TOOLS)
#include "aka_debug_tools.hh"
#endif

namespace akantu {

/* -------------------------------------------------------------------------- */
#define INIT_INTEGRATOR(type)                                                  \
  computeQuadraturePoints<type>(ghost_type);                                   \
  precomputeJacobiansOnQuadraturePoints<type>(nodes, ghost_type);              \
  checkJacobians<type>(ghost_type);

template <class IOF>
inline void
IntegratorGauss<_ek_igfem, IOF>::initIntegrator(const Array<Real> & nodes,
                                                const ElementType & type,
                                                const GhostType & ghost_type) {
  AKANTU_BOOST_IGFEM_ELEMENT_SWITCH(INIT_INTEGRATOR);
}

#undef INIT_INTEGRATOR

/* -------------------------------------------------------------------------- */
template <class IOF>
template <ElementType type>
inline UInt IntegratorGauss<_ek_igfem, IOF>::getNbIntegrationPoints(
    const GhostType &) const {
  const ElementType sub_type_1 = ElementClassProperty<type>::sub_element_type_1;
  const ElementType sub_type_2 = ElementClassProperty<type>::sub_element_type_2;
  UInt nb_quad_points_sub_1 =
      GaussIntegrationElement<sub_type_1>::getNbQuadraturePoints();
  UInt nb_quad_points_sub_2 =
      GaussIntegrationElement<sub_type_2>::getNbQuadraturePoints();
  return (nb_quad_points_sub_1 + nb_quad_points_sub_2);
}

/* -------------------------------------------------------------------------- */
template <class IOF>
template <ElementType type>
inline void IntegratorGauss<_ek_igfem, IOF>::integrateOnElement(
    const Array<Real> & f, Real * intf, UInt nb_degree_of_freedom,
    const UInt elem, const GhostType & ghost_type) const {

  Array<Real> & jac_loc = jacobians(type, ghost_type);

  UInt nb_quadrature_points = getNbIntegrationPoints<type>();
  AKANTU_DEBUG_ASSERT(f.getNbComponent() == nb_degree_of_freedom,
                      "The vector f do not have the good number of component.");

  Real * f_val = f.storage() + elem * f.getNbComponent();
  Real * jac_val = jac_loc.storage() + elem * nb_quadrature_points;

  integrate(f_val, jac_val, intf, nb_degree_of_freedom, nb_quadrature_points);
}

/* -------------------------------------------------------------------------- */
template <class IOF>
template <ElementType type>
inline Real IntegratorGauss<_ek_igfem, IOF>::integrate(
    const Vector<Real> & in_f, UInt index, const GhostType & ghost_type) const {
  const Array<Real> & jac_loc = jacobians(type, ghost_type);

  UInt nb_quadrature_points = getNbIntegrationPoints<type>();
  AKANTU_DEBUG_ASSERT(in_f.size() == nb_quadrature_points,
                      "The vector f do not have nb_quadrature_points entries.");

  Real * jac_val = jac_loc.storage() + index * nb_quadrature_points;
  Real intf;

  integrate(in_f.storage(), jac_val, &intf, 1, nb_quadrature_points);

  return intf;
  return 0.;
}

/* -------------------------------------------------------------------------- */
template <class IOF>
inline void
IntegratorGauss<_ek_igfem, IOF>::integrate(Real * f, Real * jac, Real * inte,
                                           UInt nb_degree_of_freedom,
                                           UInt nb_quadrature_points) const {
  memset(inte, 0, nb_degree_of_freedom * sizeof(Real));

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
template <class IOF>
template <ElementType type>
inline const Matrix<Real> &
IntegratorGauss<_ek_igfem, IOF>::getIntegrationPoints(
    const GhostType & ghost_type) const {
  AKANTU_DEBUG_ASSERT(
      quadrature_points.exists(type, ghost_type),
      "Quadrature points for type "
          << quadrature_points.printType(type, ghost_type)
          << " have not been initialized."
          << " Did you use 'computeQuadraturePoints' function ?");
  return quadrature_points(type, ghost_type);
}

/* -------------------------------------------------------------------------- */
template <class IOF>
template <ElementType type>
inline void IntegratorGauss<_ek_igfem, IOF>::computeQuadraturePoints(
    const GhostType & ghost_type) {
  /// typedef for the two subelement_types and the parent element type
  const ElementType sub_type_1 = ElementClassProperty<type>::sub_element_type_1;
  const ElementType sub_type_2 = ElementClassProperty<type>::sub_element_type_2;

  /// store the quadrature points on the two subelements
  Matrix<Real> & quads_sub_1 = quadrature_points(sub_type_1, ghost_type);
  Matrix<Real> & quads_sub_2 = quadrature_points(sub_type_2, ghost_type);
  quads_sub_1 = GaussIntegrationElement<sub_type_1>::getQuadraturePoints();
  quads_sub_2 = GaussIntegrationElement<sub_type_2>::getQuadraturePoints();

  /// store all quad points for the current type
  UInt nb_quad_points_sub_1 =
      GaussIntegrationElement<sub_type_1>::getNbQuadraturePoints();
  UInt nb_quad_points_sub_2 =
      GaussIntegrationElement<sub_type_2>::getNbQuadraturePoints();

  UInt spatial_dimension = mesh.getSpatialDimension();

  Matrix<Real> & quads = quadrature_points(type, ghost_type);
  quads = Matrix<Real>(spatial_dimension,
                       nb_quad_points_sub_1 + nb_quad_points_sub_2);

  Matrix<Real> quads_1(quads.storage(), quads.rows(), nb_quad_points_sub_1);
  quads_1 = quads_sub_1;

  Matrix<Real> quads_2(quads.storage() + quads.rows() * nb_quad_points_sub_1,
                       quads.rows(), nb_quad_points_sub_2);

  quads_2 = quads_sub_2;
}

/* -------------------------------------------------------------------------- */
template <class IOF>
template <ElementType type>
inline void
IntegratorGauss<_ek_igfem, IOF>::computeJacobianOnQuadPointsByElement(
    const Matrix<Real> & node_coords, Vector<Real> & jacobians) {

  /// optimize: get the matrix from the ElementTypeMap
  Matrix<Real> quad = GaussIntegrationElement<type>::getQuadraturePoints();
  // jacobian
  ElementClass<type>::computeJacobian(quad, node_coords, jacobians);
}

/* -------------------------------------------------------------------------- */
template <class IOF>
inline IntegratorGauss<_ek_igfem, IOF>::IntegratorGauss(
    const Mesh & mesh, const ID & id, const MemoryID & memory_id)
    : Integrator(mesh, id, memory_id) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class IOF>
template <ElementType type>
inline void IntegratorGauss<_ek_igfem, IOF>::checkJacobians(
    const GhostType & ghost_type) const {
  AKANTU_DEBUG_IN();
  /// typedef for the two subelement_types and the parent element type
  const ElementType sub_type_1 = ElementClassProperty<type>::sub_element_type_1;
  const ElementType sub_type_2 = ElementClassProperty<type>::sub_element_type_2;
  UInt nb_quad_points_sub_1 =
      GaussIntegrationElement<sub_type_1>::getNbQuadraturePoints();
  UInt nb_quad_points_sub_2 =
      GaussIntegrationElement<sub_type_2>::getNbQuadraturePoints();
  UInt nb_quadrature_points = nb_quad_points_sub_1 + nb_quad_points_sub_2;

  UInt nb_element;

  nb_element = mesh.getConnectivity(type, ghost_type).getSize();

  Real * jacobians_val = jacobians(type, ghost_type).storage();

  for (UInt i = 0; i < nb_element * nb_quadrature_points;
       ++i, ++jacobians_val) {
    if (*jacobians_val < 0)
      AKANTU_ERROR(
          "Negative jacobian computed,"
          << " possible problem in the element node ordering (Quadrature Point "
          << i % nb_quadrature_points << ":" << i / nb_quadrature_points << ":"
          << type << ":" << ghost_type << ")");
  }
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class IOF>
template <ElementType type>
inline void
IntegratorGauss<_ek_igfem, IOF>::precomputeJacobiansOnQuadraturePoints(
    const Array<Real> & nodes, const GhostType & ghost_type) {
  AKANTU_DEBUG_IN();

  /// typedef for the two subelement_types and the parent element type
  const ElementType sub_type_1 = ElementClassProperty<type>::sub_element_type_1;
  const ElementType sub_type_2 = ElementClassProperty<type>::sub_element_type_2;

  UInt spatial_dimension = mesh.getSpatialDimension();
  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
  /// get the number of nodes for the subelements and the parent element
  UInt nb_nodes_sub_1 =
      ElementClass<sub_type_1>::getNbNodesPerInterpolationElement();
  UInt nb_nodes_sub_2 =
      ElementClass<sub_type_2>::getNbNodesPerInterpolationElement();
  UInt nb_quadrature_points_sub_1 =
      GaussIntegrationElement<sub_type_1>::getNbQuadraturePoints();
  UInt nb_quadrature_points_sub_2 =
      GaussIntegrationElement<sub_type_2>::getNbQuadraturePoints();
  UInt nb_quadrature_points =
      nb_quadrature_points_sub_1 + nb_quadrature_points_sub_2;

  UInt nb_element = mesh.getNbElement(type, ghost_type);

  Array<Real> * jacobians_tmp;
  if (!jacobians.exists(type, ghost_type))
    jacobians_tmp = &jacobians.alloc(nb_element * nb_quadrature_points, 1, type,
                                     ghost_type);
  else {
    jacobians_tmp = &jacobians(type, ghost_type);
    jacobians_tmp->resize(nb_element * nb_quadrature_points);
  }

  Array<Real>::vector_iterator jacobians_it =
      jacobians_tmp->begin_reinterpret(nb_quadrature_points, nb_element);

  Vector<Real> weights_sub_1 =
      GaussIntegrationElement<sub_type_1>::getWeights();
  Vector<Real> weights_sub_2 =
      GaussIntegrationElement<sub_type_2>::getWeights();

  Array<Real> x_el(0, spatial_dimension * nb_nodes_per_element);
  FEEngine::extractNodalToElementField(mesh, nodes, x_el, type, ghost_type);

  Array<Real>::const_matrix_iterator x_it =
      x_el.begin(spatial_dimension, nb_nodes_per_element);

  //  Matrix<Real> local_coord(spatial_dimension, nb_nodes_per_element);
  for (UInt elem = 0; elem < nb_element; ++elem, ++jacobians_it, ++x_it) {
    const Matrix<Real> & X = *x_it;

    Matrix<Real> sub_1_coords(spatial_dimension, nb_nodes_sub_1);
    Matrix<Real> sub_2_coords(spatial_dimension, nb_nodes_sub_2);

    ElementClass<type>::getSubElementCoords(X, sub_1_coords, 0);
    ElementClass<type>::getSubElementCoords(X, sub_2_coords, 1);

    Vector<Real> & J = *jacobians_it;
    /// initialize vectors to store the jacobians for each subelement
    Vector<Real> J_sub_1(nb_quadrature_points_sub_1);
    Vector<Real> J_sub_2(nb_quadrature_points_sub_2);
    computeJacobianOnQuadPointsByElement<sub_type_1>(sub_1_coords, J_sub_1);
    computeJacobianOnQuadPointsByElement<sub_type_2>(sub_2_coords, J_sub_2);
    J_sub_1 *= weights_sub_1;
    J_sub_2 *= weights_sub_2;

    /// copy results into the jacobian vector for this element
    for (UInt i = 0; i < nb_quadrature_points_sub_1; ++i) {
      J(i) = J_sub_1(i);
    }

    for (UInt i = 0; i < nb_quadrature_points_sub_2; ++i) {
      J(i + nb_quadrature_points_sub_1) = J_sub_2(i);
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class IOF>
template <ElementType type>
inline void IntegratorGauss<_ek_igfem, IOF>::integrate(
    const Array<Real> & in_f, Array<Real> & intf, UInt nb_degree_of_freedom,
    const GhostType & ghost_type, const Array<UInt> & filter_elements) const {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_ASSERT(jacobians.exists(type, ghost_type),
                      "No jacobians for the type "
                          << jacobians.printType(type, ghost_type));

  const Matrix<Real> & quads = quadrature_points(type, ghost_type);

  UInt nb_points = quads.cols();

  const Array<Real> & jac_loc = jacobians(type, ghost_type);

  Array<Real>::const_matrix_iterator J_it;
  Array<Real>::matrix_iterator inte_it;
  Array<Real>::const_matrix_iterator f_it;

  UInt nb_element;
  Array<Real> * filtered_J = NULL;
  if (filter_elements != empty_filter) {
    nb_element = filter_elements.getSize();
    filtered_J = new Array<Real>(0, jac_loc.getNbComponent());
    FEEngine::filterElementalData(mesh, jac_loc, *filtered_J, type, ghost_type,
                                  filter_elements);
    const Array<Real> & cfiltered_J = *filtered_J; // \todo temporary patch
    J_it = cfiltered_J.begin_reinterpret(nb_points, 1, nb_element);
  } else {
    nb_element = mesh.getNbElement(type, ghost_type);
    J_it = jac_loc.begin_reinterpret(nb_points, 1, nb_element);
  }

  intf.resize(nb_element);

  f_it = in_f.begin_reinterpret(nb_degree_of_freedom, nb_points, nb_element);
  inte_it = intf.begin_reinterpret(nb_degree_of_freedom, 1, nb_element);

  for (UInt el = 0; el < nb_element; ++el, ++J_it, ++f_it, ++inte_it) {
    const Matrix<Real> & f = *f_it;
    const Matrix<Real> & J = *J_it;
    Matrix<Real> & inte_f = *inte_it;

    inte_f.mul<false, false>(f, J);
  }

  delete filtered_J;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class IOF>
template <ElementType type>
inline Real IntegratorGauss<_ek_igfem, IOF>::integrate(
    const Array<Real> & in_f, const GhostType & ghost_type,
    const Array<UInt> & filter_elements) const {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_ASSERT(jacobians.exists(type, ghost_type),
                      "No jacobians for the type "
                          << jacobians.printType(type, ghost_type));

  Array<Real> intfv(0, 1);
  integrate<type>(in_f, intfv, 1, ghost_type, filter_elements);

  UInt nb_values = intfv.getSize();
  if (nb_values == 0)
    return 0.;

  UInt nb_values_to_sum = nb_values >> 1;

  std::sort(intfv.begin(), intfv.end());

  // as long as the half is not empty
  while (nb_values_to_sum) {
    UInt remaining = (nb_values - 2 * nb_values_to_sum);
    if (remaining)
      intfv(nb_values - 2) += intfv(nb_values - 1);

    // sum to consecutive values and store the sum in the first half
    for (UInt i = 0; i < nb_values_to_sum; ++i) {
      intfv(i) = intfv(2 * i) + intfv(2 * i + 1);
    }

    nb_values = nb_values_to_sum;
    nb_values_to_sum >>= 1;
  }

  AKANTU_DEBUG_OUT();
  return intfv(0);
}

/* -------------------------------------------------------------------------- */
template <class IOF>
template <ElementType type>
inline void IntegratorGauss<_ek_igfem, IOF>::integrateOnIntegrationPoints(
    const Array<Real> & in_f, Array<Real> & intf, UInt nb_degree_of_freedom,
    const GhostType & ghost_type, const Array<UInt> & filter_elements) const {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_ASSERT(jacobians.exists(type, ghost_type),
                      "No jacobians for the type "
                          << jacobians.printType(type, ghost_type));

  UInt nb_element;
  const Matrix<Real> & quads = quadrature_points(type, ghost_type);

  UInt nb_points = quads.cols();

  const Array<Real> & jac_loc = jacobians(type, ghost_type);

  Array<Real>::const_scalar_iterator J_it;
  Array<Real>::vector_iterator inte_it;
  Array<Real>::const_vector_iterator f_it;

  Array<Real> * filtered_J = NULL;
  if (filter_elements != empty_filter) {
    nb_element = filter_elements.getSize();
    filtered_J = new Array<Real>(0, jac_loc.getNbComponent());
    FEEngine::filterElementalData(mesh, jac_loc, *filtered_J, type, ghost_type,
                                  filter_elements);
    J_it = filtered_J->begin();
  } else {
    nb_element = mesh.getNbElement(type, ghost_type);
    J_it = jac_loc.begin();
  }

  intf.resize(nb_element * nb_points);

  f_it = in_f.begin(nb_degree_of_freedom);
  inte_it = intf.begin(nb_degree_of_freedom);

  for (UInt el = 0; el < nb_element; ++el, ++J_it, ++f_it, ++inte_it) {
    const Real & J = *J_it;
    const Vector<Real> & f = *f_it;
    Vector<Real> & inte_f = *inte_it;

    inte_f = f;
    inte_f *= J;
  }

  delete filtered_J;

  AKANTU_DEBUG_OUT();
}

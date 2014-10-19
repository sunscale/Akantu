/**
 * @file   integrator_gauss_inline_impl.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Tue Feb 15 2011
 * @date last modification: Mon Jun 23 2014
 *
 * @brief  inline function of gauss integrator
 *
 * @section LICENSE
 *
 * Copyright (©) 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

__END_AKANTU__

#include "fe_engine.hh"
#if defined(AKANTU_DEBUG_TOOLS)
#  include "aka_debug_tools.hh"
#endif

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
inline void IntegratorGauss<kind>::initIntegrator(const Array<Real> & nodes,
						  const ElementType & type,
						  const GhostType & ghost_type) {
#define INIT_INTEGRATOR(type)						\
  computeQuadraturePoints<type>(ghost_type);				\
  precomputeJacobiansOnQuadraturePoints<type>(nodes, ghost_type);	\
  checkJacobians<type>(ghost_type);

  AKANTU_BOOST_ALL_ELEMENT_SWITCH(INIT_INTEGRATOR);
#undef INIT_INTEGRATOR
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <ElementType type>
inline void IntegratorGauss<kind>::integrateOnElement(const Array<Real> & f,
						      Real * intf,
						      UInt nb_degree_of_freedom,
						      const UInt elem,
						      const GhostType & ghost_type) const {
  Array<Real> & jac_loc = jacobians(type, ghost_type);

  UInt nb_quadrature_points = ElementClass<type>::getNbQuadraturePoints();
  AKANTU_DEBUG_ASSERT(f.getNbComponent() == nb_degree_of_freedom ,
		      "The vector f do not have the good number of component.");

  Real * f_val    = f.storage() + elem * f.getNbComponent();
  Real * jac_val  = jac_loc.storage() + elem * nb_quadrature_points;

  integrate(f_val, jac_val, intf, nb_degree_of_freedom, nb_quadrature_points);
}


/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <ElementType type>
inline Real IntegratorGauss<kind>::integrate(const Vector<Real> & in_f,
					     UInt index,
					     const GhostType & ghost_type) const {
  const Array<Real> & jac_loc = jacobians(type, ghost_type);

  UInt nb_quadrature_points = GaussIntegrationElement<type>::getNbQuadraturePoints();
  AKANTU_DEBUG_ASSERT(in_f.size() == nb_quadrature_points ,
		      "The vector f do not have nb_quadrature_points entries.");

  Real * jac_val  = jac_loc.storage() + index * nb_quadrature_points;
  Real intf;

  integrate(in_f.storage(), jac_val, &intf, 1, nb_quadrature_points);

  return intf;
}


/* -------------------------------------------------------------------------- */
template <ElementKind kind>
inline void IntegratorGauss<kind>::integrate(Real *f, Real *jac, Real * inte,
					     UInt nb_degree_of_freedom,
					     UInt nb_quadrature_points) const {
  memset(inte, 0, nb_degree_of_freedom * sizeof(Real));

  Real *cjac = jac;
  for (UInt q = 0; q < nb_quadrature_points; ++q) {
    for (UInt dof = 0; dof < nb_degree_of_freedom; ++dof) {
      inte[dof] += *f * *cjac;
      ++f;
    }
    ++cjac;
  }
}




/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <ElementType type>
inline const Matrix<Real> & IntegratorGauss<kind>::getQuadraturePoints(const GhostType & ghost_type) const {
  AKANTU_DEBUG_ASSERT(quadrature_points.exists(type, ghost_type),
		      "Quadrature points for type "
		      << quadrature_points.printType(type, ghost_type)
 		      << " have not been initialized."
		      << " Did you use 'computeQuadraturePoints' function ?");
  return quadrature_points(type, ghost_type);
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <ElementType type>
inline void IntegratorGauss<kind>::computeQuadraturePoints(const GhostType & ghost_type) {
  Matrix<Real> & quads = quadrature_points(type, ghost_type);
  quads = GaussIntegrationElement<type>::getQuadraturePoints();
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <ElementType type>
inline void IntegratorGauss<kind>::
computeJacobianOnQuadPointsByElement(const Matrix<Real> & node_coords,
				     Vector<Real> & jacobians) {

  Matrix<Real> quad = GaussIntegrationElement<type>::getQuadraturePoints();
  // jacobian
  ElementClass<type>::computeJacobian(quad, node_coords, jacobians);
}

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
template <ElementKind kind>
IntegratorGauss<kind>::IntegratorGauss(const Mesh & mesh,
				       const ID & id,
				       const MemoryID & memory_id) :
  Integrator(mesh, id, memory_id) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <ElementType type>
void IntegratorGauss<kind>::checkJacobians(const GhostType & ghost_type) const {
  AKANTU_DEBUG_IN();

  UInt nb_quadrature_points = GaussIntegrationElement<type>::getNbQuadraturePoints();

  UInt nb_element;

  nb_element = mesh.getConnectivity(type,ghost_type).getSize();

  Real * jacobians_val = jacobians(type, ghost_type).storage();

  for (UInt i = 0; i < nb_element*nb_quadrature_points; ++i,++jacobians_val){
    if(*jacobians_val < 0)
      AKANTU_DEBUG_ERROR("Negative jacobian computed,"
			 << " possible problem in the element node ordering (Quadrature Point "
			 << i % nb_quadrature_points << ":"
			 << i / nb_quadrature_points << ":"
			 << type << ":"
			 << ghost_type << ")");
  }
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <ElementType type>
void IntegratorGauss<kind>::precomputeJacobiansOnQuadraturePoints(const Array<Real> & nodes,
								  const GhostType & ghost_type) {
  AKANTU_DEBUG_IN();

  UInt spatial_dimension    = mesh.getSpatialDimension();
  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
  UInt nb_quadrature_points = GaussIntegrationElement<type>::getNbQuadraturePoints();

  UInt nb_element = mesh.getNbElement(type,ghost_type);

  Array<Real> * jacobians_tmp;
  if(!jacobians.exists(type, ghost_type))
    jacobians_tmp = &jacobians.alloc(nb_element*nb_quadrature_points,
                                     1,
                                     type,
                                     ghost_type);
  else {
    jacobians_tmp = &jacobians(type, ghost_type);
    jacobians_tmp->resize(nb_element*nb_quadrature_points);
  }

  Array<Real>::vector_iterator jacobians_it =
    jacobians_tmp->begin_reinterpret(nb_quadrature_points, nb_element);

  Vector<Real> weights = GaussIntegrationElement<type>::getWeights();

  Array<Real> x_el(0, spatial_dimension * nb_nodes_per_element);
  FEEngine::extractNodalToElementField(mesh, nodes, x_el, type, ghost_type);

  Array<Real>::const_matrix_iterator x_it = x_el.begin(spatial_dimension,
						       nb_nodes_per_element);

  //  Matrix<Real> local_coord(spatial_dimension, nb_nodes_per_element);
  for (UInt elem = 0; elem < nb_element; ++elem, ++jacobians_it, ++x_it) {
    const Matrix<Real> & x = *x_it;
    Vector<Real> & J = *jacobians_it;
    computeJacobianOnQuadPointsByElement<type>(x, J);
    J *= weights;
  }

  // >>>>>> DEBUG CODE >>>>>> //
#if defined(AKANTU_DEBUG_TOOLS)
#if defined(AKANTU_CORE_CXX11)
  debug::element_manager.print(debug::_dm_integrator,
                               [ghost_type, this,
                                nb_element, nb_quadrature_points](const Element & el)->std::string {
                                 std::stringstream out;
                                 if(el.ghost_type == ghost_type) {
                                   Array<Real>::vector_iterator jacobians_it =
                                     jacobians(el.type, el.ghost_type).begin(nb_quadrature_points);
                                   out << " jacobian: " << jacobians_it[el.element];
                                 }
                                 return out.str();
                               });
#endif
#endif
  // <<<<<< DEBUG CODE <<<<<< //

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
#if defined(AKANTU_COHESIVE_ELEMENT)
template <>
template <ElementType type>
void IntegratorGauss<_ek_cohesive>::precomputeJacobiansOnQuadraturePoints(const Array<Real> & nodes,
									  const GhostType & ghost_type) {
  AKANTU_DEBUG_IN();

  UInt spatial_dimension    = mesh.getSpatialDimension();
  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
  UInt nb_quadrature_points = GaussIntegrationElement<type>::getNbQuadraturePoints();

  UInt nb_element = mesh.getNbElement(type,ghost_type);


  Array<Real> * jacobians_tmp;
  if(!jacobians.exists(type, ghost_type))
    jacobians_tmp = &jacobians.alloc(nb_element*nb_quadrature_points,
                                     1,
                                     type,
                                     ghost_type);
  else {
    jacobians_tmp = &jacobians(type, ghost_type);
    jacobians_tmp->resize(nb_element*nb_quadrature_points);
  }

  Array<Real>::vector_iterator jacobians_it =
    jacobians_tmp->begin_reinterpret(nb_quadrature_points, nb_element);

  Vector<Real> weights = GaussIntegrationElement<type>::getWeights();

  Array<Real> x_el(0, spatial_dimension * nb_nodes_per_element);
  FEEngine::extractNodalToElementField(mesh, nodes, x_el, type, ghost_type);

  Array<Real>::const_matrix_iterator x_it = x_el.begin(spatial_dimension,
						       nb_nodes_per_element);

  UInt nb_nodes_per_subelement = nb_nodes_per_element / 2;
  Matrix<Real> x(spatial_dimension, nb_nodes_per_subelement);

  //  Matrix<Real> local_coord(spatial_dimension, nb_nodes_per_element);
  for (UInt elem = 0; elem < nb_element; ++elem, ++jacobians_it, ++x_it) {

    for (UInt s = 0; s < spatial_dimension; ++s)
      for (UInt n = 0; n < nb_nodes_per_subelement; ++n)
	x(s, n) = ((*x_it)(s, n) + (*x_it)(s, n + nb_nodes_per_subelement))*.5;

    Vector<Real> & J = *jacobians_it;

    if (type == _cohesive_1d_2)
      J(0) = 1;
    else
      computeJacobianOnQuadPointsByElement<type>(x, J);

    J *= weights;
  }

  AKANTU_DEBUG_OUT();
}
#endif

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <ElementType type>
void IntegratorGauss<kind>::integrate(const Array<Real> & in_f,
				      Array<Real> &intf,
				      UInt nb_degree_of_freedom,
				      const GhostType & ghost_type,
				      const Array<UInt> & filter_elements) const {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_ASSERT(jacobians.exists(type, ghost_type),
		      "No jacobians for the type "
		      << jacobians.printType(type, ghost_type));

  UInt nb_points = GaussIntegrationElement<type>::getNbQuadraturePoints();

  const Array<Real> & jac_loc = jacobians(type, ghost_type);

  Array<Real>::const_matrix_iterator J_it;
  Array<Real>::matrix_iterator inte_it;
  Array<Real>::const_matrix_iterator f_it;

  UInt nb_element;
  Array<Real> * filtered_J = NULL;
  if(filter_elements != empty_filter) {
    nb_element = filter_elements.getSize();
    filtered_J = new Array<Real>(0, jac_loc.getNbComponent());
    FEEngine::filterElementalData(mesh, jac_loc, *filtered_J, type, ghost_type, filter_elements);
    const Array<Real> & cfiltered_J = *filtered_J; // \todo temporary patch
    J_it = cfiltered_J.begin_reinterpret(nb_points, 1, nb_element);
  } else {
    nb_element = mesh.getNbElement(type,ghost_type);
    J_it = jac_loc.begin_reinterpret(nb_points, 1, nb_element);
  }

  intf.resize(nb_element);

  f_it    = in_f.begin_reinterpret(nb_degree_of_freedom, nb_points, nb_element);
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
template <ElementKind kind>
template <ElementType type>
Real IntegratorGauss<kind>::integrate(const Array<Real> & in_f,
				      const GhostType & ghost_type,
				      const Array<UInt> & filter_elements) const {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_ASSERT(jacobians.exists(type, ghost_type),
		      "No jacobians for the type "
		      << jacobians.printType(type, ghost_type));

  Array<Real> intfv(0, 1);
  integrate<type>(in_f, intfv, 1, ghost_type, filter_elements);

  UInt nb_values = intfv.getSize();
  if(nb_values == 0) return 0.;

  UInt nb_values_to_sum = nb_values >> 1;

  std::sort(intfv.begin(), intfv.end());


  // as long as the half is not empty
  while(nb_values_to_sum) {
    UInt remaining = (nb_values - 2*nb_values_to_sum);
    if(remaining)  intfv(nb_values - 2) += intfv(nb_values - 1);

    // sum to consecutive values and store the sum in the first half
    for (UInt i = 0; i < nb_values_to_sum; ++i) {
      intfv(i) = intfv(2*i) + intfv(2*i + 1);
    }

    nb_values = nb_values_to_sum;
    nb_values_to_sum >>= 1;
  }

  AKANTU_DEBUG_OUT();
  return intfv(0);
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind>
template <ElementType type>
void IntegratorGauss<kind>::integrateOnQuadraturePoints(const Array<Real> & in_f,
							Array<Real> &intf,
							UInt nb_degree_of_freedom,
							const GhostType & ghost_type,
							const Array<UInt> & filter_elements) const {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_ASSERT(jacobians.exists(type, ghost_type),
		      "No jacobians for the type "
		      << jacobians.printType(type, ghost_type));

  UInt nb_element;
  UInt nb_points = GaussIntegrationElement<type>::getNbQuadraturePoints();

  const Array<Real> & jac_loc = jacobians(type, ghost_type);

  Array<Real>::const_scalar_iterator J_it;
  Array<Real>::vector_iterator inte_it;
  Array<Real>::const_vector_iterator f_it;

  Array<Real> * filtered_J = NULL;
  if(filter_elements != empty_filter) {
    nb_element = filter_elements.getSize();
    filtered_J = new Array<Real>(0, jac_loc.getNbComponent());
    FEEngine::filterElementalData(mesh, jac_loc, *filtered_J, type, ghost_type, filter_elements);
    J_it = filtered_J->begin();
  } else {
    nb_element = mesh.getNbElement(type,ghost_type);
    J_it = jac_loc.begin();
  }

  intf.resize(nb_element*nb_points);

  f_it    = in_f.begin(nb_degree_of_freedom);
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

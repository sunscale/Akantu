/**
 * @file   fe_engine_template_tmpl_field.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Wed Aug 09 2017
 * @date last modification: Thu Dec 07 2017
 *
 * @brief  implementation of the assemble field s functions
 *
 * @section LICENSE
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
#include "fe_engine_template.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_FE_ENGINE_TEMPLATE_TMPL_FIELD_HH__
#define __AKANTU_FE_ENGINE_TEMPLATE_TMPL_FIELD_HH__

namespace akantu {

/* -------------------------------------------------------------------------- */
/* Matrix lumping functions                                                   */
/* -------------------------------------------------------------------------- */
namespace fe_engine {
  namespace details {
    namespace {
      template <class Functor>
      void fillField(const Functor & field_funct, Array<Real> & field,
                     UInt nb_element, UInt nb_integration_points,
                     const ElementType & type, const GhostType & ghost_type) {
        UInt nb_degree_of_freedom = field.getNbComponent();
        field.resize(nb_integration_points * nb_element);

        auto field_it = field.begin_reinterpret(
            nb_degree_of_freedom, nb_integration_points, nb_element);

        Element el{type, 0, ghost_type};
        for (; el.element < nb_element; ++el.element, ++field_it) {
          field_funct(*field_it, el);
        }
      }
    } // namespace
  }   // namespace details
} // namespace fe_engine

/**
 * Helper class to be able to write a partial specialization on the element kind
 */
namespace fe_engine {
  namespace details {
    template <ElementKind kind> struct AssembleLumpedTemplateHelper {};

#define ASSEMBLE_LUMPED(type)                                                  \
  fem.template assembleFieldLumped<Functor, type>(field_funct, lumped, dof_id, \
                                                  dof_manager, ghost_type)

#define AKANTU_SPECIALIZE_ASSEMBLE_HELPER(kind)                                \
  template <> struct AssembleLumpedTemplateHelper<kind> {                      \
    template <template <ElementKind, class> class I,                           \
              template <ElementKind> class S, ElementKind k, class IOF,        \
              class Functor>                                                   \
    static void call(const FEEngineTemplate<I, S, k, IOF> & fem,               \
                     const Functor & field_funct, const ID & lumped,           \
                     const ID & dof_id, DOFManager & dof_manager,              \
                     ElementType type, const GhostType & ghost_type) {         \
      AKANTU_BOOST_KIND_ELEMENT_SWITCH(ASSEMBLE_LUMPED, kind);                 \
    }                                                                          \
  };

    AKANTU_BOOST_ALL_KIND(AKANTU_SPECIALIZE_ASSEMBLE_HELPER)

#undef AKANTU_SPECIALIZE_ASSEMBLE_HELPER
#undef ASSEMBLE_LUMPED
  } // namespace details
} // namespace fe_engine

/* -------------------------------------------------------------------------- */
template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IOF>
template <class Functor>
void FEEngineTemplate<I, S, kind, IOF>::assembleFieldLumped(
    const Functor & field_funct, const ID & matrix_id, const ID & dof_id,
    DOFManager & dof_manager, ElementType type,
    const GhostType & ghost_type) const {
  AKANTU_DEBUG_IN();

  fe_engine::details::AssembleLumpedTemplateHelper<kind>::call(
      *this, field_funct, matrix_id, dof_id, dof_manager, type, ghost_type);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
template <class Functor, ElementType type>
void FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::assembleFieldLumped(
    const Functor & field_funct, const ID & matrix_id, const ID & dof_id,
    DOFManager & dof_manager, const GhostType & ghost_type) const {
  AKANTU_DEBUG_IN();

  UInt nb_degree_of_freedom = dof_manager.getDOFs(dof_id).getNbComponent();
  UInt nb_element = mesh.getNbElement(type, ghost_type);
  UInt nb_integration_points = this->getNbIntegrationPoints(type);

  Array<Real> field(0, nb_degree_of_freedom);
  fe_engine::details::fillField(field_funct, field, nb_element,
                                nb_integration_points, type, ghost_type);

  switch (type) {
  case _triangle_6:
  case _quadrangle_8:
  case _tetrahedron_10:
  case _hexahedron_20:
  case _pentahedron_15:
    this->template assembleLumpedDiagonalScaling<type>(field, matrix_id, dof_id,
                                                       dof_manager, ghost_type);
    break;
  default:
    this->template assembleLumpedRowSum<type>(field, matrix_id, dof_id,
                                              dof_manager, ghost_type);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/**
 * @f$ \tilde{M}_{i} = \sum_j M_{ij} = \sum_j \int \rho \varphi_i \varphi_j dV =
 * \int \rho \varphi_i dV @f$
 */
template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
template <ElementType type>
void FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::
    assembleLumpedRowSum(const Array<Real> & field, const ID & matrix_id,
                         const ID & dof_id, DOFManager & dof_manager,
                         const GhostType & ghost_type) const {
  AKANTU_DEBUG_IN();

  UInt shapes_size = ElementClass<type>::getShapeSize();
  UInt nb_degree_of_freedom = field.getNbComponent();

  Array<Real> * field_times_shapes =
      new Array<Real>(0, shapes_size * nb_degree_of_freedom);

  shape_functions.template computeNtb<type>(field, *field_times_shapes,
                                            ghost_type);

  UInt nb_element = mesh.getNbElement(type, ghost_type);
  Array<Real> * int_field_times_shapes = new Array<Real>(
      nb_element, shapes_size * nb_degree_of_freedom, "inte_rho_x_shapes");

  integrator.template integrate<type>(
      *field_times_shapes, *int_field_times_shapes,
      nb_degree_of_freedom * shapes_size, ghost_type, empty_filter);

  delete field_times_shapes;

  dof_manager.assembleElementalArrayToLumpedMatrix(
      dof_id, *int_field_times_shapes, matrix_id, type, ghost_type);

  delete int_field_times_shapes;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/**
 * @f$ \tilde{M}_{i} = c * M_{ii} = \int_{V_e} \rho dV @f$
 */
template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
template <ElementType type>
void FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::
    assembleLumpedDiagonalScaling(const Array<Real> & field,
                                  const ID & matrix_id, const ID & dof_id,
                                  DOFManager & dof_manager,
                                  const GhostType & ghost_type) const {
  AKANTU_DEBUG_IN();

  const ElementType & type_p1 = ElementClass<type>::getP1ElementType();
  UInt nb_nodes_per_element_p1 = Mesh::getNbNodesPerElement(type_p1);
  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
  UInt nb_degree_of_freedom = field.getNbComponent();
  UInt nb_element = mesh.getNbElement(type, ghost_type);

  Vector<Real> nodal_factor(nb_nodes_per_element);

#define ASSIGN_WEIGHT_TO_NODES(corner, mid)                                    \
  {                                                                            \
    for (UInt n = 0; n < nb_nodes_per_element_p1; n++)                         \
      nodal_factor(n) = corner;                                                \
    for (UInt n = nb_nodes_per_element_p1; n < nb_nodes_per_element; n++)      \
      nodal_factor(n) = mid;                                                   \
  }

  if (type == _triangle_6)
    ASSIGN_WEIGHT_TO_NODES(1. / 12., 1. / 4.);
  if (type == _tetrahedron_10)
    ASSIGN_WEIGHT_TO_NODES(1. / 32., 7. / 48.);
  if (type == _quadrangle_8)
    ASSIGN_WEIGHT_TO_NODES(
        3. / 76.,
        16. / 76.); /** coeff. derived by scaling
                     * the diagonal terms of the corresponding
                     * consistent mass computed with 3x3 gauss points;
                     * coeff. are (1./36., 8./36.) for 2x2 gauss points */
  if (type == _hexahedron_20)
    ASSIGN_WEIGHT_TO_NODES(
        7. / 248., 16. / 248.); /** coeff. derived by scaling
                                 * the diagonal terms of the corresponding
                                 * consistent mass computed with 3x3x3 gauss
                                 * points; coeff. are (1./40.,
                                 * 1./15.) for 2x2x2 gauss points */
  if (type == _pentahedron_15) {
    // coefficients derived by scaling the diagonal terms of the corresponding
    // consistent mass computed with 8 gauss points;
    for (UInt n = 0; n < nb_nodes_per_element_p1; n++)
      nodal_factor(n) = 51. / 2358.;

    Real mid_triangle = 192. / 2358.;
    Real mid_quadrangle = 300. / 2358.;

    nodal_factor(6) = mid_triangle;
    nodal_factor(7) = mid_triangle;
    nodal_factor(8) = mid_triangle;
    nodal_factor(9) = mid_quadrangle;
    nodal_factor(10) = mid_quadrangle;
    nodal_factor(11) = mid_quadrangle;
    nodal_factor(12) = mid_triangle;
    nodal_factor(13) = mid_triangle;
    nodal_factor(14) = mid_triangle;
  }

  if (nb_element == 0) {
    AKANTU_DEBUG_OUT();
    return;
  }

#undef ASSIGN_WEIGHT_TO_NODES
  /// compute @f$ \int \rho dV = \rho V @f$ for each element
  auto int_field = std::make_unique<Array<Real>>(
      field.size(), nb_degree_of_freedom, "inte_rho_x");
  integrator.template integrate<type>(field, *int_field, nb_degree_of_freedom,
                                      ghost_type, empty_filter);

  /// distribute the mass of the element to the nodes
  auto lumped_per_node = std::make_unique<Array<Real>>(
      nb_element, nb_degree_of_freedom * nb_nodes_per_element, "mass_per_node");
  auto int_field_it = int_field->begin(nb_degree_of_freedom);
  auto lumped_per_node_it =
      lumped_per_node->begin(nb_degree_of_freedom, nb_nodes_per_element);

  for (UInt e = 0; e < nb_element; ++e) {
    for (UInt n = 0; n < nb_nodes_per_element; ++n) {
      Vector<Real> l = (*lumped_per_node_it)(n);
      l = *int_field_it;
      l *= nodal_factor(n);
    }
    ++int_field_it;
    ++lumped_per_node_it;
  }

  dof_manager.assembleElementalArrayToLumpedMatrix(dof_id, *lumped_per_node,
                                                   matrix_id, type, ghost_type);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/**
 * Helper class to be able to write a partial specialization on the element kind
 */
namespace fe_engine {
  namespace details {
    template <ElementKind kind> struct AssembleFieldMatrixHelper {};

#define ASSEMBLE_MATRIX(type)                                                  \
  fem.template assembleFieldMatrix<Functor, type>(                             \
      field_funct, matrix_id, dof_id, dof_manager, ghost_type)

#define AKANTU_SPECIALIZE_ASSEMBLE_FIELD_MATRIX_HELPER(kind)                   \
  template <> struct AssembleFieldMatrixHelper<kind> {                         \
    template <template <ElementKind, class> class I,                           \
              template <ElementKind> class S, ElementKind k, class IOF,        \
              class Functor>                                                   \
    static void call(const FEEngineTemplate<I, S, k, IOF> & fem,               \
                     const Functor & field_funct, const ID & matrix_id,        \
                     const ID & dof_id, DOFManager & dof_manager,              \
                     ElementType type, const GhostType & ghost_type) {         \
      AKANTU_BOOST_KIND_ELEMENT_SWITCH(ASSEMBLE_MATRIX, kind);                 \
    }                                                                          \
  };

    AKANTU_BOOST_ALL_KIND(AKANTU_SPECIALIZE_ASSEMBLE_FIELD_MATRIX_HELPER)

#undef AKANTU_SPECIALIZE_ASSEMBLE_FIELD_MATRIX_HELPER
#undef ASSEMBLE_MATRIX
  } // namespace details
} // namespace fe_engine

/* -------------------------------------------------------------------------- */
template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IOF>
template <class Functor>
void FEEngineTemplate<I, S, kind, IOF>::assembleFieldMatrix(
    const Functor & field_funct, const ID & matrix_id, const ID & dof_id,
    DOFManager & dof_manager, ElementType type,
    const GhostType & ghost_type) const {
  AKANTU_DEBUG_IN();
  fe_engine::details::AssembleFieldMatrixHelper<kind>::template call(
      *this, field_funct, matrix_id, dof_id, dof_manager, type, ghost_type);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/**
 * @f$ \tilde{M}_{i} = \sum_j M_{ij} = \sum_j \int \rho \varphi_i \varphi_j dV =
 * \int \rho \varphi_i dV @f$
 */
template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
template <class Functor, ElementType type>
void FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::assembleFieldMatrix(
    const Functor & field_funct, const ID & matrix_id, const ID & dof_id,
    DOFManager & dof_manager, const GhostType & ghost_type) const {
  AKANTU_DEBUG_IN();

  UInt shapes_size = ElementClass<type>::getShapeSize();
  UInt nb_degree_of_freedom = dof_manager.getDOFs(dof_id).getNbComponent();
  UInt lmat_size = nb_degree_of_freedom * shapes_size;
  UInt nb_element = mesh.getNbElement(type, ghost_type);

  // \int N * N  so degree 2 * degree of N
  const UInt polynomial_degree =
      2 * ElementClassProperty<type>::polynomial_degree;

  // getting the integration points
  Matrix<Real> integration_points =
      integrator.template getIntegrationPoints<type, polynomial_degree>();

  UInt nb_integration_points = integration_points.cols();
  UInt vect_size = nb_integration_points * nb_element;

  // getting the shapes on the integration points
  Array<Real> shapes(0, shapes_size);
  shape_functions.template computeShapesOnIntegrationPoints<type>(
      mesh.getNodes(), integration_points, shapes, ghost_type);

  // Extending the shape functions
  /// \TODO move this in the shape functions as Voigt format shapes to have the
  /// code in common with the structural elements
  Array<Real> modified_shapes(vect_size, lmat_size * nb_degree_of_freedom, 0.);
  Array<Real> local_mat(vect_size, lmat_size * lmat_size);
  auto mshapes_it = modified_shapes.begin(nb_degree_of_freedom, lmat_size);
  auto shapes_it = shapes.begin(shapes_size);

  for (UInt q = 0; q < vect_size; ++q, ++mshapes_it, ++shapes_it) {
    for (UInt d = 0; d < nb_degree_of_freedom; ++d) {
      for (UInt s = 0; s < shapes_size; ++s) {
        (*mshapes_it)(d, s * nb_degree_of_freedom + d) = (*shapes_it)(s);
      }
    }
  }

  // getting the value to assemble on the integration points
  Array<Real> field(vect_size, nb_degree_of_freedom);
  fe_engine::details::fillField(field_funct, field, nb_element,
                                nb_integration_points, type, ghost_type);

  // computing \rho * N
  mshapes_it = modified_shapes.begin(nb_degree_of_freedom, lmat_size);
  auto lmat = local_mat.begin(lmat_size, lmat_size);
  auto field_it = field.begin_reinterpret(nb_degree_of_freedom, field.size());

  for (UInt q = 0; q < vect_size; ++q, ++lmat, ++mshapes_it, ++field_it) {
    const auto & rho = *field_it;
    const auto & N = *mshapes_it;
    auto & mat = *lmat;

    Matrix<Real> Nt = N.transpose();
    for (UInt d = 0; d < Nt.cols(); ++d) {
      Nt(d) *= rho(d);
    }

    mat.template mul<false, false>(Nt, N);
  }

  // integrate the elemental values
  Array<Real> int_field_times_shapes(nb_element, lmat_size * lmat_size,
                                     "inte_rho_x_shapes");
  this->integrator.template integrate<type, polynomial_degree>(
      local_mat, int_field_times_shapes, lmat_size * lmat_size, ghost_type);

  // assemble the elemental values to the matrix
  dof_manager.assembleElementalMatricesToMatrix(
      matrix_id, dof_id, int_field_times_shapes, type, ghost_type);

  AKANTU_DEBUG_OUT();
}

} // namespace akantu

#endif /* __AKANTU_FE_ENGINE_TEMPLATE_TMPL_FIELD_HH__ */

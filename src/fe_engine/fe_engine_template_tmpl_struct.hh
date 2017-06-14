/**
 * @file   fe_engine_template_tmpl_struct.hh
 *
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @author Sébastien Hartmann <sebastien.hartmann@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Mon Jul 07 2014
 * @date last modification: Thu Oct 15 2015
 *
 * @brief  Template implementation of FEEngineTemplate for Structural Element
 * Kinds
 *
 * @section LICENSE
 *
 * Copyright  (©)  2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de Lausanne)
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

/* -------------------------------------------------------------------------- */
#include "shape_linked.hh"

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template <>
inline const Array<Real> &
FEEngineTemplate<IntegratorGauss, ShapeLinked, _ek_structural,
                 DefaultIntegrationOrderFunctor>::
    getShapesDerivatives(const ElementType & type, const GhostType & ghost_type,
                         UInt id) const {

  AKANTU_DEBUG_IN();
  const Array<Real> * ret = NULL;

#define GET_SHAPES(type)                                                       \
  ret = &(shape_functions.getShapesDerivatives(type, ghost_type, id));

  AKANTU_BOOST_STRUCTURAL_ELEMENT_SWITCH(GET_SHAPES);
#undef GET_SHAPES

  AKANTU_DEBUG_OUT();
  return *ret;
}

/* -------------------------------------------------------------------------- */
template <>
inline const Array<Real> & FEEngineTemplate<
    IntegratorGauss, ShapeLinked, _ek_structural,
    DefaultIntegrationOrderFunctor>::getShapes(const ElementType & type,
                                               const GhostType & ghost_type,
                                               UInt id) const {
  AKANTU_DEBUG_IN();
  const Array<Real> * ret = NULL;

#define GET_SHAPES(type)                                                       \
  ret = &(shape_functions.getShapes(type, ghost_type, id));

  AKANTU_BOOST_STRUCTURAL_ELEMENT_SWITCH(GET_SHAPES);
#undef GET_SHAPES

  AKANTU_DEBUG_OUT();
  return *ret;
}

/* -------------------------------------------------------------------------- */
template <ElementKind kind, typename = void>
struct AssembleFieldMatrixStructHelper {};

template <ElementKind kind>
struct AssembleFieldMatrixStructHelper<
    kind, typename std::enable_if<kind == _ek_structural>::type> {
  template <template <ElementKind, class> class I,
            template <ElementKind> class S, ElementKind k, class IOF>
  static void call(const FEEngineTemplate<I, S, k, IOF> & fem,
                   const Array<Real> & field_1, UInt nb_degree_of_freedom,
                   SparseMatrix & M, Array<Real> * n,
                   ElementTypeMapArray<Real> & rotation_mat,
                   const ElementType & type, const GhostType & ghost_type) {
#define ASSEMBLE_MATRIX(type)                                                  \
  fem.template assembleFieldMatrix<type>(field_1, nb_degree_of_freedom, M, n,  \
                                         rotation_mat, ghost_type)

    AKANTU_BOOST_KIND_ELEMENT_SWITCH(ASSEMBLE_MATRIX, _ek_structural);
#undef ASSEMBLE_MATRIX
  }
};

template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
inline void
FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::assembleFieldMatrix(
    const Array<Real> & field_1, UInt nb_degree_of_freedom, SparseMatrix & M,
    Array<Real> * n, ElementTypeMapArray<Real> & rotation_mat,
    const ElementType & type, const GhostType & ghost_type) const {
  AKANTU_DEBUG_IN();

  AssembleFieldMatrixStructHelper<kind>::template call(
      *this, field_1, nb_degree_of_freedom, M, n, rotation_mat, type,
      ghost_type);

  AKANTU_DEBUG_OUT();
}
/* -------------------------------------------------------------------------- */
template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
inline void
FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::computeShapesMatrix(
    const ElementType &, UInt, UInt, Array<Real> *,
    __attribute__((unused)) UInt, UInt, UInt, const bool,
    const GhostType &) const {

  AKANTU_DEBUG_TO_IMPLEMENT();
}

/* -------------------------------------------------------------------------- */
template <>
inline void FEEngineTemplate<IntegratorGauss, ShapeLinked, _ek_structural,
                             DefaultIntegrationOrderFunctor>::
    computeShapesMatrix(const ElementType & type, UInt nb_degree_of_freedom,
                        UInt nb_nodes_per_element, Array<Real> * n, UInt id,
                        UInt degree_to_interpolate, UInt degree_interpolated,
                        const bool sign, // true +, false -
                        const GhostType & ghost_type) const {
  AKANTU_DEBUG_IN();

  UInt nb_element = mesh.getNbElement(type);
  UInt nb_quadrature_points = getNbIntegrationPoints(type);

  UInt nt_n_field_size = nb_degree_of_freedom * nb_nodes_per_element;
  UInt n_size = n->getNbComponent() / nt_n_field_size;

  Array<Real>::const_vector_iterator shape =
      getShapes(type, ghost_type, id).begin(nb_nodes_per_element);
  Array<Real>::matrix_iterator N_it = n->begin(n_size, nt_n_field_size);

  int c;
  if (sign == true) {
    c = 1;
  } else {
    c = -1;
  }

  UInt line = degree_interpolated;
  UInt coll = degree_to_interpolate;

  for (UInt e = 0; e < nb_element; ++e) {
    for (UInt q = 0; q < nb_quadrature_points; ++q, ++N_it, ++shape) {
      const Vector<Real> & shapes = *shape;
      Matrix<Real> & N = *N_it;
      N(line, coll) = shapes(0) * c;
      N(line, coll + nb_degree_of_freedom) = shapes(1) * c;
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <>
template <ElementType type>
inline void FEEngineTemplate<IntegratorGauss, ShapeLinked, _ek_structural,
                             DefaultIntegrationOrderFunctor>::
    assembleFieldMatrix(const Array<Real> & field_1, UInt nb_degree_of_freedom,
                        SparseMatrix & M, Array<Real> * n,
                        ElementTypeMapArray<Real> & rotation_mat,
                        const GhostType & ghost_type) const {
  AKANTU_DEBUG_IN();

  UInt nb_element = mesh.getNbElement(type);
  UInt nb_nodes_per_element = mesh.getNbNodesPerElement(type);
  UInt nb_quadrature_points = getNbIntegrationPoints(type);

  UInt nt_n_field_size = nb_degree_of_freedom * nb_nodes_per_element;
  UInt n_size = n->getNbComponent() / nt_n_field_size;

  Array<Real> * nt_n_field = new Array<Real>(
      nb_element * nb_quadrature_points, // nt_n_size * nt_n_size, nb_elem *
                                         // nb_quad_points?
      nt_n_field_size * nt_n_field_size, "NT*N*field");
  Array<Real> * nt = new Array<Real>(nb_element * nb_quadrature_points,
                                     n_size * nt_n_field_size, "N*T");
  Array<Real> t = rotation_mat(type);
  nt_n_field->clear();
  nt->clear();

  Array<Real>::matrix_iterator N = n->begin(n_size, nt_n_field_size);
  Array<Real>::matrix_iterator Nt_N_field =
      nt_n_field->begin(nt_n_field_size, nt_n_field_size);
  Array<Real>::matrix_iterator T =
      rotation_mat(type).begin(nt_n_field_size, nt_n_field_size);
  Array<Real>::matrix_iterator NT = nt->begin(n_size, nt_n_field_size);
  Real * field_val = field_1.storage();

  for (UInt e = 0; e < nb_element; ++e, ++T) {
    for (UInt q = 0; q < nb_quadrature_points;
         ++q, ++N, ++NT, ++Nt_N_field, /*++T,*/ ++field_val) {
      NT->mul<false, false>(*N, *T);
      Nt_N_field->mul<true, false>(*NT, *NT, *field_val);
    }
  }

  Array<Real> * int_nt_n_field = new Array<Real>(
      nb_element, nt_n_field_size * nt_n_field_size, "NT*N*field");
  int_nt_n_field->clear();

  integrate(*nt_n_field, *int_nt_n_field, nt_n_field_size * nt_n_field_size,
            type);
  //  integrate(*nt_n_field, *int_nt_n_field, nb_degree_of_freedom, type);

  assembleMatrix(*int_nt_n_field, M, nb_degree_of_freedom, type);

  delete nt;
  delete nt_n_field;
  delete int_nt_n_field;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
template <ElementType type>
inline void
FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::assembleFieldMatrix(
    const Array<Real> & field_1, UInt nb_degree_of_freedom, SparseMatrix & M,
    Array<Real> * n, ElementTypeMapArray<Real> & rotation_mat,
    const GhostType & ghost_type) const {
  AKANTU_DEBUG_TO_IMPLEMENT();
}

__END_AKANTU__

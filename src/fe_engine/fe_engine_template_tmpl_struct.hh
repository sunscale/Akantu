/**
 * @file   fe_engine_template_tmpl_struct.hh
 *
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @author Sébastien Hartmann <sebastien.hartmann@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Mon Jul 07 2014
 * @date last modification: Tue Feb 20 2018
 *
 * @brief  Template implementation of FEEngineTemplate for Structural Element
 * Kinds
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

/* -------------------------------------------------------------------------- */
#include "shape_structural.hh"

namespace akantu {

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

// template <template <ElementKind, class> class I, template <ElementKind> class
// S,
//           ElementKind kind, class IntegrationOrderFunctor>
// inline void
// FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::assembleFieldMatrix(
//     const Array<Real> & field_1, UInt nb_degree_of_freedom, SparseMatrix & M,
//     Array<Real> * n, ElementTypeMapArray<Real> & rotation_mat,
//     const ElementType & type, const GhostType & ghost_type) const {
//   AKANTU_DEBUG_IN();

//   AssembleFieldMatrixStructHelper<kind>::template call(
//       *this, field_1, nb_degree_of_freedom, M, n, rotation_mat, type,
//       ghost_type);

//   AKANTU_DEBUG_OUT();
// }
// /* --------------------------------------------------------------------------
// */ template <template <ElementKind, class> class I, template <ElementKind>
// class S,
//           ElementKind kind, class IntegrationOrderFunctor>
// inline void
// FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::computeShapesMatrix(
//     const ElementType &, UInt, UInt, Array<Real> *, UInt, UInt, UInt,
//     const bool, const GhostType &) const {

//   AKANTU_TO_IMPLEMENT();
// }

/* -------------------------------------------------------------------------- */
template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
template <ElementType type>
inline void
FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::assembleFieldMatrix(
    const Array<Real> &, UInt, SparseMatrix &, Array<Real> *,
    ElementTypeMapArray<Real> &, const GhostType &) const {
  AKANTU_TO_IMPLEMENT();
}

} // namespace akantu

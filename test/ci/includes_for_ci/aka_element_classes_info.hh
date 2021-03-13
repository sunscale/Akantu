/**
 * @file   aka_element_classes_info.hh.in
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Sun Jul 19 2015
 * @date last modification: Tue Feb 20 2018
 *
 * @brief  Declaration of the enums for the element classes
 *
 *
 * Copyright (©) 2015-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "aka_safe_enum.hh"
/* -------------------------------------------------------------------------- */
#include <boost/preprocessor.hpp>
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_AKA_ELEMENT_CLASSES_INFO_HH_
#define AKANTU_AKA_ELEMENT_CLASSES_INFO_HH_

namespace akantu {

/* -------------------------------------------------------------------------- */
/* Element Types                                                              */
/* -------------------------------------------------------------------------- */

/// @enum ElementType type of elements
enum ElementType {
  _not_defined,
  _cohesive_1d_2,
  _cohesive_2d_4,
  _cohesive_2d_6,
  _cohesive_3d_12,
  _cohesive_3d_16,
  _cohesive_3d_6,
  _cohesive_3d_8,
  _point_1,
  _segment_2,
  _segment_3,
  _triangle_3,
  _triangle_6,
  _quadrangle_4,
  _quadrangle_8,
  _tetrahedron_4,
  _tetrahedron_10,
  _pentahedron_6,
  _pentahedron_15,
  _hexahedron_8,
  _hexahedron_20,
  _bernoulli_beam_2,
  _bernoulli_beam_3,
  _discrete_kirchhoff_triangle_18,
  _max_element_type
};


#define AKANTU_ek_cohesive_ELEMENT_TYPE	\
  (_cohesive_1d_2)			\
  (_cohesive_2d_4)			\
  (_cohesive_2d_6)			\
  (_cohesive_3d_12)			\
  (_cohesive_3d_16)			\
  (_cohesive_3d_6)			\
  (_cohesive_3d_8)

#define AKANTU_ek_regular_ELEMENT_TYPE	\
  (_point_1)				\
  (_segment_2)				\
  (_segment_3)				\
  (_triangle_3)				\
  (_triangle_6)				\
  (_quadrangle_4)			\
  (_quadrangle_8)			\
  (_tetrahedron_4)			\
  (_tetrahedron_10)			\
  (_pentahedron_6)			\
  (_pentahedron_15)			\
  (_hexahedron_8)			\
  (_hexahedron_20)

#define AKANTU_ek_structural_ELEMENT_TYPE	\
  (_bernoulli_beam_2)			\
  (_bernoulli_beam_3)			\
  (_discrete_kirchhoff_triangle_18)


#define AKANTU_ALL_ELEMENT_TYPE		\
  AKANTU_ek_cohesive_ELEMENT_TYPE	\
  AKANTU_ek_regular_ELEMENT_TYPE	\
  AKANTU_ek_structural_ELEMENT_TYPE

/* -------------------------------------------------------------------------- */
/* Element Kinds                                                              */
/* -------------------------------------------------------------------------- */

#define AKANTU_COHESIVE_KIND	(_ek_cohesive)
#define AKANTU_REGULAR_KIND	(_ek_regular)
#define AKANTU_STRUCTURAL_KIND	(_ek_structural)

#define AKANTU_ELEMENT_KIND		\
  AKANTU_COHESIVE_KIND			\
  AKANTU_REGULAR_KIND			\
  AKANTU_STRUCTURAL_KIND

enum ElementKind {
  BOOST_PP_SEQ_ENUM(AKANTU_ELEMENT_KIND),
  _ek_not_defined
};


/* -------------------------------------------------------------------------- */
struct ElementKind_def {
  using type = ElementKind;
  static const type _begin_ = BOOST_PP_SEQ_HEAD(AKANTU_ELEMENT_KIND);
  static const type _end_   = _ek_not_defined;
};

using element_kind_t = safe_enum<ElementKind_def> ;

/* -------------------------------------------------------------------------- */
/// @enum GeometricalType type of element potentially contained in a Mesh
enum GeometricalType {
  _gt_cohesive_1d_2,
  _gt_cohesive_2d_4,
  _gt_cohesive_2d_6,
  _gt_cohesive_3d_12,
  _gt_cohesive_3d_16,
  _gt_cohesive_3d_6,
  _gt_cohesive_3d_8,
  _gt_point,
  _gt_segment_2,
  _gt_segment_3,
  _gt_triangle_3,
  _gt_triangle_6,
  _gt_quadrangle_4,
  _gt_quadrangle_8,
  _gt_tetrahedron_4,
  _gt_tetrahedron_10,
  _gt_hexahedron_8,
  _gt_hexahedron_20,
  _gt_pentahedron_6,
  _gt_pentahedron_15,
  _gt_not_defined
};

/* -------------------------------------------------------------------------- */
/* Interpolation Types                                                        */
/* -------------------------------------------------------------------------- */
#define AKANTU_INTERPOLATION_TYPES		\
  (_itp_lagrange_point_1)			\
  (_itp_lagrange_segment_2)			\
  (_itp_lagrange_segment_3)			\
  (_itp_lagrange_triangle_3)			\
  (_itp_lagrange_triangle_6)			\
  (_itp_lagrange_quadrangle_4)			\
  (_itp_serendip_quadrangle_8)			\
  (_itp_lagrange_tetrahedron_4)			\
  (_itp_lagrange_tetrahedron_10)		\
  (_itp_lagrange_hexahedron_8)			\
  (_itp_serendip_hexahedron_20)			\
  (_itp_lagrange_pentahedron_6)			\
  (_itp_lagrange_pentahedron_15)\
  (_itp_hermite_2)			\
  (_itp_bernoulli_beam_2)			\
  (_itp_bernoulli_beam_3)			\
  (_itp_discrete_kirchhoff_triangle_6)		\
  (_itp_discrete_kirchhoff_triangle_18)

/// @enum InterpolationType type of elements
enum InterpolationType {
  BOOST_PP_SEQ_ENUM(AKANTU_INTERPOLATION_TYPES),
  _itp_not_defined
};

/* -------------------------------------------------------------------------- */
/* Some sub types less probable to change                                     */
/* -------------------------------------------------------------------------- */
/// @enum GeometricalShapeType types of shapes to define the contains
/// function in the element classes
enum GeometricalShapeType {
  _gst_point,
  _gst_triangle,
  _gst_square,
  _gst_prism,
  _gst_not_defined
};

/* -------------------------------------------------------------------------- */
/// @enum GaussIntegrationType classes of types using common
/// description of the gauss point position and weights
enum GaussIntegrationType {
  _git_point,
  _git_segment,
  _git_triangle,
  _git_tetrahedron,
  _git_pentahedron,
  _git_not_defined
};

/* -------------------------------------------------------------------------- */
/// @enum InterpolationKind the family of interpolation types
enum InterpolationKind {
  _itk_lagrangian,
  _itk_structural,
  _itk_not_defined
};

/* -------------------------------------------------------------------------- */
// BOOST PART: TOUCH ONLY IF YOU KNOW WHAT YOU ARE DOING
#define AKANTU_BOOST_CASE_MACRO(r, macro, _type)                               \
  case _type: {                                                                \
    macro(_type);                                                              \
    break;                                                                     \
  }

#define AKANTU_BOOST_LIST_SWITCH(macro1, list1, var)                           \
  do {                                                                         \
    switch (var) {                                                             \
      BOOST_PP_SEQ_FOR_EACH(AKANTU_BOOST_CASE_MACRO, macro1, list1)            \
    default: {                                                                 \
      AKANTU_ERROR("Type ("                                                    \
                   << var /* NOLINT */ << ") not handled by this function");   \
    }                                                                          \
    }                                                                          \
  } while (0)

#define AKANTU_BOOST_LIST_SWITCH_NO_DEFAULT(macro1, list1, var)                \
  do {                                                                         \
    switch (var) {                                                             \
      BOOST_PP_SEQ_FOR_EACH(AKANTU_BOOST_CASE_MACRO, macro1, list1)            \
    case _not_defined: /* FALLTHRU */                                          \
    case _max_element_type:                                                    \
      break;                                                                   \
    }                                                                          \
  } while (0)

#define AKANTU_BOOST_ELEMENT_SWITCH(macro1, list1)                             \
  AKANTU_BOOST_LIST_SWITCH(macro1, list1, type)

#define AKANTU_BOOST_ELEMENT_SWITCH_NO_DEFAULT(macro1, list1)                  \
  AKANTU_BOOST_LIST_SWITCH_NO_DEFAULT(macro1, list1, type)

#define AKANTU_BOOST_ALL_ELEMENT_SWITCH(macro)                                 \
  AKANTU_BOOST_ELEMENT_SWITCH(macro, AKANTU_ALL_ELEMENT_TYPE)

#define AKANTU_BOOST_ALL_ELEMENT_SWITCH_NO_DEFAULT(macro)                      \
  AKANTU_BOOST_ELEMENT_SWITCH_NO_DEFAULT(macro, AKANTU_ALL_ELEMENT_TYPE)

#define AKANTU_BOOST_LIST_MACRO(r, macro, type) macro(type)

#define AKANTU_BOOST_APPLY_ON_LIST(macro, list)                                \
  BOOST_PP_SEQ_FOR_EACH(AKANTU_BOOST_LIST_MACRO, macro, list)

#define AKANTU_BOOST_ALL_ELEMENT_LIST(macro)                                   \
  AKANTU_BOOST_APPLY_ON_LIST(macro, AKANTU_ALL_ELEMENT_TYPE)

#define AKANTU_GET_ELEMENT_LIST(kind) AKANTU##kind##_ELEMENT_TYPE

#define AKANTU_BOOST_KIND_ELEMENT_SWITCH(macro, kind)                          \
  AKANTU_BOOST_ELEMENT_SWITCH(macro, AKANTU_GET_ELEMENT_LIST(kind))

// BOOST_PP_SEQ_TO_LIST does not exists in Boost < 1.49
#define AKANTU_GENERATE_KIND_LIST(seq)                                         \
  BOOST_PP_TUPLE_TO_LIST(BOOST_PP_SEQ_SIZE(seq), BOOST_PP_SEQ_TO_TUPLE(seq))

#define AKANTU_ELEMENT_KIND_BOOST_LIST                                         \
  AKANTU_GENERATE_KIND_LIST(AKANTU_ELEMENT_KIND)

#define AKANTU_BOOST_ALL_KIND_LIST(macro, list)                                \
  BOOST_PP_LIST_FOR_EACH(AKANTU_BOOST_LIST_MACRO, macro, list)

#define AKANTU_BOOST_ALL_KIND(macro)                                           \
  AKANTU_BOOST_ALL_KIND_LIST(macro, AKANTU_ELEMENT_KIND_BOOST_LIST)

#define AKANTU_BOOST_ALL_KIND_SWITCH(macro)                                    \
  AKANTU_BOOST_LIST_SWITCH(macro, AKANTU_ELEMENT_KIND, kind)


#define AKANTU_BOOST_COHESIVE_ELEMENT_SWITCH(macro)                   \
 AKANTU_BOOST_ELEMENT_SWITCH(macro,                                     \
                              AKANTU_ek_cohesive_ELEMENT_TYPE)

#define AKANTU_BOOST_COHESIVE_ELEMENT_LIST(macro)                     \
  AKANTU_BOOST_APPLY_ON_LIST(macro,                                     \
                             AKANTU_ek_cohesive_ELEMENT_TYPE)

#define AKANTU_BOOST_REGULAR_ELEMENT_SWITCH(macro)                   \
 AKANTU_BOOST_ELEMENT_SWITCH(macro,                                     \
                              AKANTU_ek_regular_ELEMENT_TYPE)

#define AKANTU_BOOST_REGULAR_ELEMENT_LIST(macro)                     \
  AKANTU_BOOST_APPLY_ON_LIST(macro,                                     \
                             AKANTU_ek_regular_ELEMENT_TYPE)

#define AKANTU_BOOST_STRUCTURAL_ELEMENT_SWITCH(macro)                   \
 AKANTU_BOOST_ELEMENT_SWITCH(macro,                                     \
                              AKANTU_ek_structural_ELEMENT_TYPE)

#define AKANTU_BOOST_STRUCTURAL_ELEMENT_LIST(macro)                     \
  AKANTU_BOOST_APPLY_ON_LIST(macro,                                     \
                             AKANTU_ek_structural_ELEMENT_TYPE)


// /// define kept for compatibility reasons (they are most probably not needed
// /// anymore) \todo check if they can be removed
// #define AKANTU_REGULAR_ELEMENT_TYPE	AKANTU_ek_regular_ELEMENT_TYPE
// #define AKANTU_COHESIVE_ELEMENT_TYPE	AKANTU_ek_cohesive_ELEMENT_TYPE
// #define AKANTU_STRUCTURAL_ELEMENT_TYPE  AKANTU_ek_structural_ELEMENT_TYPE
// #define AKANTU_IGFEM_ELEMENT_TYPE       AKANTU_ek_igfem_ELEMENT_TYPE

/* -------------------------------------------------------------------------- */
/* Lists of interests for FEEngineTemplate functions                          */
/* -------------------------------------------------------------------------- */
#define AKANTU_FE_ENGINE_LIST_ASSEMBLE_FIELDS                           \
  AKANTU_GENERATE_KIND_LIST((_ek_regular)				\
                            (_ek_structural))
#define AKANTU_FE_ENGINE_LIST_COMPUTE_SHAPES_DERIVATIVES                \
  AKANTU_GENERATE_KIND_LIST((_ek_regular)				\
                            (_ek_structural))
#define AKANTU_FE_ENGINE_LIST_COMPUTE_SHAPES                            \
  AKANTU_GENERATE_KIND_LIST((_ek_regular)				\
                            (_ek_structural))
#define AKANTU_FE_ENGINE_LIST_INTERPOLATE                               \
  AKANTU_GENERATE_KIND_LIST((_ek_regular))
#define AKANTU_FE_ENGINE_LIST_LAGRANGE_BASE                             \
  AKANTU_GENERATE_KIND_LIST((_ek_cohesive)				\
                            (_ek_regular))
#define AKANTU_FE_ENGINE_LIST_INVERSE_MAP                               \
  AKANTU_GENERATE_KIND_LIST((_ek_cohesive)				\
                            (_ek_regular))
#define AKANTU_FE_ENGINE_LIST_INTERPOLATE_ON_INTEGRATION_POINTS         \
  AKANTU_GENERATE_KIND_LIST((_ek_cohesive)				\
                            (_ek_regular)				\
                            (_ek_structural))
#define AKANTU_FE_ENGINE_LIST_GRADIENT_ON_INTEGRATION_POINTS            \
  AKANTU_GENERATE_KIND_LIST((_ek_cohesive)				\
                            (_ek_regular)				\
                            (_ek_structural))
#define AKANTU_FE_ENGINE_LIST_GET_SHAPES_DERIVATIVES                    \
  AKANTU_GENERATE_KIND_LIST((_ek_cohesive)				\
                            (_ek_regular)				\
                            (_ek_structural))
#define AKANTU_FE_ENGINE_LIST_CONTAINS                                  \
  AKANTU_GENERATE_KIND_LIST((_ek_cohesive)				\
                            (_ek_regular))
#define AKANTU_FE_ENGINE_LIST_COMPUTE_NORMALS_ON_INTEGRATION_POINTS     \
  AKANTU_GENERATE_KIND_LIST((_ek_cohesive)				\
                            (_ek_regular))


} // akantu

#endif /* AKANTU_AKA_ELEMENT_CLASSES_INFO_HH_ */

#include "aka_element_classes_info_inline_impl.hh"

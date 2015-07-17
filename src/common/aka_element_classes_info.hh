/**
 * @file   aka_element_classes_info.hh
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @date   Tue May 19 11:43:07 2015
 *
 * @brief  Declaration of the enums for the element classes
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
/// @boost sequence of element to loop on in global tasks
#define AKANTU_ek_regular_ELEMENT_TYPE		\
  (_point_1)					\
  (_segment_2)					\
  (_segment_3)					\
  (_triangle_3)					\
  (_triangle_6)					\
  (_quadrangle_4)				\
  (_quadrangle_8)				\
  (_tetrahedron_4)				\
  (_tetrahedron_10)				\
  (_pentahedron_6)				\
  (_pentahedron_15)				\
  (_hexahedron_8)				\
  (_hexahedron_20)				\

#if defined(AKANTU_STRUCTURAL_MECHANICS)
#define AKANTU_ek_structural_ELEMENT_TYPE	\
  (_bernoulli_beam_2)				\
  (_bernoulli_beam_3)					\
  (_kirchhoff_shell)
#else
#define AKANTU_ek_structural_ELEMENT_TYPE
#endif

#if defined(AKANTU_COHESIVE_ELEMENT)
#  define AKANTU_ek_cohesive_ELEMENT_TYPE	\
  (_cohesive_2d_4)				\
  (_cohesive_2d_6)				\
  (_cohesive_1d_2)				\
  (_cohesive_3d_6)				\
  (_cohesive_3d_12)
#else
#  define AKANTU_ek_cohesive_ELEMENT_TYPE
#endif

#if defined(AKANTU_IGFEM)
#define AKANTU_ek_igfem_ELEMENT_TYPE		\
  (_igfem_segment_3)				\
  (_igfem_triangle_4)				\
  (_igfem_triangle_5)
#else
#  define AKANTU_ek_igfem_ELEMENT_TYPE
#endif

// This enum cannot be boosted dy to swig
/// @enum ElementType type of elements
enum ElementType {
  _not_defined,
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

#if defined(AKANTU_STRUCTURAL_MECHANICS)
  _bernoulli_beam_2,
  _bernoulli_beam_3,
  _kirchhoff_shell,
#endif

#if defined(AKANTU_COHESIVE_ELEMENT)
  _cohesive_2d_4,
  _cohesive_2d_6,
  _cohesive_1d_2,
  _cohesive_3d_6,
  _cohesive_3d_12,
#endif
#if defined(AKANTU_IGFEM)
  _igfem_segment_3,
  _igfem_triangle_4,
  _igfem_triangle_5,
#endif
  _max_element_type
};



/* -------------------------------------------------------------------------- */
/// @enum GeometricalType type of element potentially contained in a Mesh
enum GeometricalType {
  _gt_point,             ///< point @remark only for some algorithm to be generic like mesh partitioning */
  _gt_segment_2,         ///< 2 nodes segment
  _gt_segment_3,         ///< 3 nodes segment
  _gt_triangle_3,        ///< 3 nodes triangle
  _gt_triangle_6,        ///< 6 nodes triangle
  _gt_quadrangle_4,      ///< 4 nodes quadrangle
  _gt_quadrangle_8,      ///< 8 nodes quadrangle
  _gt_tetrahedron_4,     ///< 4 nodes tetrahedron
  _gt_tetrahedron_10,    ///< 10 nodes tetrahedron
  _gt_hexahedron_8,      ///< 8 nodes hexahedron
  _gt_hexahedron_20,     ///< 20 nodes hexahedron
  _gt_pentahedron_6,     ///< 6 nodes pentahedron
  _gt_pentahedron_15,    ///< 15 nodes pentahedron
#if defined(AKANTU_COHESIVE_ELEMENT)
  _gt_cohesive_2d_4,     ///< 4 nodes 2D cohesive
  _gt_cohesive_2d_6,     ///< 6 nodes 2D cohesive
  _gt_cohesive_1d_2,     ///< 2 nodes 1D cohesive
  _gt_cohesive_3d_6,     ///< 6 nodes 3D cohesive
  _gt_cohesive_3d_12,    ///< 12 nodes 3D cohesive
#endif
#if defined(AKANTU_IGFEM)
  _gt_igfem_segment_3,     ///< 3 nodes igfem segment
  _gt_igfem_triangle_4,     ///< 4 nodes igfem triangle
  _gt_igfem_triangle_5,     ///< 5 nodes igfem triangle
#endif
  _gt_not_defined
};

/* -------------------------------------------------------------------------- */
/// @enum InterpolationType type of elements
enum InterpolationType {
  _itp_lagrange_point_1,           ///< zeroth (!) order lagrangian point (for compatibility purposes)
  _itp_lagrange_segment_2,         ///< first order lagrangian segment
  _itp_lagrange_segment_3,         ///< second order lagrangian segment
  _itp_lagrange_triangle_3,        ///< first order lagrangian triangle
  _itp_lagrange_triangle_6,        ///< second order lagrangian triangle
  _itp_lagrange_quadrangle_4,      ///< first order lagrangian quadrangle
  _itp_serendip_quadrangle_8,      /**< second order serendipian quadrangle
				      @remark used insted of the 9 node
				      lagrangian element */
  _itp_lagrange_tetrahedron_4,     ///< first order lagrangian tetrahedron
  _itp_lagrange_tetrahedron_10,    ///< second order lagrangian tetrahedron
  _itp_lagrange_hexahedron_8,      ///< first order lagrangian hexahedron
  _itp_serendip_hexahedron_20,     ///< second order serendipian hexahedron
  _itp_lagrange_pentahedron_6,     ///< first order lagrangian pentahedron
  _itp_lagrange_pentahedron_15,    ///< second order lagrangian pentahedron
#if defined(AKANTU_STRUCTURAL_MECHANICS)
  _itp_bernoulli_beam,             ///< Bernoulli beam
  _itp_kirchhoff_shell,            ///< Kirchhoff shell
#endif
#if defined(AKANTU_IGFEM)
  _itp_igfem_segment_3,            ///< first oder igfem segment with one enriched node
  _itp_igfem_triangle_4,           ///<first order igfem triangle with one enriched node
  _itp_igfem_triangle_5,           ///<first order igfem triangle with two enriched nodes
#endif
  _itp_not_defined
};

/* -------------------------------------------------------------------------- */
/* Element Kinds                                                              */
/* -------------------------------------------------------------------------- */
#define AKANTU_REGULAR_KIND      (_ek_regular)

#ifdef AKANTU_COHESIVE_ELEMENT
#  define AKANTU_COHESIVE_KIND   (_ek_cohesive)
#else
#  define AKANTU_COHESIVE_KIND
#endif

#ifdef AKANTU_STRUCTURAL_MECHANICS
#  define AKANTU_STRUCTURAL_KIND (_ek_structural)
#else
#  define AKANTU_STRUCTURAL_KIND
#endif

#ifdef AKANTU_IGFEM
#  define AKANTU_IGFEM_KIND      (_ek_igfem)
#else
#  define AKANTU_IGFEM_KIND
#endif

/* -------------------------------------------------------------------------- */
/* Some sub types less probable to change                                     */
/* -------------------------------------------------------------------------- */
/// @enum GeometricalShapeType types of shapes to define the contains
/// function in the element classes
enum GeometricalShapeType {
  _gst_not_defined,
  _gst_point,
  _gst_triangle,
  _gst_square,
  _gst_prism
};

/* -------------------------------------------------------------------------- */
/// @enum GaussIntergrationType classes of types using common
/// description of the gauss point position and weights
enum GaussIntergrationType {
  _git_not_defined,
  _git_point,
  _git_segment,
  _git_triangle,
  _git_tetrahedron,
  _git_pentahedron
};

/* -------------------------------------------------------------------------- */
/// @enum InterpolationKind the family of interpolation types
enum InterpolationKind {
  _itk_not_defined,
  _itk_lagrangian,
  _itk_structural,
  _itk_igfem
};

/* -------------------------------------------------------------------------- */
/* Autmomatic type generation                                                 */
/* -------------------------------------------------------------------------- */
#define AKANTU_ELEMENT_KIND			\
  AKANTU_REGULAR_KIND                           \
  AKANTU_COHESIVE_KIND				\
  AKANTU_STRUCTURAL_KIND			\
  AKANTU_IGFEM_KIND

#define AKANTU_ALL_ELEMENT_TYPE			\
  AKANTU_ek_regular_ELEMENT_TYPE		\
  AKANTU_ek_cohesive_ELEMENT_TYPE		\
  AKANTU_ek_structural_ELEMENT_TYPE		\
  AKANTU_ek_igfem_ELEMENT_TYPE

#define AKANTU_NOT_STRUCTURAL_ELEMENT_TYPE	\
  AKANTU_ek_regular_ELEMENT_TYPE		\
  AKANTU_ek_cohesive_ELEMENT_TYPE		\
  AKANTU_ek_igfem_ELEMENT_TYPE

#ifndef SWIG
enum ElementKind {
  BOOST_PP_SEQ_ENUM(AKANTU_ELEMENT_KIND),
  _ek_not_defined
};
#else
enum ElementKind;
#endif

/* -------------------------------------------------------------------------- */
// BOOST PART: TOUCH ONLY IF YOU KNOW WHAT YOU ARE DOING
#define AKANTU_BOOST_CASE_MACRO(r,macro,_type)				\
  case _type : { macro(_type); break; }

#define AKANTU_BOOST_LIST_SWITCH(macro1, list1, var)			\
  do {									\
    switch(var) {							\
      BOOST_PP_SEQ_FOR_EACH(AKANTU_BOOST_CASE_MACRO, macro1, list1)	\
    default: {								\
      AKANTU_DEBUG_ERROR("Type (" << var << ") not handled by this function"); \
    }									\
    }									\
  } while(0)

#define AKANTU_BOOST_ELEMENT_SWITCH(macro1, list1)			\
  AKANTU_BOOST_LIST_SWITCH(macro1, list1, type)

#define AKANTU_BOOST_ALL_ELEMENT_SWITCH(macro)				\
  AKANTU_BOOST_ELEMENT_SWITCH(macro,					\
			      AKANTU_ALL_ELEMENT_TYPE)

#define AKANTU_BOOST_REGULAR_ELEMENT_SWITCH(macro)			\
  AKANTU_BOOST_ELEMENT_SWITCH(macro,					\
			      AKANTU_ek_regular_ELEMENT_TYPE)

#define AKANTU_BOOST_COHESIVE_ELEMENT_SWITCH(macro)			\
  AKANTU_BOOST_ELEMENT_SWITCH(macro,					\
			      AKANTU_ek_cohesive_ELEMENT_TYPE)

#define AKANTU_BOOST_STRUCTURAL_ELEMENT_SWITCH(macro)			\
  AKANTU_BOOST_ELEMENT_SWITCH(macro,					\
			      AKANTU_ek_structural_ELEMENT_TYPE)

#define AKANTU_BOOST_IGFEM_ELEMENT_SWITCH(macro)			\
  AKANTU_BOOST_ELEMENT_SWITCH(macro,					\
			      AKANTU_ek_igfem_ELEMENT_TYPE)

#define AKANTU_BOOST_LIST_MACRO(r, macro, type)				\
  macro(type)

#define AKANTU_BOOST_APPLY_ON_LIST(macro, list)				\
  BOOST_PP_SEQ_FOR_EACH(AKANTU_BOOST_LIST_MACRO, macro, list)

#define AKANTU_BOOST_ALL_ELEMENT_LIST(macro)				\
  AKANTU_BOOST_APPLY_ON_LIST(macro,					\
			     AKANTU_ALL_ELEMENT_TYPE)

#define AKANTU_BOOST_REGULAR_ELEMENT_LIST(macro)			\
  AKANTU_BOOST_APPLY_ON_LIST(macro,					\
			     AKANTU_ek_regular_ELEMENT_TYPE)

#define AKANTU_BOOST_STRUCTURAL_ELEMENT_LIST(macro)			\
  AKANTU_BOOST_APPLY_ON_LIST(macro,					\
			     AKANTU_ek_structural_ELEMENT_TYPE)

#define AKANTU_BOOST_COHESIVE_ELEMENT_LIST(macro)			\
  AKANTU_BOOST_APPLY_ON_LIST(macro,					\
			     AKANTU_ek_cohesive_ELEMENT_TYPE)

#define AKANTU_BOOST_IGFEM_ELEMENT_LIST(macro)				\
  AKANTU_BOOST_APPLY_ON_LIST(macro,					\
			     AKANTU_ek_igfem_ELEMENT_TYPE)

#define AKANTU_GET_ELEMENT_LIST(kind)					\
  AKANTU##kind##_ELEMENT_TYPE

#define AKANTU_BOOST_KIND_ELEMENT_SWITCH(macro, kind)			\
  AKANTU_BOOST_ELEMENT_SWITCH(macro,					\
			      AKANTU_GET_ELEMENT_LIST(kind))

// BOOST_PP_SEQ_TO_LIST does not exists in Boost < 1.49
#define AKANTU_GENERATE_KIND_LIST(seq)                                  \
  BOOST_PP_TUPLE_TO_LIST(BOOST_PP_SEQ_SIZE(seq),                        \
			 BOOST_PP_SEQ_TO_TUPLE(seq))

#define AKANTU_ELEMENT_KIND_BOOST_LIST AKANTU_GENERATE_KIND_LIST(AKANTU_ELEMENT_KIND)

#define AKANTU_BOOST_ALL_KIND_LIST(macro, list)				\
  BOOST_PP_LIST_FOR_EACH(AKANTU_BOOST_LIST_MACRO, macro, list)

#define AKANTU_BOOST_ALL_KIND(macro)					\
  AKANTU_BOOST_ALL_KIND_LIST(macro, AKANTU_ELEMENT_KIND_BOOST_LIST)

#define AKANTU_BOOST_ALL_KIND_SWITCH(macro)				\
  AKANTU_BOOST_LIST_SWITCH(macro,					\
			   AKANTU_ELEMENT_KIND,				\
			   kind)

/// define kept for compatibility reasons (they are most probably not needed
/// anymore) \todo check if they can be removed
#define AKANTU_REGULAR_ELEMENT_TYPE	AKANTU_ek_regular_ELEMENT_TYPE
#define AKANTU_COHESIVE_ELEMENT_TYPE	AKANTU_ek_cohesive_ELEMENT_TYPE
#define AKANTU_STRUCTURAL_ELEMENT_TYPE  AKANTU_ek_structural_ELEMENT_TYPE
#define AKANTU_IGFEM_ELEMENT_TYPE       AKANTU_ek_igfem_ELEMENT_TYPE

__END_AKANTU__

#include "aka_element_classes_info_inline_impl.cc"

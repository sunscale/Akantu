/**
 * @file   aka_common.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Mon Jun 14 2010
 * @date last modification: Mon Sep 15 2014
 *
 * @brief  common type descriptions for akantu
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
 * @section DESCRIPTION
 *
 * All common things to be included in the projects files
 *
 */

/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_COMMON_HH__
#define __AKANTU_COMMON_HH__

/* -------------------------------------------------------------------------- */
#include <list>
#include <limits>
#include <boost/preprocessor.hpp>

/* -------------------------------------------------------------------------- */
#define __BEGIN_AKANTU__ namespace akantu {
#define __END_AKANTU__ };
/* -------------------------------------------------------------------------- */
#define __BEGIN_AKANTU_DUMPER__ namespace dumper {
#define __END_AKANTU_DUMPER__ }
/* -------------------------------------------------------------------------- */
#if defined(WIN32)
#  define __attribute__(x)
#endif

/* -------------------------------------------------------------------------- */
#include "aka_config.hh"
#include "aka_error.hh"
#include "aka_safe_enum.hh"
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
/* Common types                                                               */
/* -------------------------------------------------------------------------- */

typedef double Real;
typedef unsigned int UInt;
typedef unsigned long long UInt64;
typedef signed int Int;

typedef std::string ID;

static const Real UINT_INIT_VALUE = 0;
#ifdef AKANTU_NDEBUG
  static const Real REAL_INIT_VALUE = 0;
#else
  static const Real REAL_INIT_VALUE = std::numeric_limits<Real>::quiet_NaN();
#endif

/* -------------------------------------------------------------------------- */
/* Memory types                                                               */
/* -------------------------------------------------------------------------- */

typedef UInt MemoryID;

/* -------------------------------------------------------------------------- */
/* Mesh/FEM/Model types                                                       */
/* -------------------------------------------------------------------------- */

typedef std::string Surface;
typedef std::pair<Surface, Surface> SurfacePair;
typedef std::list< SurfacePair > SurfacePairList;

/* -------------------------------------------------------------------------- */

extern const UInt _all_dimensions;

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
  (_hexahedron_8)

#if defined(AKANTU_STRUCTURAL_MECHANICS)
#define AKANTU_ek_structural_ELEMENT_TYPE	\
  (_bernoulli_beam_2)				\
  (_bernoulli_beam_3)  				\
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
   (_igfem_triangle_3)
#else
#define AKANTU_ek_igfem_ELEMENT_TYPE
#endif

#define AKANTU_ALL_ELEMENT_TYPE			\
  AKANTU_ek_regular_ELEMENT_TYPE		\
  AKANTU_ek_cohesive_ELEMENT_TYPE		\
  AKANTU_ek_structural_ELEMENT_TYPE		\
  AKANTU_ek_igfem_ELEMENT_TYPE

#define AKANTU_NOT_STRUCTURAL_ELEMENT_TYPE	\
  AKANTU_ek_regular_ELEMENT_TYPE		\
  AKANTU_ek_cohesive_ELEMENT_TYPE		\
  AKANTU_ek_igfem_ELEMENT_TYPE

/// @enum ElementType type of elements
enum ElementType {
  _not_defined,
  _point_1,
  _segment_2,         ///< first order segment
  _segment_3,         ///< second order segment
  _triangle_3,        ///< first order triangle
  _triangle_6,        ///< second order triangle
  _tetrahedron_4,     ///< first order tetrahedron
  _tetrahedron_10,    ///< second order tetrahedron
  _quadrangle_4,      ///< first order quadrangle
  _quadrangle_8,      ///< second order quadrangle
  _hexahedron_8,      ///< first order hexahedron
  _pentahedron_6,     ///< first order pentahedron
#if defined (AKANTU_STRUCTURAL_MECHANICS)
  _bernoulli_beam_2,  ///< Bernoulli beam 2D
  _bernoulli_beam_3,  ///< Bernoulli beam 3D
  _kirchhoff_shell,   ///< Kirchhoff shell
#endif

#if defined(AKANTU_COHESIVE_ELEMENT)
  _cohesive_2d_4,     ///< first order 2D cohesive
  _cohesive_2d_6,     ///< second order 2D cohesive
  _cohesive_1d_2,     ///< first order 1D cohesive
  _cohesive_3d_6,     ///< first order 3D cohesive
  _cohesive_3d_12,     ///< second order 3D cohesive
#endif
#if defined(AKANTU_IGFEM)
  _igfem_triangle_3,  ///< first order triangle for IGFEM
#endif
  _max_element_type
};

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
  _gt_pentahedron_6,     ///< 6 nodes pentahedron
#if defined(AKANTU_COHESIVE_ELEMENT)
  _gt_cohesive_2d_4,     ///< 4 nodes 2D cohesive
  _gt_cohesive_2d_6,     ///< 6 nodes 2D cohesive
  _gt_cohesive_1d_2,     ///< 2 nodes 1D cohesive
  _gt_cohesive_3d_6,     ///< 6 nodes 3D cohesive
  _gt_cohesive_3d_12,     ///< 12 nodes 3D cohesive
#endif
  _gt_not_defined
};

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
  _itp_lagrange_pentahedron_6,      ///< first order lagrangian pentahedron
#if defined(AKANTU_STRUCTURAL_MECHANICS)
  _itp_bernoulli_beam,             ///< Bernoulli beam
  _itp_kirchhoff_shell,            ///< Kirchhoff shell
#endif

  _itp_not_defined
};

//! standard output stream operator for ElementType
inline std::ostream & operator <<(std::ostream & stream, ElementType type);


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

#define AKANTU_ELEMENT_KIND			\
  AKANTU_REGULAR_KIND                           \
  AKANTU_COHESIVE_KIND				\
  AKANTU_STRUCTURAL_KIND			\
  AKANTU_IGFEM_KIND

enum ElementKind {
  BOOST_PP_SEQ_ENUM(AKANTU_ELEMENT_KIND),
  _ek_not_defined
};

/// small help to use names for directions
enum SpacialDirection {
  _x = 0,
  _y = 1,
  _z = 2
};

/// enum MeshIOType type of mesh reader/writer
enum MeshIOType {
  _miot_auto,        ///< Auto guess of the reader to use based on the extension
  _miot_gmsh,        ///< Gmsh files
  _miot_gmsh_struct, ///< Gsmh reader with reintpretation of elements has structures elements
  _miot_diana,       ///< TNO Diana mesh format
  _miot_abaqus       ///< Abaqus mesh format
};

/// enum AnalysisMethod type of solving method used to solve the equation of motion
enum AnalysisMethod {
  _static,
  _implicit_dynamic,
  _explicit_lumped_mass,
  _explicit_lumped_capacity,
  _explicit_consistent_mass
};

//! enum ContactResolutionMethod types of solving for the contact
enum ContactResolutionMethod {
  _penalty,
  _lagrangian,
  _augmented_lagrangian,
  _nitsche,
  _mortar
};

//! enum ContactImplementationMethod types for different contact implementations
enum ContactImplementationMethod {
  _none,
  _uzawa,
  _generalized_newton
};




/// enum SolveConvergenceMethod different resolution algorithms
enum SolveConvergenceMethod {
  _scm_newton_raphson_tangent,             ///< Newton-Raphson with tangent matrix
  _scm_newton_raphson_tangent_modified,    ///< Newton-Raphson with constant tangent matrix
  _scm_newton_raphson_tangent_not_computed ///< Newton-Raphson with constant tangent matrix (user has to assemble it)
};

/// enum SolveConvergenceCriteria different convergence criteria
enum SolveConvergenceCriteria {
  _scc_residual, ///< Use residual to test the convergence
  _scc_increment, ///< Use increment to test the convergence
  _scc_residual_mass_wgh ///< Use residual weighted by inv. nodal mass to testb
};

/// enum CohesiveMethod type of insertion of cohesive elements
enum CohesiveMethod {
  _intrinsic,
  _extrinsic
};

/// myfunction(double * position, double * stress/force, double * normal, unsigned int material_id)
typedef void (*BoundaryFunction)(double *,double *, double*, unsigned int);

/// @enum BoundaryFunctionType type of function passed for boundary conditions
enum BoundaryFunctionType {
  _bft_stress,
  _bft_traction,
  _bft_traction_local
};

/// @enum SparseMatrixType type of sparse matrix used
enum SparseMatrixType {
  _unsymmetric,
  _symmetric
};

/* -------------------------------------------------------------------------- */
/* Contact                                                                    */
/* -------------------------------------------------------------------------- */
typedef ID ContactID;
typedef ID ContactSearchID;
typedef ID ContactNeighborStructureID;

enum ContactType {
  _ct_not_defined = 0,
  _ct_2d_expli    = 1,
  _ct_3d_expli    = 2,
  _ct_rigid       = 3
};

enum ContactSearchType {
  _cst_not_defined = 0,
  _cst_2d_expli    = 1,
  _cst_expli       = 2
};

enum ContactNeighborStructureType {
  _cnst_not_defined  = 0,
  _cnst_regular_grid = 1,
  _cnst_2d_grid      = 2
};

/* -------------------------------------------------------------------------- */
/* Friction                                                                   */
/* -------------------------------------------------------------------------- */
typedef ID FrictionID;

/* -------------------------------------------------------------------------- */
/* Ghosts handling                                                            */
/* -------------------------------------------------------------------------- */

typedef ID SynchronizerID;

/// @enum CommunicatorType type of communication method to use
enum CommunicatorType {
  _communicator_mpi,
  _communicator_dummy
};

/// @enum SynchronizationTag type of synchronizations
enum SynchronizationTag {
  //--- SolidMechanicsModel tags ---
  _gst_smm_mass,         //< synchronization of the SolidMechanicsModel.mass
  _gst_smm_for_gradu,    //< synchronization of the SolidMechanicsModel.displacement
  _gst_smm_boundary,     //< synchronization of the boundary, forces, velocities and displacement
  _gst_smm_uv,           //< synchronization of the nodal velocities and displacement
  _gst_smm_res,          //< synchronization of the nodal residual
  _gst_smm_init_mat,     //< synchronization of the data to initialize materials
  _gst_smm_stress,       //< synchronization of the stresses to compute the internal forces
  _gst_smmc_facets,      //< synchronization of facet data to setup facet synch
  _gst_smmc_facets_conn, //< synchronization of facet global connectivity
  _gst_smmc_facets_stress, //< synchronization of facets' stress to setup facet synch
  _gst_smmc_damage,      //< synchronization of damage
  //--- CohesiveElementInserter tags ---
  _gst_ce_inserter,      //< synchronization of global nodes id of newly inserted cohesive elements
  //--- GroupManager tags ---
  _gst_gm_clusters,      //< synchronization of clusters
  //--- HeatTransfer tags ---
  _gst_htm_capacity,     //< synchronization of the nodal heat capacity
  _gst_htm_temperature,  //< synchronization of the nodal temperature
  _gst_htm_gradient_temperature,  //< synchronization of the element gradient temperature
  //--- LevelSet tags ---
  /// synchronization of the nodal level set value phi
  _gst_htm_phi,
  /// synchronization of the element gradient phi
  _gst_htm_gradient_phi,
  //--- Material non local ---
  _gst_mnl_for_average,  //< synchronization of data to average in non local material
  _gst_mnl_weight,       //< synchronization of data for the weight computations
  //--- General tags ---
  _gst_test,             //< Test tag
  _gst_material_id,      //< synchronization of the material ids
  _gst_for_dump,         //< everything that needs to be synch before dump
  //--- Contact & Friction ---
  _gst_cf_nodal,         //< synchronization of disp, velo, and current position
  _gst_cf_incr           //< synchronization of increment
};

/// standard output stream operator for SynchronizationTag
inline std::ostream & operator <<(std::ostream & stream, SynchronizationTag type);

/// @enum GhostType type of ghost
enum GhostType {
  _not_ghost,
  _ghost,
  _casper  // not used but a real cute ghost
};

/* -------------------------------------------------------------------------- */
struct GhostType_def {
  typedef GhostType type;
  static const type _begin_ = _not_ghost;
  static const type _end_   = _casper;
};

typedef safe_enum<GhostType_def> ghost_type_t;

/// standard output stream operator for GhostType
inline std::ostream & operator <<(std::ostream & stream, GhostType type);


/// @enum SynchronizerOperation reduce operation that the synchronizer can perform
enum SynchronizerOperation {
  _so_sum,
  _so_min,
  _so_max,
  _so_null
};


/* -------------------------------------------------------------------------- */
/* Global defines                                                             */
/* -------------------------------------------------------------------------- */
#define AKANTU_MIN_ALLOCATION 2000

#define AKANTU_INDENT " "
#define AKANTU_INCLUDE_INLINE_IMPL

/* -------------------------------------------------------------------------- */
// int 2 type construct
template <int d>
struct Int2Type {
  enum { result = d };
};

// type 2 type construct
template <class T>
class Type2Type {
  typedef T OriginalType;
};

/* -------------------------------------------------------------------------- */
template<class T>
struct is_scalar {
  enum{ value = false };
};

#define AKANTU_SPECIFY_IS_SCALAR(type)		\
  template<>					\
  struct is_scalar<type> {			\
    enum { value = true };			\
  }


AKANTU_SPECIFY_IS_SCALAR(Real);
AKANTU_SPECIFY_IS_SCALAR(UInt);
AKANTU_SPECIFY_IS_SCALAR(Int);
AKANTU_SPECIFY_IS_SCALAR(bool);

template < typename T1, typename T2 >
struct is_same {
  enum { value = false };      // is_same represents a bool.
};

template < typename T >
struct is_same<T, T> {
  enum { value = true };
};

/* -------------------------------------------------------------------------- */
#define AKANTU_SET_MACRO(name, variable, type)	\
  inline void set##name (type variable) {	\
    this->variable = variable;			\
  }

#define AKANTU_GET_MACRO(name, variable, type)	\
  inline type get##name () const {		\
    return variable;				\
  }

#define AKANTU_GET_MACRO_NOT_CONST(name, variable, type)	\
  inline type get##name () {					\
    return variable;						\
  }

#define AKANTU_GET_MACRO_BY_SUPPORT_TYPE(name, variable, type,		\
					 support, con)			\
  inline con Array<type> &						\
  get##name (const support & el_type,					\
	     const GhostType & ghost_type = _not_ghost) con {		\
    return variable(el_type, ghost_type);				\
  }

#define AKANTU_GET_MACRO_BY_ELEMENT_TYPE(name, variable, type)		\
  AKANTU_GET_MACRO_BY_SUPPORT_TYPE(name, variable, type, ElementType,)
#define AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(name, variable, type)	\
  AKANTU_GET_MACRO_BY_SUPPORT_TYPE(name, variable, type, ElementType, const)

#define AKANTU_GET_MACRO_BY_GEOMETRIE_TYPE(name, variable, type)	\
  AKANTU_GET_MACRO_BY_SUPPORT_TYPE(name, variable, type, GeometricalType,)
#define AKANTU_GET_MACRO_BY_GEOMETRIE_TYPE_CONST(name, variable, type)	\
  AKANTU_GET_MACRO_BY_SUPPORT_TYPE(name, variable, type, GeometricalType, const)

/* -------------------------------------------------------------------------- */
/// initialize the static part of akantu
void initialize(int & argc, char ** & argv);
/// initialize the static part of akantu and read the global input_file
void initialize(const std::string & input_file, int & argc, char ** & argv);
/* -------------------------------------------------------------------------- */
/// finilize correctly akantu and clean the memory
void finalize ();
/* -------------------------------------------------------------------------- */

/*
 * For intel compiler annoying remark
 */
#if defined(__INTEL_COMPILER)
/// remark #981: operands are evaluated in unspecified order
#pragma warning ( disable : 981 )

/// remark #383: value copied to temporary, reference to temporary used
#pragma warning ( disable : 383 )

#endif //defined(__INTEL_COMPILER)

/* -------------------------------------------------------------------------- */
/* string manipulation                                                        */
/* -------------------------------------------------------------------------- */
inline std::string to_lower(const std::string & str);
/* -------------------------------------------------------------------------- */
inline std::string trim(const std::string & to_trim);
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
/// give a string representation of the a human readable size in bit
template<typename T>
std::string printMemorySize(UInt size);
/* -------------------------------------------------------------------------- */


__END_AKANTU__

#include "aka_fwd.hh"

__BEGIN_AKANTU__

/// get access to the internal argument parser
cppargparse::ArgumentParser & getStaticArgumentParser();

/// get access to the internal input file parser
Parser & getStaticParser();

/// get access to the user part of the internal input file parser
const ParserSection & getUserParser();

__END_AKANTU__

/* -------------------------------------------------------------------------- */
// BOOST PART: TOUCH ONLY IF YOU KNOW WHAT YOU ARE DOING

#define AKANTU_BOOST_CASE_MACRO(r,macro,type)				\
  case type : { macro(type); break; }

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

#define AKANTU_BOOST_LIST_MACRO(r, macro, type)	                        \
  macro(type)

#define AKANTU_BOOST_APPLY_ON_LIST(macro, list)			        \
  BOOST_PP_SEQ_FOR_EACH(AKANTU_BOOST_LIST_MACRO, macro, list)

#define AKANTU_BOOST_ALL_ELEMENT_LIST(macro)				\
  AKANTU_BOOST_APPLY_ON_LIST(macro,					\
			     AKANTU_ALL_ELEMENT_TYPE)

#define AKANTU_BOOST_REGULAR_ELEMENT_LIST(macro)		        \
  AKANTU_BOOST_APPLY_ON_LIST(macro,				        \
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

#define AKANTU_GET_ELEMENT_LIST(kind)	                                \
  AKANTU##kind##_ELEMENT_TYPE

#define AKANTU_BOOST_KIND_ELEMENT_SWITCH(macro, kind)			\
  AKANTU_BOOST_ELEMENT_SWITCH(macro,					\
			      AKANTU_GET_ELEMENT_LIST(kind))

// BOOST_PP_SEQ_TO_LIST does not exists in Boost < 1.49
#define AKANTU_GENERATE_KIND_LIST(seq)                                  \
  BOOST_PP_TUPLE_TO_LIST(BOOST_PP_SEQ_SIZE(seq),                        \
                         BOOST_PP_SEQ_TO_TUPLE(seq))

#define AKANTU_ELEMENT_KIND_BOOST_LIST AKANTU_GENERATE_KIND_LIST(AKANTU_ELEMENT_KIND)

#define AKANTU_BOOST_ALL_KIND_LIST(macro, list)			        \
  BOOST_PP_LIST_FOR_EACH(AKANTU_BOOST_LIST_MACRO, macro, list)

#define AKANTU_BOOST_ALL_KIND(macro)					\
  AKANTU_BOOST_ALL_KIND_LIST(macro, AKANTU_ELEMENT_KIND_BOOST_LIST)

#define AKANTU_BOOST_ALL_KIND_SWITCH(macro)			        \
  AKANTU_BOOST_LIST_SWITCH(macro,				        \
			   AKANTU_ELEMENT_KIND,			        \
			   kind)

/// define kept for compatibility reasons (they are most probably not needed
/// anymore) \todo check if they can be removed
#define AKANTU_REGULAR_ELEMENT_TYPE	AKANTU_ek_regular_ELEMENT_TYPE
#define AKANTU_COHESIVE_ELEMENT_TYPE	AKANTU_ek_cohesive_ELEMENT_TYPE
#define AKANTU_STRUCTURAL_ELEMENT_TYPE  AKANTU_ek_structural_ELEMENT_TYPE
#define AKANTU_IGFEM_ELEMENT_TYPE       AKANTU_ek_igfem_ELEMENT_TYPE


#include "aka_common_inline_impl.cc"

#endif /* __AKANTU_COMMON_HH__ */

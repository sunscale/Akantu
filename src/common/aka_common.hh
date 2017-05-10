/**
 * @file   aka_common.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Mon Jun 14 2010
 * @date last modification: Thu Jan 21 2016
 *
 * @brief  common type descriptions for akantu
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2012, 2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de
 * Lausanne)  Laboratory (LSMS  -  Laboratoire de  Simulation  en Mécanique  des
 * Solides)
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
#include <type_traits>

#if __cplusplus < 201402L
namespace std {
template< bool B, class T = void >
using enable_if_t = typename enable_if<B,T>::type;
}
#endif
/* -------------------------------------------------------------------------- */
#define __BEGIN_AKANTU__ namespace akantu {
#define __END_AKANTU__ }
/* -------------------------------------------------------------------------- */
#define __BEGIN_AKANTU_DUMPER__ namespace dumper {
#define __END_AKANTU_DUMPER__ }
/* -------------------------------------------------------------------------- */
#if defined(WIN32)
#define __attribute__(x)
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
typedef std::string ID;

#ifdef AKANTU_NDEBUG
static const Real REAL_INIT_VALUE = Real(0.);
#else
static const Real REAL_INIT_VALUE = std::numeric_limits<Real>::quiet_NaN();
#endif

/* -------------------------------------------------------------------------- */
/* Memory types                                                               */
/* -------------------------------------------------------------------------- */

typedef UInt MemoryID;

typedef std::string Surface;
typedef std::pair<Surface, Surface> SurfacePair;
typedef std::list<SurfacePair> SurfacePairList;

/* -------------------------------------------------------------------------- */
extern const UInt _all_dimensions;

/* -------------------------------------------------------------------------- */
/* Mesh/FEM/Model types                                                       */
/* -------------------------------------------------------------------------- */
__END_AKANTU__

#include "aka_element_classes_info.hh"

__BEGIN_AKANTU__

/// small help to use names for directions
enum SpacialDirection { _x = 0, _y = 1, _z = 2 };

/// enum MeshIOType type of mesh reader/writer
enum MeshIOType {
  _miot_auto,        ///< Auto guess of the reader to use based on the extension
  _miot_gmsh,        ///< Gmsh files
  _miot_gmsh_struct, ///< Gsmh reader with reintpretation of elements has
                     /// structures elements
  _miot_diana,       ///< TNO Diana mesh format
  _miot_abaqus       ///< Abaqus mesh format
};

/// enum AnalysisMethod type of solving method used to solve the equation of
/// motion
enum AnalysisMethod {
  _static = 0,
  _implicit_dynamic = 1,
  _explicit_lumped_mass = 2,
  _explicit_lumped_capacity = 2,
  _explicit_consistent_mass = 3
};

/// enum DOFSupportType defines which kind of dof that can exists
enum DOFSupportType { _dst_nodal, _dst_generic };

/// Type of non linear resolution available in akantu
enum NonLinearSolverType {
  _nls_linear,                  ///< No non linear convergence loop
  _nls_newton_raphson,          ///< Regular Newton-Raphson
  _nls_newton_raphson_modified, ///< Newton-Raphson with initial tangent
  _nls_lumped,                  ///< Case of lumped mass or equivalent matrix
  _nls_auto                     ///< This will take a default value that make sense in case of
                                ///  model::getNewSolver
};

/// Define the node/dof type
enum NodeType : Int {
  _nt_pure_gost = -3,
  _nt_master = -2,
  _nt_normal = -1
};

/// Type of time stepping solver
enum TimeStepSolverType {
  _tsst_static,         ///< Static solution
  _tsst_dynamic,        ///< Dynamic solver
  _tsst_dynamic_lumped, ///< Dynamic solver with lumped mass
};

/// Type of integration scheme
enum IntegrationSchemeType {
  _ist_pseudo_time,             ///< Pseudo Time
  _ist_forward_euler,           ///< GeneralizedTrapezoidal(0)
  _ist_trapezoidal_rule_1,      ///< GeneralizedTrapezoidal(1/2)
  _ist_backward_euler,          ///< GeneralizedTrapezoidal(1)
  _ist_central_difference,      ///< NewmarkBeta(0, 1/2)
  _ist_fox_goodwin,             ///< NewmarkBeta(1/6, 1/2)
  _ist_trapezoidal_rule_2,      ///< NewmarkBeta(1/2, 1/2)
  _ist_linear_acceleration,     ///< NewmarkBeta(1/3, 1/2)
  _ist_newmark_beta,            ///< generic NewmarkBeta with user defined
                                /// alpha and beta
  _ist_generalized_trapezoidal  ///< generic GeneralizedTrapezoidal with user
                                ///  defined alpha
};

/// enum SolveConvergenceCriteria different convergence criteria
enum SolveConvergenceCriteria {
  _scc_residual,         ///< Use residual to test the convergence
  _scc_solution,         ///< Use solution to test the convergence
  _scc_residual_mass_wgh ///< Use residual weighted by inv. nodal mass to testb
};

/// enum CohesiveMethod type of insertion of cohesive elements
enum CohesiveMethod { _intrinsic, _extrinsic };

/// @enum SparseMatrixType type of sparse matrix used
enum MatrixType { _unsymmetric, _symmetric };

/* -------------------------------------------------------------------------- */
/* Ghosts handling                                                            */
/* -------------------------------------------------------------------------- */

typedef ID SynchronizerID;

/// @enum CommunicatorType type of communication method to use
enum CommunicatorType { _communicator_mpi, _communicator_dummy };

/// @enum SynchronizationTag type of synchronizations
enum SynchronizationTag {
  //--- Generic tags ---
  _gst_whatever,
  _gst_update,
  _gst_size,

  //--- SolidMechanicsModel tags ---
  _gst_smm_mass,      ///< synchronization of the SolidMechanicsModel.mass
  _gst_smm_for_gradu, ///< synchronization of the
                      /// SolidMechanicsModel.displacement
  _gst_smm_boundary,  ///< synchronization of the boundary, forces, velocities
                      /// and displacement
  _gst_smm_uv,  ///< synchronization of the nodal velocities and displacement
  _gst_smm_res, ///< synchronization of the nodal residual
  _gst_smm_init_mat, ///< synchronization of the data to initialize materials
  _gst_smm_stress,  ///< synchronization of the stresses to compute the internal
                    /// forces
  _gst_smmc_facets, ///< synchronization of facet data to setup facet synch
  _gst_smmc_facets_conn,   ///< synchronization of facet global connectivity
  _gst_smmc_facets_stress, ///< synchronization of facets' stress to setup facet
                           /// synch
  _gst_smmc_damage,        ///< synchronization of damage

  // --- GlobalIdsUpdater tags ---
  _gst_giu_global_conn, ///< synchronization of global connectivities

  // --- CohesiveElementInserter tags ---
  _gst_ce_groups, ///< synchronization of cohesive element insertion depending
                  /// on facet groups

  // --- GroupManager tags ---
  _gst_gm_clusters, ///< synchronization of clusters

  // --- HeatTransfer tags ---
  _gst_htm_capacity,             ///< synchronization of the nodal heat capacity
  _gst_htm_temperature,          ///< synchronization of the nodal temperature
  _gst_htm_gradient_temperature, ///< synchronization of the element gradient
                                 /// temperature
  // --- LevelSet tags ---
  _gst_htm_phi,          ///< synchronization of the nodal level set value phi
  _gst_htm_gradient_phi, ///< synchronization of the element gradient phi

  //--- Material non local ---
  _gst_mnl_for_average, ///< synchronization of data to average in non local
                        /// material
  _gst_mnl_weight,      ///< synchronization of data for the weight computations

  // --- NeighborhoodSynchronization tags ---
  _gst_nh_criterion,

  // --- General tags ---
  _gst_test,        ///< Test tag
  _gst_user_1,      ///< tag for user simulations
  _gst_user_2,      ///< tag for user simulations
  _gst_material_id, ///< synchronization of the material ids
  _gst_for_dump,    ///< everything that needs to be synch before dump

  // --- Contact & Friction ---
  _gst_cf_nodal, ///< synchronization of disp, velo, and current position
  _gst_cf_incr,  ///< synchronization of increment

  // --- Solver tags ---
  _gst_solver_solution ///< synchronization of the solution obained with the
                       /// PETSc solver
};

/// standard output stream operator for SynchronizationTag
inline std::ostream & operator<<(std::ostream & stream,
                                 SynchronizationTag type);

/// @enum GhostType type of ghost
enum GhostType {
  _not_ghost,
  _ghost,
  _casper // not used but a real cute ghost
};

/* -------------------------------------------------------------------------- */
struct GhostType_def {
  typedef GhostType type;
  static const type _begin_ = _not_ghost;
  static const type _end_ = _casper;
};

typedef safe_enum<GhostType_def> ghost_type_t;

/// standard output stream operator for GhostType
inline std::ostream & operator<<(std::ostream & stream, GhostType type);

/// @enum SynchronizerOperation reduce operation that the synchronizer can
/// perform
enum SynchronizerOperation {
  _so_sum,
  _so_min,
  _so_max,
  _so_prod,
  _so_land,
  _so_band,
  _so_lor,
  _so_bor,
  _so_lxor,
  _so_bxor,
  _so_min_loc,
  _so_max_loc,
  _so_null
};

/* -------------------------------------------------------------------------- */
/* Global defines                                                             */
/* -------------------------------------------------------------------------- */
#define AKANTU_MIN_ALLOCATION 2000

#define AKANTU_INDENT " "
#define AKANTU_INCLUDE_INLINE_IMPL

/* -------------------------------------------------------------------------- */
template <class T> struct is_scalar {
  enum { value = false };
};

#define AKANTU_SPECIFY_IS_SCALAR(type)                                         \
  template <> struct is_scalar<type> {                                         \
    enum { value = true };                                                     \
  }

AKANTU_SPECIFY_IS_SCALAR(Real);
AKANTU_SPECIFY_IS_SCALAR(UInt);
AKANTU_SPECIFY_IS_SCALAR(Int);
AKANTU_SPECIFY_IS_SCALAR(bool);

template <typename T1, typename T2> struct is_same {
  enum { value = false }; // is_same represents a bool.
};

template <typename T> struct is_same<T, T> {
  enum { value = true };
};

/* -------------------------------------------------------------------------- */
#define AKANTU_SET_MACRO(name, variable, type)                                 \
  inline void set##name(type variable) { this->variable = variable; }

#define AKANTU_GET_MACRO(name, variable, type)                                 \
  inline type get##name() const { return variable; }

#define AKANTU_GET_MACRO_NOT_CONST(name, variable, type)                       \
  inline type get##name() { return variable; }

#define AKANTU_GET_MACRO_BY_SUPPORT_TYPE(name, variable, type, support, con)   \
  inline con Array<type> & get##name(                                          \
      const support & el_type, const GhostType & ghost_type = _not_ghost)      \
      con {                                                                    \
    return variable(el_type, ghost_type);                                      \
  }

#define AKANTU_GET_MACRO_BY_ELEMENT_TYPE(name, variable, type)                 \
  AKANTU_GET_MACRO_BY_SUPPORT_TYPE(name, variable, type, ElementType, )
#define AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(name, variable, type)           \
  AKANTU_GET_MACRO_BY_SUPPORT_TYPE(name, variable, type, ElementType, const)

#define AKANTU_GET_MACRO_BY_GEOMETRIE_TYPE(name, variable, type)               \
  AKANTU_GET_MACRO_BY_SUPPORT_TYPE(name, variable, type, GeometricalType, )
#define AKANTU_GET_MACRO_BY_GEOMETRIE_TYPE_CONST(name, variable, type)         \
  AKANTU_GET_MACRO_BY_SUPPORT_TYPE(name, variable, type, GeometricalType, const)

/* -------------------------------------------------------------------------- */
/// initialize the static part of akantu
void initialize(int & argc, char **& argv);
/// initialize the static part of akantu and read the global input_file
void initialize(const std::string & input_file, int & argc, char **& argv);
/* -------------------------------------------------------------------------- */
/// finilize correctly akantu and clean the memory
void finalize();
/* -------------------------------------------------------------------------- */
/// Read an new input file
void readInputFile(const std::string & input_file);
/* -------------------------------------------------------------------------- */

/*
 * For intel compiler annoying remark
 */
// #if defined(__INTEL_COMPILER)
// /// remark #981: operands are evaluated in unspecified order
// #pragma warning(disable : 981)
// /// remark #383: value copied to temporary, reference to temporary used
// #pragma warning(disable : 383)
// #endif // defined(__INTEL_COMPILER)

/* -------------------------------------------------------------------------- */
/* string manipulation                                                        */
/* -------------------------------------------------------------------------- */
inline std::string to_lower(const std::string & str);
/* -------------------------------------------------------------------------- */
inline std::string trim(const std::string & to_trim);
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
/// give a string representation of the a human readable size in bit
template <typename T> std::string printMemorySize(UInt size);
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

#include "aka_common_inline_impl.cc"

/* -------------------------------------------------------------------------- */

#if defined(AKANTU_UNORDERED_MAP_IS_CXX11)

__BEGIN_AKANTU_UNORDERED_MAP__
#if AKANTU_INTEGER_SIZE == 4
#define AKANTU_HASH_COMBINE_MAGIC_NUMBER 0x9e3779b9
#elif AKANTU_INTEGER_SIZE == 8
#define AKANTU_HASH_COMBINE_MAGIC_NUMBER 0x9e3779b97f4a7c13LL
#endif

/**
 * Hashing function for pairs based on hash_combine from boost The magic number
 * is coming from the golden number @f[\phi = \frac{1 + \sqrt5}{2}@f]
 * @f[\frac{2^32}{\phi} = 0x9e3779b9@f]
 * http://stackoverflow.com/questions/4948780/magic-number-in-boosthash-combine
 * http://burtleburtle.net/bob/hash/doobs.html
 */
template <typename a, typename b> struct hash<std::pair<a, b> > {
public:
  hash() : ah(), bh() {}
  size_t operator()(const std::pair<a, b> & p) const {
    size_t seed = ah(p.first);
    return bh(p.second) + AKANTU_HASH_COMBINE_MAGIC_NUMBER + (seed << 6) +
           (seed >> 2);
  }

private:
  const hash<a> ah;
  const hash<b> bh;
};

__END_AKANTU_UNORDERED_MAP__

#endif

#endif /* __AKANTU_COMMON_HH__ */

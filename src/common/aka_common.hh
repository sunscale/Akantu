/**
 * @file   aka_common.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Mon Jun 14 2010
 * @date last modification: Mon Feb 12 2018
 *
 * @brief  common type descriptions for akantu
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

#ifndef __AKANTU_COMMON_HH__
#define __AKANTU_COMMON_HH__

#include "aka_compatibilty_with_cpp_standard.hh"

/* -------------------------------------------------------------------------- */
#if defined(WIN32)
#define __attribute__(x)
#endif

/* -------------------------------------------------------------------------- */
#include "aka_config.hh"
#include "aka_error.hh"
#include "aka_safe_enum.hh"
/* -------------------------------------------------------------------------- */
#include <boost/preprocessor.hpp>
#include <limits>
#include <list>
#include <memory>
#include <string>
#include <type_traits>
#include <unordered_map>
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
/* Constants                                                                  */
/* -------------------------------------------------------------------------- */
namespace {
  [[gnu::unused]] constexpr UInt _all_dimensions{
      std::numeric_limits<UInt>::max()};
#ifdef AKANTU_NDEBUG
  [[gnu::unused]] constexpr Real REAL_INIT_VALUE{0.};
#else
  [[gnu::unused]] constexpr Real REAL_INIT_VALUE{
      std::numeric_limits<Real>::quiet_NaN()};
#endif
} // namespace

/* -------------------------------------------------------------------------- */
/* Common types                                                               */
/* -------------------------------------------------------------------------- */
using ID = std::string;
using MemoryID = UInt;
} // namespace akantu

/* -------------------------------------------------------------------------- */
#include "aka_enum_macros.hh"
/* -------------------------------------------------------------------------- */
#include "aka_element_classes_info.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {
/* -------------------------------------------------------------------------- */
/* Mesh/FEM/Model types                                                       */
/* -------------------------------------------------------------------------- */

/// small help to use names for directions
enum SpatialDirection { _x = 0, _y = 1, _z = 2 };

/// enum MeshIOType type of mesh reader/writer
enum MeshIOType {
  _miot_auto,        ///< Auto guess of the reader to use based on the extension
  _miot_gmsh,        ///< Gmsh files
  _miot_gmsh_struct, ///< Gsmh reader with reintpretation of elements has
                     /// structures elements
  _miot_diana,       ///< TNO Diana mesh format
  _miot_abaqus       ///< Abaqus mesh format
};

/// enum MeshEventHandlerPriority defines relative order of execution of
/// events
enum EventHandlerPriority {
  _ehp_highest = 0,
  _ehp_mesh = 5,
  _ehp_fe_engine = 9,
  _ehp_synchronizer = 10,
  _ehp_dof_manager = 20,
  _ehp_model = 94,
  _ehp_non_local_manager = 100,
  _ehp_lowest = 100
};

#if !defined(DOXYGEN)
// clang-format off
#define AKANTU_MODEL_TYPES                                              \
  (model)                                                               \
  (solid_mechanics_model)                                               \
  (solid_mechanics_model_cohesive)                                      \
  (heat_transfer_model)                                                 \
  (structural_mechanics_model)                                          \
  (embedded_model)
// clang-format on

/// enum ModelType defines which type of physics is solved
AKANTU_CLASS_ENUM_DECLARE(ModelType, AKANTU_MODEL_TYPES)
AKANTU_CLASS_ENUM_OUTPUT_STREAM(ModelType, AKANTU_MODEL_TYPES)
AKANTU_CLASS_ENUM_INPUT_STREAM(ModelType, AKANTU_MODEL_TYPES)
#else
enum class ModelType {
  model,
  solid_mechanics_model,
  solid_mechanics_model_cohesive,
  heat_transfer_model,
  structural_mechanics_model,
  embedded_model,
};
#endif
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

#if !defined(DOXYGEN)
// clang-format off
#define AKANTU_NON_LINEAR_SOLVER_TYPES                                 \
  (linear)                                                             \
  (newton_raphson)                                                     \
  (newton_raphson_modified)                                            \
  (lumped)                                                             \
  (gmres)                                                              \
  (bfgs)                                                               \
  (cg)                                                                 \
  (auto)
// clang-format on
AKANTU_CLASS_ENUM_DECLARE(NonLinearSolverType, AKANTU_NON_LINEAR_SOLVER_TYPES)
AKANTU_CLASS_ENUM_OUTPUT_STREAM(NonLinearSolverType,
                                AKANTU_NON_LINEAR_SOLVER_TYPES)
AKANTU_CLASS_ENUM_INPUT_STREAM(NonLinearSolverType,
                               AKANTU_NON_LINEAR_SOLVER_TYPES)
#else
/// Type of non linear resolution available in akantu
enum class NonLinearSolverType {
  _linear,                  ///< No non linear convergence loop
  _newton_raphson,          ///< Regular Newton-Raphson
  _newton_raphson_modified, ///< Newton-Raphson with initial tangent
  _lumped,                  ///< Case of lumped mass or equivalent matrix
  _gmres,
  _bfgs,
  _cg,
  _auto ///< This will take a default value that make sense in case of
        ///  model::getNewSolver
};
#endif

#if !defined(DOXYGEN)
// clang-format off
#define AKANTU_TIME_STEP_SOLVER_TYPE                                    \
  (static)                                                             \
  (dynamic)                                                            \
  (dynamic_lumped)                                                     \
  (not_defined)
// clang-format on
AKANTU_CLASS_ENUM_DECLARE(TimeStepSolverType, AKANTU_TIME_STEP_SOLVER_TYPE)
AKANTU_CLASS_ENUM_OUTPUT_STREAM(TimeStepSolverType,
                                AKANTU_TIME_STEP_SOLVER_TYPE)
AKANTU_CLASS_ENUM_INPUT_STREAM(TimeStepSolverType, AKANTU_TIME_STEP_SOLVER_TYPE)
#else
/// Type of time stepping solver
enum class TimeStepSolverType {
  _static,         ///< Static solution
  _dynamic,        ///< Dynamic solver
  _dynamic_lumped, ///< Dynamic solver with lumped mass
  _not_defined,    ///< For not defined cases
};
#endif

#if !defined(DOXYGEN)
// clang-format off
#define AKANTU_INTEGRATION_SCHEME_TYPE                                  \
  (pseudo_time)                                                        \
  (forward_euler)                                                      \
  (trapezoidal_rule_1)                                                 \
  (backward_euler)                                                     \
  (central_difference)                                                 \
  (fox_goodwin)                                                        \
  (trapezoidal_rule_2)                                                 \
  (linear_acceleration)                                                \
  (newmark_beta)                                                       \
  (generalized_trapezoidal)
// clang-format on
AKANTU_CLASS_ENUM_DECLARE(IntegrationSchemeType, AKANTU_INTEGRATION_SCHEME_TYPE)
AKANTU_CLASS_ENUM_OUTPUT_STREAM(IntegrationSchemeType,
                                AKANTU_INTEGRATION_SCHEME_TYPE)
AKANTU_CLASS_ENUM_INPUT_STREAM(IntegrationSchemeType,
                               AKANTU_INTEGRATION_SCHEME_TYPE)
#else
/// Type of integration scheme
enum class IntegrationSchemeType {
  _pseudo_time,            ///< Pseudo Time
  _forward_euler,          ///< GeneralizedTrapezoidal(0)
  _trapezoidal_rule_1,     ///< GeneralizedTrapezoidal(1/2)
  _backward_euler,         ///< GeneralizedTrapezoidal(1)
  _central_difference,     ///< NewmarkBeta(0, 1/2)
  _fox_goodwin,            ///< NewmarkBeta(1/6, 1/2)
  _trapezoidal_rule_2,     ///< NewmarkBeta(1/2, 1/2)
  _linear_acceleration,    ///< NewmarkBeta(1/3, 1/2)
  _newmark_beta,           ///< generic NewmarkBeta with user defined
                           /// alpha and beta
  _generalized_trapezoidal ///< generic GeneralizedTrapezoidal with user
                           ///  defined alpha
};
#endif

#if !defined(DOXYGEN)
// clang-format off
#define AKANTU_SOLVE_CONVERGENCE_CRITERIA       \
  (residual)                                    \
  (solution)                                    \
  (residual_mass_wgh)
// clang-format on
AKANTU_CLASS_ENUM_DECLARE(SolveConvergenceCriteria,
                          AKANTU_SOLVE_CONVERGENCE_CRITERIA)
AKANTU_CLASS_ENUM_OUTPUT_STREAM(SolveConvergenceCriteria,
                                AKANTU_SOLVE_CONVERGENCE_CRITERIA)
AKANTU_CLASS_ENUM_INPUT_STREAM(SolveConvergenceCriteria,
                               AKANTU_SOLVE_CONVERGENCE_CRITERIA)
#else
/// enum SolveConvergenceCriteria different convergence criteria
enum class SolveConvergenceCriteria {
  _residual,         ///< Use residual to test the convergence
  _solution,         ///< Use solution to test the convergence
  _residual_mass_wgh ///< Use residual weighted by inv. nodal mass to
                     ///< testb
};
#endif

/// enum CohesiveMethod type of insertion of cohesive elements
enum CohesiveMethod { _intrinsic, _extrinsic };

/// @enum MatrixType type of sparse matrix used
enum MatrixType { _unsymmetric, _symmetric, _mt_not_defined };

/* -------------------------------------------------------------------------- */
/* Ghosts handling                                                            */
/* -------------------------------------------------------------------------- */
/// @enum CommunicatorType type of communication method to use
enum CommunicatorType { _communicator_mpi, _communicator_dummy };

#if !defined(DOXYGEN)
// clang-format off
#define AKANTU_SYNCHRONIZATION_TAG              \
  (whatever)                                    \
  (update)                                      \
  (ask_nodes)                                   \
  (size)                                        \
  (smm_mass)                                    \
  (smm_for_gradu)                               \
  (smm_boundary)                                \
  (smm_uv)                                      \
  (smm_res)                                     \
  (smm_init_mat)                                \
  (smm_stress)                                  \
  (smmc_facets)                                 \
  (smmc_facets_conn)                            \
  (smmc_facets_stress)                          \
  (smmc_damage)                                 \
  (giu_global_conn)                             \
  (ce_groups)                                   \
  (gm_clusters)                                 \
  (htm_temperature)                             \
  (htm_gradient_temperature)                    \
  (htm_phi)                                     \
  (htm_gradient_phi)                            \
  (mnl_for_average)                             \
  (mnl_weight)                                  \
  (nh_criterion)                                \
  (test)                                        \
  (user_1)                                      \
  (user_2)                                      \
  (material_id)                                 \
  (for_dump)                                    \
  (cf_nodal)                                    \
  (cf_incr)                                     \
  (solver_solution)
// clang-format on
AKANTU_CLASS_ENUM_DECLARE(SynchronizationTag, AKANTU_SYNCHRONIZATION_TAG)
AKANTU_CLASS_ENUM_OUTPUT_STREAM(SynchronizationTag, AKANTU_SYNCHRONIZATION_TAG)
#else
/// @enum SynchronizationTag type of synchronizations
enum class SynchronizationTag {
  //--- Generic tags ---
  _whatever,
  _update,
  _ask_nodes,
  _size,

  //--- SolidMechanicsModel tags ---
  _smm_mass,      ///< synchronization of the SolidMechanicsModel.mass
  _smm_for_gradu, ///< synchronization of the
                  /// SolidMechanicsModel.displacement
  _smm_boundary,  ///< synchronization of the boundary, forces, velocities
                  /// and displacement
  _smm_uv,        ///< synchronization of the nodal velocities and displacement
  _smm_res,       ///< synchronization of the nodal residual
  _smm_init_mat,  ///< synchronization of the data to initialize materials
  _smm_stress,    ///< synchronization of the stresses to compute the
                  ///< internal
                  /// forces
  _smmc_facets,   ///< synchronization of facet data to setup facet synch
  _smmc_facets_conn,   ///< synchronization of facet global connectivity
  _smmc_facets_stress, ///< synchronization of facets' stress to setup
                       ///< facet
                       /// synch
  _smmc_damage,        ///< synchronization of damage

  // --- GlobalIdsUpdater tags ---
  _giu_global_conn, ///< synchronization of global connectivities

  // --- CohesiveElementInserter tags ---
  _ce_groups, ///< synchronization of cohesive element insertion depending
              /// on facet groups

  // --- GroupManager tags ---
  _gm_clusters, ///< synchronization of clusters

  // --- HeatTransfer tags ---
  _htm_temperature,          ///< synchronization of the nodal temperature
  _htm_gradient_temperature, ///< synchronization of the element gradient
                             /// temperature
  // --- LevelSet tags ---
  _htm_phi,          ///< synchronization of the nodal level set value phi
  _htm_gradient_phi, ///< synchronization of the element gradient phi

  //--- Material non local ---
  _mnl_for_average, ///< synchronization of data to average in non local
                    /// material
  _mnl_weight,      ///< synchronization of data for the weight computations

  // --- NeighborhoodSynchronization tags ---
  _nh_criterion,

  // --- General tags ---
  _test,        ///< Test tag
  _user_1,      ///< tag for user simulations
  _user_2,      ///< tag for user simulations
  _material_id, ///< synchronization of the material ids
  _for_dump,    ///< everything that needs to be synch before dump

  // --- Contact & Friction ---
  _cf_nodal, ///< synchronization of disp, velo, and current position
  _cf_incr,  ///< synchronization of increment

  // --- Solver tags ---
  _solver_solution ///< synchronization of the solution obained with the
                   /// PETSc solver
};
#endif

/// @enum GhostType type of ghost
enum GhostType {
  _not_ghost = 0,
  _ghost = 1,
  _casper // not used but a real cute ghost
};

/// Define the flag that can be set to a node
enum class NodeFlag : std::uint8_t {
  _normal = 0x00,
  _distributed = 0x01,
  _master = 0x03,
  _slave = 0x05,
  _pure_ghost = 0x09,
  _shared_mask = 0x0F,
  _periodic = 0x10,
  _periodic_master = 0x30,
  _periodic_slave = 0x50,
  _periodic_mask = 0xF0,
  _local_master_mask = 0xCC, // ~(_master & _periodic_mask)
};

inline NodeFlag operator&(const NodeFlag & a, const NodeFlag & b) {
  using under = std::underlying_type_t<NodeFlag>;
  return NodeFlag(under(a) & under(b));
}

inline NodeFlag operator|(const NodeFlag & a, const NodeFlag & b) {
  using under = std::underlying_type_t<NodeFlag>;
  return NodeFlag(under(a) | under(b));
}

inline NodeFlag & operator|=(NodeFlag & a, const NodeFlag & b) {
  a = a | b;
  return a;
}

inline NodeFlag & operator&=(NodeFlag & a, const NodeFlag & b) {
  a = a & b;
  return a;
}

inline NodeFlag operator~(const NodeFlag & a) {
  using under = std::underlying_type_t<NodeFlag>;
  return NodeFlag(~under(a));
}

std::ostream & operator<<(std::ostream & stream, NodeFlag flag);

} // namespace akantu

AKANTU_ENUM_HASH(GhostType)

namespace akantu {
/* -------------------------------------------------------------------------- */
struct GhostType_def {
  using type = GhostType;
  static const type _begin_ = _not_ghost;
  static const type _end_ = _casper;
};

using ghost_type_t = safe_enum<GhostType_def>;
extern ghost_type_t ghost_types;

/// standard output stream operator for GhostType
inline std::ostream & operator<<(std::ostream & stream, GhostType type);

/* -------------------------------------------------------------------------- */
/* Global defines                                                             */
/* -------------------------------------------------------------------------- */
#define AKANTU_MIN_ALLOCATION 2000

#define AKANTU_INDENT ' '
#define AKANTU_INCLUDE_INLINE_IMPL

/* -------------------------------------------------------------------------- */
#define AKANTU_SET_MACRO(name, variable, type)                                 \
  inline void set##name(type variable) { this->variable = variable; }

#define AKANTU_GET_MACRO(name, variable, type)                                 \
  inline type get##name() const { return variable; }

#define AKANTU_GET_MACRO_NOT_CONST(name, variable, type)                       \
  inline type get##name() { return variable; }

#define AKANTU_GET_MACRO_DEREF_PTR(name, ptr)                                  \
  inline decltype(auto) get##name() const {                                    \
    if (not ptr) {                                                             \
      AKANTU_EXCEPTION("The member " << #ptr << " is not initialized");        \
    }                                                                          \
    return (*ptr);                                                             \
  }

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

/* -------------------------------------------------------------------------- */
/* string manipulation */
/* -------------------------------------------------------------------------- */
inline std::string to_lower(const std::string & str);
/* -------------------------------------------------------------------------- */
inline std::string trim(const std::string & to_trim);
inline std::string trim(const std::string & to_trim, char c);
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
/// give a string representation of the a human readable size in bit
template <typename T> std::string printMemorySize(UInt size);
/* -------------------------------------------------------------------------- */

struct TensorTrait {};
struct TensorProxyTrait {};

} // namespace akantu

/* -------------------------------------------------------------------------- */
/* Type traits                                                                */
/* -------------------------------------------------------------------------- */
namespace aka {

/* ------------------------------------------------------------------------ */
template <typename T> using is_tensor = std::is_base_of<akantu::TensorTrait, T>;
template <typename T>
using is_tensor_proxy = std::is_base_of<akantu::TensorProxyTrait, T>;
/* ------------------------------------------------------------------------ */
template <typename T> using is_scalar = std::is_arithmetic<T>;
/* ------------------------------------------------------------------------ */
template <typename R, typename T,
          std::enable_if_t<std::is_reference<T>::value> * = nullptr>
bool is_of_type(T && t) {
  return (
      dynamic_cast<std::add_pointer_t<
          std::conditional_t<std::is_const<std::remove_reference_t<T>>::value,
                             std::add_const_t<R>, R>>>(&t) != nullptr);
}

/* -------------------------------------------------------------------------- */
template <typename R, typename T> bool is_of_type(std::unique_ptr<T> & t) {
  return (
      dynamic_cast<std::add_pointer_t<
          std::conditional_t<std::is_const<T>::value, std::add_const_t<R>, R>>>(
          t.get()) != nullptr);
}

/* ------------------------------------------------------------------------ */
template <typename R, typename T,
          std::enable_if_t<std::is_reference<T>::value> * = nullptr>
decltype(auto) as_type(T && t) {
  static_assert(
      disjunction<
          std::is_base_of<std::decay_t<T>, std::decay_t<R>>, // down-cast
          std::is_base_of<std::decay_t<R>, std::decay_t<T>>  // up-cast
          >::value,
      "Type T and R are not valid for a as_type conversion");
  return dynamic_cast<std::add_lvalue_reference_t<
      std::conditional_t<std::is_const<std::remove_reference_t<T>>::value,
                         std::add_const_t<R>, R>>>(t);
}

/* -------------------------------------------------------------------------- */
template <typename R, typename T,
          std::enable_if_t<std::is_pointer<T>::value> * = nullptr>
decltype(auto) as_type(T && t) {
  return &as_type<R>(*t);
}

/* -------------------------------------------------------------------------- */
template <typename R, typename T>
decltype(auto) as_type(const std::shared_ptr<T> & t) {
  return std::dynamic_pointer_cast<R>(t);
}

} // namespace aka

#include "aka_common_inline_impl.hh"

#include "aka_fwd.hh"

namespace akantu {
/// get access to the internal argument parser
cppargparse::ArgumentParser & getStaticArgumentParser();

/// get access to the internal input file parser
Parser & getStaticParser();

/// get access to the user part of the internal input file parser
const ParserSection & getUserParser();

#define AKANTU_CURRENT_FUNCTION                                                \
  (std::string(__func__) + "():" + std::to_string(__LINE__))

} // namespace akantu

/* -------------------------------------------------------------------------- */
#if AKANTU_INTEGER_SIZE == 4
#define AKANTU_HASH_COMBINE_MAGIC_NUMBER 0x9e3779b9
#elif AKANTU_INTEGER_SIZE == 8
#define AKANTU_HASH_COMBINE_MAGIC_NUMBER 0x9e3779b97f4a7c13LL
#endif

namespace std {
/**
 * Hashing function for pairs based on hash_combine from boost The magic
 * number is coming from the golden number @f[\phi = \frac{1 + \sqrt5}{2}@f]
 * @f[\frac{2^32}{\phi} = 0x9e3779b9@f]
 * http://stackoverflow.com/questions/4948780/magic-number-in-boosthash-combine
 * http://burtleburtle.net/bob/hash/doobs.html
 */
template <typename a, typename b> struct hash<std::pair<a, b>> {
  hash() = default;
  size_t operator()(const std::pair<a, b> & p) const {
    size_t seed = ah(p.first);
    return bh(p.second) + AKANTU_HASH_COMBINE_MAGIC_NUMBER + (seed << 6) +
           (seed >> 2);
  }

private:
  const hash<a> ah{};
  const hash<b> bh{};
};

} // namespace std


#endif /* __AKANTU_COMMON_HH__ */

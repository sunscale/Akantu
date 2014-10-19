/**
 * @file   structural_mechanics_model.hh
 *
 * @author Sébastien Hartmann <sebastien.hartmann@epfl.ch>
 * @author Damien Spielmann <damien.spielmann@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Fabian Barras <fabian.barras@epfl.ch>
 *
 * @date creation: Fri Jul 15 2011
 * @date last modification: Tue Sep 02 2014
 *
 * @brief  Particular implementation of the structural elements in the StructuralMechanicsModel
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

/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_STRUCTURAL_MECHANICS_MODEL_HH__
#define __AKANTU_STRUCTURAL_MECHANICS_MODEL_HH__

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "model.hh"
#include "integrator_gauss.hh"
#include "shape_linked.hh"
#include "aka_types.hh"
#include "dumpable.hh"
#include "solver.hh"
#include "integration_scheme_2nd_order.hh"

/* -------------------------------------------------------------------------- */
namespace akantu {
  class SparseMatrix;
}

__BEGIN_AKANTU__


struct StructuralMaterial {
  Real E;
  Real A;
  Real I;
  Real Iz;
  Real Iy;
  Real GJ;
  Real rho;
  Real t;
  Real nu;
};

struct StructuralMechanicsModelOptions : public ModelOptions {
  StructuralMechanicsModelOptions(AnalysisMethod analysis_method = _static) :
    analysis_method(analysis_method)  {}
    AnalysisMethod analysis_method;
};

extern const StructuralMechanicsModelOptions default_structural_mechanics_model_options;

class StructuralMechanicsModel : public Model {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  typedef FEEngineTemplate<IntegratorGauss, ShapeLinked, _ek_structural> MyFEEngineType;

  StructuralMechanicsModel(Mesh & mesh,
			   UInt spatial_dimension = _all_dimensions,
			   const ID & id = "structural_mechanics_model",
			   const MemoryID & memory_id = 0);

  virtual ~StructuralMechanicsModel();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  /// initialize fully the model
  void initFull(const ModelOptions & options = default_structural_mechanics_model_options);

  /// initialize the internal vectors
  void initArrays();

  /// initialize the model
  void initModel();

  /// initialize the solver
  void initSolver(SolverOptions & options = _solver_no_options);

  /// initialize the stuff for the implicit solver
  void initImplicit(bool dynamic = false,
		    SolverOptions & solver_options = _solver_no_options);


  /// compute the stresses per elements
  void computeStresses();

  /// assemble the stiffness matrix
  void assembleStiffnessMatrix();

  /// assemble the mass matrix for consistent mass resolutions
  void assembleMass();

  /// implicit time integration predictor
  void implicitPred();

  /// implicit time integration corrector
  void implicitCorr();

  /// update the residual vector
  void updateResidual();

  /// solve the system
  void solve();
 
  bool testConvergenceIncrement(Real tolerance);
  bool testConvergenceIncrement(Real tolerance, Real & error);

  virtual void printself(std::ostream & stream, int indent = 0) const {};

  void computeRotationMatrix(const ElementType & type);

protected:
  UInt getTangentStiffnessVoigtSize(const ElementType & type);

  /// compute Rotation Matrices
  template<const ElementType type>
  void computeRotationMatrix(Array<Real> & rotations) {};

  /// compute A and solve @f[ A\delta u = f_ext - f_int @f]
  template<NewmarkBeta::IntegrationSchemeCorrectorType type>
  void solve(Array<Real> & increment, Real block_val = 1.);

 /* ------------------------------------------------------------------------ */
 /* Mass (structural_mechanics_model_mass.cc)                                     */
 /* ------------------------------------------------------------------------ */

  /// assemble the mass matrix for either _ghost or _not_ghost elements
  void assembleMass(GhostType ghost_type);

  /// computes rho
  void computeRho(Array<Real> & rho,
		  ElementType type,
		  GhostType ghost_type);

  /// finish the computation of residual to solve in increment
  void updateResidualInternal();

  /* ------------------------------------------------------------------------ */

private:
  template<ElementType type>
  inline UInt getTangentStiffnessVoigtSize();

  template <ElementType type>
  void assembleStiffnessMatrix();

  template <ElementType type>
  void assembleMass();

  template<ElementType type>
  void computeStressOnQuad();

  template <ElementType type>
  void computeTangentModuli(Array<Real> & tangent_moduli);

  template <ElementType type>
  void transferBMatrixToSymVoigtBMatrix(Array<Real> & B, bool local = false);

  template <ElementType type>
  void transferNMatrixToSymVoigtNMatrix(Array<Real> & N_matrix);


  /* ------------------------------------------------------------------------ */
  /* Dumpable interface                                                       */
  /* ------------------------------------------------------------------------ */
public:

  virtual dumper::Field * createNodalFieldReal(const std::string & field_name,
					       const std::string & group_name,
					       bool padding_flag);

  virtual dumper::Field * createNodalFieldBool(const std::string & field_name,
					       const std::string & group_name,
					       bool padding_flag);


  virtual dumper::Field * createElementalField(const std::string & field_name, 
					       const std::string & group_name,
					       bool padding_flag,
					       const ElementKind & kind);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  /// set the value of the time step
  void setTimeStep(Real time_step);

  /// return the dimension of the system space
  AKANTU_GET_MACRO(SpatialDimension, spatial_dimension, UInt);

  /// get the StructuralMechanicsModel::displacement vector
  AKANTU_GET_MACRO(Displacement, *displacement_rotation, Array<Real> &);

  /// get the StructuralMechanicsModel::velocity vector
  AKANTU_GET_MACRO(Velocity,        *velocity,               Array<Real> &);

  /// get    the    StructuralMechanicsModel::acceleration    vector,   updated    by
  /// StructuralMechanicsModel::updateAcceleration
  AKANTU_GET_MACRO(Acceleration,    *acceleration,           Array<Real> &);

  /// get the StructuralMechanicsModel::force vector (boundary forces)
  AKANTU_GET_MACRO(Force,        *force_momentum,        Array<Real> &);

 /// get the StructuralMechanicsModel::residual vector, computed by StructuralMechanicsModel::updateResidual
  AKANTU_GET_MACRO(Residual,     *residual,        const Array<Real> &);
  /// get the StructuralMechanicsModel::boundary vector
  AKANTU_GET_MACRO(BlockedDOFs,  *blocked_dofs,    Array<bool> &);

  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(RotationMatrix, rotation_matrix, Real);

  AKANTU_GET_MACRO(StiffnessMatrix, *stiffness_matrix, const SparseMatrix &);

  AKANTU_GET_MACRO(MassMatrix, *mass_matrix, const SparseMatrix &);

  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(Stress, stress, Real);

  AKANTU_GET_MACRO_BY_ELEMENT_TYPE(ElementMaterial, element_material, UInt);

  AKANTU_GET_MACRO_BY_ELEMENT_TYPE(Set_ID, set_ID, UInt);

  void addMaterial(StructuralMaterial & material) { materials.push_back(material); }

  /**
   * @brief set the StructuralMechanicsModel::increment_flag  to on, the activate the
   * update of the StructuralMechanicsModel::increment vector
   */
  void setIncrementFlagOn();

  /* ------------------------------------------------------------------------ */
  /* Boundaries (structural_mechanics_model_boundary.cc)                      */
  /* ------------------------------------------------------------------------ */
public:
  /// Compute Linear load function set in global axis
  template <ElementType type>
  void computeForcesByGlobalTractionArray(const Array<Real> & tractions);

  /// Compute Linear load function set in local axis
  template <ElementType type>
  void computeForcesByLocalTractionArray(const Array<Real> & tractions);

  /// compute force vector from a function(x,y,momentum) that describe stresses
  template <ElementType type>
  void computeForcesFromFunction(BoundaryFunction in_function,
				 BoundaryFunctionType function_type);

 /**
   * solve a step (predictor + convergence loop + corrector) using the
   * the given convergence method (see akantu::SolveConvergenceMethod)
   * and the given convergence criteria (see
   * akantu::SolveConvergenceCriteria)
   **/
  template<SolveConvergenceMethod method, SolveConvergenceCriteria criteria>
  bool solveStep(Real tolerance,
		 UInt max_iteration = 100);

  template<SolveConvergenceMethod method, SolveConvergenceCriteria criteria>
  bool solveStep(Real tolerance,
		 Real & error,
		 UInt max_iteration = 100);

  /// test if the system is converged
  template<SolveConvergenceCriteria criteria>
  bool testConvergence(Real tolerance, Real & error);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// time step
  Real time_step;

  /// conversion coefficient form force/mass to acceleration
  Real f_m2a;

  /// displacements array
  Array<Real> * displacement_rotation;

  /// displacements array at the previous time step (used in finite deformation)
  Array<Real> * previous_displacement;

  /// velocities array
  Array<Real> * velocity;

  /// accelerations array
  Array<Real> * acceleration;

  /// forces array
  Array<Real> * force_momentum;

  /// lumped mass array
  Array<Real> * mass;

  /// stress arraz

  ElementTypeMapArray<Real> stress;

  /// residuals array
  Array<Real> * residual;

  /// boundaries array
  Array<bool> * blocked_dofs;

  /// position of a dof in the K matrix
  Array<Int> * equation_number;

  ElementTypeMapArray<UInt> element_material;

  // Define sets of beams
  ElementTypeMapArray<UInt> set_ID;

  /// local equation_number to global
  unordered_map<UInt, UInt>::type local_eq_num_to_global;

  /// stiffness matrix
  SparseMatrix * stiffness_matrix;

  /// mass matrix
  SparseMatrix * mass_matrix;

  /// velocity damping matrix
  SparseMatrix * velocity_damping_matrix;

  /// jacobian matrix
  SparseMatrix * jacobian_matrix;

  /// increment of displacement
  Array<Real> * increment;

  /// solver for implicit
  Solver * solver;

  /// number of degre of freedom
  UInt nb_degree_of_freedom;

  // Rotation matrix
  ElementTypeMapArray<Real> rotation_matrix;

  /// analysis method check the list in akantu::AnalysisMethod
  AnalysisMethod method;

  /// flag defining if the increment must be computed or not
  bool increment_flag;

  /// integration scheme of second order used
  IntegrationScheme2ndOrder * integrator;

  /* -------------------------------------------------------------------------- */
  std::vector<StructuralMaterial> materials;
};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "structural_mechanics_model_inline_impl.cc"

/// standard output stream operator
inline std::ostream & operator <<(std::ostream & stream, const StructuralMechanicsModel & _this)
{
  _this.printself(stream);
  return stream;
}


__END_AKANTU__

#endif /* __AKANTU_STRUCTURAL_MECHANICS_MODEL_HH__ */

/**
 * @file   structural_mechanics_model.hh
 *
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 * @author Sébastien Hartmann <sebastien.hartmann@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Damien Spielmann <damien.spielmann@epfl.ch>
 *
 * @date creation: Fri Jul 15 2011
 * @date last modification: Tue Feb 20 2018
 *
 * @brief  Particular implementation of the structural elements in the
 * StructuralMechanicsModel
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
#include "aka_named_argument.hh"
#include "boundary_condition.hh"
#include "model.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_STRUCTURAL_MECHANICS_MODEL_HH_
#define AKANTU_STRUCTURAL_MECHANICS_MODEL_HH_

/* -------------------------------------------------------------------------- */
namespace akantu {
class Material;
class MaterialSelector;
class DumperIOHelper;
class NonLocalManager;
template <ElementKind kind, class IntegrationOrderFunctor>
class IntegratorGauss;
template <ElementKind kind> class ShapeStructural;
} // namespace akantu

namespace akantu {

struct StructuralMaterial {
  Real E{0};
  Real A{1};
  Real I{0};
  Real Iz{0};
  Real Iy{0};
  Real GJ{0};
  Real rho{0};
  Real t{0};
  Real nu{0};
};

class StructuralMechanicsModel : public Model {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  using MyFEEngineType =
      FEEngineTemplate<IntegratorGauss, ShapeStructural, _ek_structural>;

  StructuralMechanicsModel(Mesh & mesh, UInt dim = _all_dimensions,
                           const ID & id = "structural_mechanics_model");

  ~StructuralMechanicsModel() override;

  /// Init full model
  void initFullImpl(const ModelOptions & options) override;

  /// Init boundary FEEngine
  void initFEEngineBoundary() override;

  /* ------------------------------------------------------------------------ */
  /* Virtual methods from SolverCallback                                      */
  /* ------------------------------------------------------------------------ */
  /// get the type of matrix needed
  MatrixType getMatrixType(const ID & matrix_id) override;

  /// callback to assemble a Matrix
  void assembleMatrix(const ID & matrix_id) override;

  /// callback to assemble a lumped Matrix
  void assembleLumpedMatrix(const ID & matrix_id) override;

  /// callback to assemble the residual (rhs)
  void assembleResidual() override;

  void assembleResidual(const ID & residual_part) override;

  bool canSplitResidual() override { return false; }

  /// compute kinetic energy
  Real getKineticEnergy();

  /// compute potential energy
  Real getPotentialEnergy();

  /// compute the specified energy
  Real getEnergy(const ID & energy);

  /* ------------------------------------------------------------------------ */
  /* Virtual methods from Model                                               */
  /* ------------------------------------------------------------------------ */
protected:
  /// get some default values for derived classes
  std::tuple<ID, TimeStepSolverType>
  getDefaultSolverID(const AnalysisMethod & method) override;

  ModelSolverOptions
  getDefaultSolverOptions(const TimeStepSolverType & type) const override;

  static UInt getNbDegreeOfFreedom(ElementType type);

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
  void initSolver(TimeStepSolverType time_step_solver_type,
                  NonLinearSolverType non_linear_solver_type) override;

  /// initialize the model
  void initModel() override;

  /// compute the stresses per elements
  void computeStresses();

  /// compute the nodal forces
  void assembleInternalForce();

  /// compute the nodal forces for an element type
  void assembleInternalForce(ElementType type, GhostType gt);

  /// assemble the stiffness matrix
  void assembleStiffnessMatrix();

  /// assemble the mass matrix for consistent mass resolutions
  void assembleMassMatrix();

protected:
  /// assemble the mass matrix for either _ghost or _not_ghost elements
  void assembleMassMatrix(GhostType ghost_type);

  /// computes rho
  void computeRho(Array<Real> & rho, ElementType type, GhostType ghost_type);

  /// finish the computation of residual to solve in increment
  void updateResidualInternal();

  /* ------------------------------------------------------------------------ */
private:
  template <ElementType type> void assembleStiffnessMatrix();
  template <ElementType type> void computeStressOnQuad();
  template <ElementType type>
  void computeTangentModuli(Array<Real> & tangent_moduli);

  /* ------------------------------------------------------------------------ */
  /* Dumpable interface                                                       */
  /* ------------------------------------------------------------------------ */
public:
  std::shared_ptr<dumpers::Field>
  createNodalFieldReal(const std::string & field_name,
                       const std::string & group_name,
                       bool padding_flag) override;

  std::shared_ptr<dumpers::Field>
  createNodalFieldBool(const std::string & field_name,
                       const std::string & group_name,
                       bool padding_flag) override;

  std::shared_ptr<dumpers::Field>
  createElementalField(const std::string & field_name,
                       const std::string & group_name, bool padding_flag,
                       UInt spatial_dimension, ElementKind kind) override;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// set the value of the time step
  void setTimeStep(Real time_step, const ID & solver_id = "") override;

  /// return the dimension of the system space
  AKANTU_GET_MACRO(SpatialDimension, spatial_dimension, UInt);

  /// get the StructuralMechanicsModel::displacement vector
  AKANTU_GET_MACRO(Displacement, *displacement_rotation, Array<Real> &);

  /// get the StructuralMechanicsModel::velocity vector
  AKANTU_GET_MACRO(Velocity, *velocity, Array<Real> &);

  /// get    the    StructuralMechanicsModel::acceleration    vector,   updated
  /// by
  /// StructuralMechanicsModel::updateAcceleration
  AKANTU_GET_MACRO(Acceleration, *acceleration, Array<Real> &);

  /// get the StructuralMechanicsModel::external_force vector
  AKANTU_GET_MACRO(ExternalForce, *external_force, Array<Real> &);

  /// get the StructuralMechanicsModel::internal_force vector (boundary forces)
  AKANTU_GET_MACRO(InternalForce, *internal_force, Array<Real> &);

  /// get the StructuralMechanicsModel::boundary vector
  AKANTU_GET_MACRO(BlockedDOFs, *blocked_dofs, Array<bool> &);

  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(RotationMatrix, rotation_matrix, Real);

  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(Stress, stress, Real);

  AKANTU_GET_MACRO_BY_ELEMENT_TYPE(ElementMaterial, element_material, UInt);

  AKANTU_GET_MACRO_BY_ELEMENT_TYPE(Set_ID, set_ID, UInt);

  /**
   * \brief	This function adds the `StructuralMaterial` material to the list of
   * materials managed by *this.
   *
   * It is important that this function might invalidate all references to
   * structural materials, that were previously optained by `getMaterial()`.
   *
   * \param  material The new material.
   *
   * \return The ID of the material that was added.
   *
   * \note The return type is is new.
   */
  UInt addMaterial(StructuralMaterial & material, const ID & name = "");

  const StructuralMaterial &
  getMaterialByElement(const Element & element) const;

  /**
   * \brief	Returns the ith material of *this.
   * \param  i		The ith material
   */
  const StructuralMaterial & getMaterial(UInt material_index) const;

  const StructuralMaterial & getMaterial(const ID & name) const;

  /**
   * \brief	Returns the number of the different materials inside *this.
   */
  UInt getNbMaterials() const { return materials.size(); }

  /* ------------------------------------------------------------------------ */
  /* Boundaries (structural_mechanics_model_boundary.cc)                      */
  /* ------------------------------------------------------------------------ */
public:
  /// Compute Linear load function set in global axis
  void computeForcesByGlobalTractionArray(const Array<Real> & traction_global,
                                          ElementType type);

  /// Compute Linear load function set in local axis
  void computeForcesByLocalTractionArray(const Array<Real> & tractions,
                                         ElementType type);

  /// compute force vector from a function(x,y,momentum) that describe stresses
  // template <ElementType type>
  // void computeForcesFromFunction(BoundaryFunction in_function,
  //                                BoundaryFunctionType function_type);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// time step
  Real time_step;

  /// conversion coefficient form force/mass to acceleration
  Real f_m2a;

  /// displacements array
  std::unique_ptr<Array<Real>> displacement_rotation;

  /// velocities array
  std::unique_ptr<Array<Real>> velocity;

  /// accelerations array
  std::unique_ptr<Array<Real>> acceleration;

  /// forces array
  std::unique_ptr<Array<Real>> internal_force;

  /// forces array
  std::unique_ptr<Array<Real>> external_force;

  /// lumped mass array
  std::unique_ptr<Array<Real>> mass;

  /// boundaries array
  std::unique_ptr<Array<bool>> blocked_dofs;

  /// stress array
  ElementTypeMapArray<Real> stress;

  ElementTypeMapArray<UInt> element_material;

  // Define sets of beams
  ElementTypeMapArray<UInt> set_ID;

  /// number of degre of freedom
  UInt nb_degree_of_freedom;

  // Rotation matrix
  ElementTypeMapArray<Real> rotation_matrix;

  // /// analysis method check the list in akantu::AnalysisMethod
  // AnalysisMethod method;

  /// flag defining if the increment must be computed or not
  bool increment_flag;

  bool need_to_reassemble_mass{true};
  bool need_to_reassemble_stiffness{true};

  /* ------------------------------------------------------------------------ */
  std::vector<StructuralMaterial> materials;
  std::map<std::string, UInt> materials_names_to_id;
};

} // namespace akantu

#include "structural_mechanics_model_inline_impl.hh"

#endif /* AKANTU_STRUCTURAL_MECHANICS_MODEL_HH_ */

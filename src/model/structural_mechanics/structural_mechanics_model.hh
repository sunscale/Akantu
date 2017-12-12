/**
 * @file   structural_mechanics_model.hh
 *
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @author Sébastien Hartmann <sebastien.hartmann@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Damien Spielmann <damien.spielmann@epfl.ch>
 *
 * @date creation: Fri Jul 15 2011
 * @date last modification: Thu Jan 21 2016
 *
 * @brief  Particular implementation of the structural elements in the
 * StructuralMechanicsModel
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
 */
/* -------------------------------------------------------------------------- */
#include "aka_named_argument.hh"
#include "model.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_STRUCTURAL_MECHANICS_MODEL_HH__
#define __AKANTU_STRUCTURAL_MECHANICS_MODEL_HH__

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

  StructuralMechanicsModel(Mesh & mesh,
                           UInt spatial_dimension = _all_dimensions,
                           const ID & id = "structural_mechanics_model",
                           const MemoryID & memory_id = 0);

  virtual ~StructuralMechanicsModel();

  /* ------------------------------------------------------------------------ */
  /* Virtual methods from SolverCallback */
  /* ------------------------------------------------------------------------ */
  /// get the type of matrix needed
  MatrixType getMatrixType(const ID &) override;

  /// callback to assemble a Matrix
  void assembleMatrix(const ID &) override;

  /// callback to assemble a lumped Matrix
  void assembleLumpedMatrix(const ID &) override;

  /// callback to assemble the residual (rhs)
  void assembleResidual() override;

  /* ------------------------------------------------------------------------ */
  /* Virtual methods from MeshEventHandler */
  /* ------------------------------------------------------------------------ */

  /// function to implement to react on  akantu::NewNodesEvent
  void onNodesAdded(const Array<UInt> & nodes_list,
                    const NewNodesEvent & event) override;
  /// function to implement to react on  akantu::RemovedNodesEvent
  void onNodesRemoved(const Array<UInt> & nodes_list,
                      const Array<UInt> & new_numbering,
                      const RemovedNodesEvent & event) override;
  /// function to implement to react on  akantu::NewElementsEvent
  void onElementsAdded(const Array<Element> & elements_list,
                       const NewElementsEvent & event) override;
  /// function to implement to react on  akantu::RemovedElementsEvent
  void onElementsRemoved(const Array<Element> & elements_list,
                         const ElementTypeMapArray<UInt> & new_numbering,
                         const RemovedElementsEvent & event) override;
  /// function to implement to react on  akantu::ChangedElementsEvent
  void onElementsChanged(const Array<Element> & old_elements_list,
                         const Array<Element> & new_elements_list,
                         const ElementTypeMapArray<UInt> & new_numbering,
                         const ChangedElementsEvent & event) override;

  /* ------------------------------------------------------------------------ */
  /* Virtual methods from Model */
  /* ------------------------------------------------------------------------ */
  /// get some default values for derived classes
  std::tuple<ID, TimeStepSolverType>
  getDefaultSolverID(const AnalysisMethod & method) override;

  ModelSolverOptions
  getDefaultSolverOptions(const TimeStepSolverType & type) const;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  void initSolver(TimeStepSolverType, NonLinearSolverType) override;

  /// initialize the model
  void initModel() override;

  /// compute the stresses per elements
  void computeStresses();

  /// assemble the stiffness matrix
  void assembleStiffnessMatrix();

  /// assemble the mass matrix for consistent mass resolutions
  void assembleMass();

  /// update the residual vector
  void updateResidual();

  void computeRotationMatrix(const ElementType & type);

protected:
  /// compute Rotation Matrices
  template <const ElementType type>
  void computeRotationMatrix(__attribute__((unused)) Array<Real> & rotations) {}

  /* ------------------------------------------------------------------------ */
  /* Mass (structural_mechanics_model_mass.cc) */
  /* ------------------------------------------------------------------------ */

  /// assemble the mass matrix for either _ghost or _not_ghost elements
  void assembleMass(GhostType ghost_type);

  /// computes rho
  void computeRho(Array<Real> & rho, ElementType type, GhostType ghost_type);

  /// finish the computation of residual to solve in increment
  void updateResidualInternal();

  /* ------------------------------------------------------------------------ */
private:
  template <ElementType type> void assembleStiffnessMatrix();
  template <ElementType type> void assembleMass();
  template <ElementType type> void computeStressOnQuad();
  template <ElementType type>
  void computeTangentModuli(Array<Real> & tangent_moduli);

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

  virtual dumper::Field *
  createElementalField(const std::string & field_name,
                       const std::string & group_name, bool padding_flag,
                       const ElementKind & kind,
                       const std::string & fe_engine_id = "");

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// set the value of the time step
  // void setTimeStep(Real time_step, const ID & solver_id = "") override;

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
  AKANTU_GET_MACRO(ExternalForce, *external_force_momentum, Array<Real> &);

  /// get the StructuralMechanicsModel::internal_force vector (boundary forces)
  AKANTU_GET_MACRO(InternalForce, *internal_force_momentum, Array<Real> &);

  /// get the StructuralMechanicsModel::boundary vector
  AKANTU_GET_MACRO(BlockedDOFs, *blocked_dofs, Array<bool> &);

  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(RotationMatrix, rotation_matrix, Real);

  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(Stress, stress, Real);

  AKANTU_GET_MACRO_BY_ELEMENT_TYPE(ElementMaterial, element_material, UInt);

  AKANTU_GET_MACRO_BY_ELEMENT_TYPE(Set_ID, set_ID, UInt);

  void addMaterial(StructuralMaterial & material) {
    materials.push_back(material);
  }

  /**
   * @brief set the StructuralMechanicsModel::increment_flag  to on, the
   * activate the
   * update of the StructuralMechanicsModel::increment vector
   */
  void setIncrementFlagOn();

  const StructuralMaterial & getMaterial(const Element & element) const {
    return materials[element_material(element)];
  }

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
  Array<Real> * displacement_rotation{nullptr};

  /// velocities array
  Array<Real> * velocity{nullptr};

  /// accelerations array
  Array<Real> * acceleration{nullptr};

  /// forces array
  Array<Real> * internal_force_momentum{nullptr};

  /// forces array
  Array<Real> * external_force_momentum{nullptr};

  /// lumped mass array
  Array<Real> * mass{nullptr};

  /// boundaries array
  Array<bool> * blocked_dofs{nullptr};

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

  /* ------------------------------------------------------------------------ */
  std::vector<StructuralMaterial> materials;
};

} // namespace akantu

#endif /* __AKANTU_STRUCTURAL_MECHANICS_MODEL_HH__ */

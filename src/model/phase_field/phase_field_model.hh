/**
 * @file   phase_field_model.hh
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 * 
 * @date creation: Sun Jul 30 2018
 * @date last modification: Mon Feb 05 2018
 *
 * @brief  Model of Phase Field
 *
 * @section LICENSE
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
#include "data_accessor.hh"
#include "fe_engine.hh"
#include "model.hh"
/* -------------------------------------------------------------------------- */
#include <array>
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_PHASE_FIELD_MODEL_HH__
#define __AKANTU_PHASE_FIELD_MODEL_HH__

namespace akantu {
template <ElementKind kind, class IntegrationOrderFuntor>
class IntegratorGauss;
template <ElementKind kind> class Shapelagrange;  
} // namespace akantu

namespace akantu {

class PhaseFieldModel : public Model,
			public DataAccessor<Element>,
			public DataAccessor<UInt> {

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  using FEEngineType = FEEngineTemplate<IntegratorGauss, Shapelagrange>;

  PhaseFieldModel(Mesh & mesh, UInt spatial_dimension = _all_dimensions,
		  const ID & id = "phase_field_model",
		  const MemoryID & memory_id = 0);

  virtual ~PhaseFieldModel();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
protected:
  /// generic function to initialize everything ready for explicit dynamics
  void initFullImpl(const ModelOptions & options) override;

  /// read one material file to instantiate all the materials
  void readMaterials();

  /// allocate all vectors
  void initSolver(TimeStepSolverType, NonLinearSolverType) override;

  /// initialize the model
  void initModel() override;

  void predictor() override;

  /// compute the heat flux
  void assembleResidual() override;

  /// get the type of matrix needed
  MatrixType getMatrixType(const ID &) override;

  /// callback to assemble a Matrix
  void assembleMatrix(const ID &) override;

  std::tuple<ID, TimeStepSolverType>
  getDefaultSolverID(const AnalysisMethod & method) override;

  ModelSolverOptions
  getDefaultSolverOptions(const TimeStepSolverType & type) const;
  
  /// compute the internal forces 
  void assembleInternalForces();

  
  /* ------------------------------------------------------------------------ */
  /* Methods for static                                                       */
  /* ------------------------------------------------------------------------ */
public:
  /// assemble damage matrix
  void assembleDamageMatrix();

  /// assemble damage gradient matrix
  void assembleDamageGradMatrix();
  
  /// compute the damage on quadrature points
  void computeDamage(Array<Real> & damage, ElementType type, GhostType ghost_type);

private:
  /// assemble the damage matrix (w/ ghost type)
  template<UInt dim>
  void assembleDamageMatrix(const GhostType & ghost_type);

   /// assemble the damage grad matrix (w/ ghost type)
  template<UInt dim>
  void assembleDamageGradMatrix(const GhostType & ghost_type);

  /// compute vector strain history field for each quadrature point
  void computeStrainHistoryOnQuadPoints(const GhostType & ghost_type);

  /// compute the fracture energy
  Real computeFractureEnergyByNode();

  
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  ///number of iterations
  UInt n_iter;

  /// damage array
  Array<Real> * damage{nullptr};

  /// damage field on quadrature points
  ElementTypeMapArray<Real> damage_on_qpoints;

  /// critical local fracture energy density on quadrature points
  ElementTypeMapArray<Real> gc_on_qpoints;

  /// vector \phi plus on quadrature points
  ElementTypeMapArray<Real> strain_history_on_qpoints;

  /// displacement field on quadrature points
  ElementTypeMapArray<Real> displacement_on_qpoints;

  /// boundary vector
  Array<bool> * blocked_dofs{nullptr};

  /// external force vector
  Array<Real> * external_force{nullptr};

  /// residuals array
  Array<Real> * internal_force{nullptr};

  /// lengthscale parameter
  Real l_0;

  /// critical energy release rate
  Real g_c;

  /// Lame's parameter
  Real lame_lambda;

  /// Lame's paramter
  Real lame_mu;

};

}

#endif

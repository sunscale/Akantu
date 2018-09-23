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
  void computeDamage(ElementType type, GhostType ghost_type);

  /// coupling parameters damage and strains from solid mechanics model
  void setCouplingParameters(ElementTypeMapArray<Real> & strain_on_qpoints,
			     Array<Real> & damage);
  
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
  /* Data Accessor inherited members                                          */
  /* ------------------------------------------------------------------------ */
public:
  UInt getNbData(const Array<Element> & elements,
		 const SynchronizationTag & tag) const override;

  void packData(CommunicationBuffer & buffer, const Array<Element> & elements,
		const SynchronizationTag & tag) const override;

  void unpackData(CommunicationBuffer & buffer, const Array<Element> & elements,
		  const SynchronizationTag & tag) override;

  UInt getNbData(const Array<UInt> & indexes,
		 const SynchronizationTag & tag) const override;

  void packData(CommunicationBuffer & buffer, const Array<UInt> & dofs,
		const SynchronizationTag & tag) const override;

  void unpackData(CommunicationBuffer & buffer, const Array<UInt> & dofs,
		  const SynchronizationTag & tag) override;

  /* ------------------------------------------------------------------------ */
  /* Dumpable Interface                                                       */
  /* ------------------------------------------------------------------------ */
public:
  dumper::Field * createNodalFieldReal(const std::string & field_name,
				       const std::string & group_name,
				       bool padding_flag) override;

  dumper::Field * createNodalFieldBool(const std::string & field_name,
				       const std::string & group_name,
				       bool padding_flag) override;
  
  dumper::Field * createElementalField(const std::string & field_name,
				       const std::string & group_name,
				       bool padding_flag,
				       const UInt & spatial_dimension,
				       const ElementKind & kind) override;

  virtual void dump(const std::string & dumper_name);

  virtual void dump(const std::string & dumper_name, UInt step);

  virtual void dump(const std::string & dumper_name, Real time, UInt step);

  void dump() override;

  virtual void dump(UInt step);

  virtual void dump(Real time, UInt step);
  
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  ///number of iterations
  UInt n_iter;

  /// damage array
  Array<Real> * damage{nullptr};

  /// the speed of change in damage
  ElementTypeMapArray<Real> damage_gradient;

  /// damage field on quadrature points
  ElementTypeMapArray<Real> damage_on_qpoints;

  /// critical local fracture energy density on quadrature points
  ElementTypeMapArray<Real> gc_on_qpoints;

  /// critical local damage energy density on quadrature points
  ElementTypeMapArray<Real> damage_energy_density_on_qpoints;

  /// critical local damage energy on quadrature points
  ElementTypeMapArray<Real> damage_energy_on_qpoints;

  /// vector \phi plus on quadrature points
  ElementTypeMapArray<Real> strain_history_on_qpoints;

  /// strain on quadrature points
  ElementTypeMapArray<Real> strain_on_qpoints;

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

} // akantu

/* ------------------------------------------------------------------------ */
/* inline functions                                                         */
/* ------------------------------------------------------------------------ */
#include "phase_field_model_inline_impl.cc"

#endif

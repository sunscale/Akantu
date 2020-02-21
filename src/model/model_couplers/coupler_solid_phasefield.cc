/**
 * @file   coupler_solid_phasefield.cc
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Fri Sep 28 2018
 * @date last modification: Thu Jun 20 2019
 *
 * @brief  class for coupling of solid mechancis and phasefield model
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
#include "coupler_solid_phasefield.hh"
#include "dumpable_inline_impl.hh"
#include "integrator_gauss.hh"
#include "shape_lagrange.hh"

#ifdef AKANTU_USE_IOHELPER
#include "dumper_iohelper_paraview.hh"
#endif
/* -------------------------------------------------------------------------- */

namespace akantu {

CouplerSolidPhaseField::CouplerSolidPhaseField(Mesh & mesh, UInt dim,
                                               const ID & id,
                                               const MemoryID & memory_id,
                                               const ModelType model_type)
    : Model(mesh, model_type, dim, id, memory_id) {

  AKANTU_DEBUG_IN();

  this->registerFEEngineObject<MyFEEngineType>("CouplerSolidPhaseField", mesh,
                                               Model::spatial_dimension);

#if defined(AKANTU_USE_IOHELPER)
  this->mesh.registerDumper<DumperParaview>("coupler_solid_phasefield", id,
                                            true);
  this->mesh.addDumpMeshToDumper("coupler_solid_phasefield", mesh,
                                 Model::spatial_dimension, _not_ghost,
                                 _ek_regular);
#endif

  this->registerDataAccessor(*this);

  solid = new SolidMechanicsModel(mesh, Model::spatial_dimension,
                                  "solid_mechanics_model", 0);
  phasefield = new PhaseFieldModel(mesh, Model::spatial_dimension,
                                   "phase_field_model", 0);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
CouplerSolidPhaseField::~CouplerSolidPhaseField() {}

/* -------------------------------------------------------------------------- */
void CouplerSolidPhaseField::initFullImpl(const ModelOptions & options) {

  Model::initFullImpl(options);

  this->initBC(*this, *displacement, *displacement_increment, *external_force);
}

/* -------------------------------------------------------------------------- */
void CouplerSolidPhaseField::initModel() {

  getFEEngine().initShapeFunctions(_not_ghost);
  getFEEngine().initShapeFunctions(_ghost);
}

/* -------------------------------------------------------------------------- */
FEEngine & CouplerSolidPhaseField::getFEEngineBoundary(const ID & name) {
  return dynamic_cast<FEEngine &>(
      getFEEngineClassBoundary<MyFEEngineType>(name));
}

/* -------------------------------------------------------------------------- */
void CouplerSolidPhaseField::initSolver(TimeStepSolverType,
                                        NonLinearSolverType) {}

/* -------------------------------------------------------------------------- */
std::tuple<ID, TimeStepSolverType>
CouplerSolidPhaseField::getDefaultSolverID(const AnalysisMethod & method) {

  return std::make_tuple("unkown", TimeStepSolverType::_not_defined);
}

/* -------------------------------------------------------------------------- */
TimeStepSolverType CouplerSolidPhaseField::getDefaultSolverType() const {
  return TimeStepSolverType::_dynamic_lumped;
}

/* -------------------------------------------------------------------------- */
ModelSolverOptions CouplerSolidPhaseField::getDefaultSolverOptions(
    const TimeStepSolverType & type) const {
  ModelSolverOptions options;

  switch (type) {
  case TimeStepSolverType::_dynamic_lumped: {
    options.non_linear_solver_type = NonLinearSolverType::_lumped;
    options.integration_scheme_type["displacement"] =
        IntegrationSchemeType::_central_difference;
    options.solution_type["displacement"] = IntegrationScheme::_acceleration;
    break;
  }
  case TimeStepSolverType::_dynamic: {
    options.non_linear_solver_type = NonLinearSolverType::_lumped;
    options.integration_scheme_type["displacement"] =
        IntegrationSchemeType::_central_difference;
    options.solution_type["displacement"] = IntegrationScheme::_acceleration;
    break;
  }
  case TimeStepSolverType::_static: {
    options.non_linear_solver_type = NonLinearSolverType::_linear;
    options.integration_scheme_type["displacement"] =
        IntegrationSchemeType::_pseudo_time;
    options.solution_type["displacement"] = IntegrationScheme::_not_defined;
    break;
  }
  default:
    AKANTU_EXCEPTION(type << " is not a valid time step solver type");
    break;
  }

  return options;
}

/* -------------------------------------------------------------------------- */
void CouplerSolidPhaseField::assembleResidual() {

  // computes the internal forces
  this->assembleInternalForces();

  auto & solid_internal_force = solid->getInternalForce();
  auto & solid_external_force = solid->getExternalForce();

  auto & phasefield_internal_force = phasefield->getInternalForce();
  auto & phasefield_external_force = phasefield->getExternalForce();

  /* ------------------------------------------------------------------------ */
  this->getDOFManager().assembleToResidual("displacement", solid_external_force,
                                           1);
  this->getDOFManager().assembleToResidual("displacement", solid_internal_force,
                                           1);
  this->getDOFManager().assembleToResidual("damage", phasefield_external_force,
                                           1);
  this->getDOFManager().assembleToResidual("damage", phasefield_internal_force,
                                           1);
}

/* -------------------------------------------------------------------------- */
void CouplerSolidPhaseField::assembleResidual(const ID & residual_part) {
  AKANTU_DEBUG_IN();

  auto & solid_internal_force = solid->getInternalForce();
  auto & solid_external_force = solid->getExternalForce();

  auto & phasefield_internal_force = phasefield->getInternalForce();
  auto & phasefield_external_force = phasefield->getExternalForce();

  if ("external" == residual_part) {
    this->getDOFManager().assembleToResidual("displacement",
                                             solid_external_force, 1);
    this->getDOFManager().assembleToResidual("displacement",
                                             solid_internal_force, 1);
    AKANTU_DEBUG_OUT();
    return;
  }

  if ("internal" == residual_part) {
    this->getDOFManager().assembleToResidual("damage",
                                             phasefield_external_force, 1);
    this->getDOFManager().assembleToResidual("damage",
                                             phasefield_internal_force, 1);
    AKANTU_DEBUG_OUT();
    return;
  }

  AKANTU_CUSTOM_EXCEPTION(
      debug::SolverCallbackResidualPartUnknown(residual_part));

  AKANTU_DEBUG_OUT();
}



/* -------------------------------------------------------------------------- */
MatrixType CouplerSolidPhaseField::getMatrixType(const ID & matrix_id) {

  if (matrix_id == "K")
    return _symmetric;
  if (matrix_id == "M") {
    return _symmetric;
  }
  return _mt_not_defined;
}

/* -------------------------------------------------------------------------- */
void CouplerSolidPhaseField::assembleMatrix(const ID & matrix_id) {

  if (matrix_id == "K") {
    this->assembleStiffnessMatrix();
  } else if (matrix_id == "M") {
    solid->assembleMass();
  }
}

/* -------------------------------------------------------------------------- */
void CouplerSolidPhaseField::assembleLumpedMatrix(const ID & matrix_id) {
  if (matrix_id == "M") {
    solid->assembleMassLumped();
  }
}

/* -------------------------------------------------------------------------- */
void CouplerSolidPhaseField::assembleInternalForces() {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_INFO("Assemble the internal forces");

  solid->assembleInternalForces();
  phasefield->assembleInternalForces();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void CouplerSolidPhaseField::assembleStiffnessMatrix() {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_INFO("Assemble the new stiffness matrix");

  solid->assembleStiffnessMatrix();
  phasefield->assembleStiffnessMatrix();
  
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void CouplerSolidPhaseField::assembleMassLumped() {
  solid->assembleMassLumped();
}

/* -------------------------------------------------------------------------- */
void CouplerSolidPhaseField::assembleMass() { solid->assembleMass(); }

/* -------------------------------------------------------------------------- */
void CouplerSolidPhaseField::assembleMassLumped(GhostType ghost_type) {
  //solid->assembleMassLumped(ghost_type);
}

/* -------------------------------------------------------------------------- */
void CouplerSolidPhaseField::assembleMass(GhostType ghost_type) {
  //solid->assembleMass(ghost_type);
}

/* ------------------------------------------------------------------------- */
void CouplerSolidPhaseField::computeDamageOnQuadPoints(
    const GhostType & ghost_type) {

  AKANTU_DEBUG_IN();

  auto & fem = phasefield->getFEEngine();
  auto & mesh = phasefield->getMesh();

  switch (spatial_dimension) {
  case 1: {
    auto & mat = static_cast<MaterialPhaseField<1> &>(solid->getMaterial(0));
    auto & damage = mat.getDamage();
    for (auto & type :
         mesh.elementTypes(Model::spatial_dimension, ghost_type)) {
      auto & damage_on_qpoints_vect = damage(type, ghost_type);
      fem.interpolateOnIntegrationPoints(
          phasefield->getDamage(), damage_on_qpoints_vect, 1, type, ghost_type);
    }

    break;
  }
  case 2: {
    auto & mat = static_cast<MaterialPhaseField<2> &>(solid->getMaterial(0));
    auto & damage = mat.getDamage();

    for (auto & type :
         mesh.elementTypes(Model::spatial_dimension, ghost_type)) {
      auto & damage_on_qpoints_vect = damage(type, ghost_type);
      fem.interpolateOnIntegrationPoints(
          phasefield->getDamage(), damage_on_qpoints_vect, 1, type, ghost_type);
    }
    break;
  }
  default:
    auto & mat = static_cast<MaterialPhaseField<3> &>(solid->getMaterial(0));
    break;
  }

  AKANTU_DEBUG_OUT();
}

/* ------------------------------------------------------------------------- */
void CouplerSolidPhaseField::computeStrainOnQuadPoints(
    const GhostType & ghost_type) {
  AKANTU_DEBUG_IN();

  auto & strain_on_qpoints = phasefield->getStrain();
  auto & gradu_on_qpoints = solid->getMaterial(0).getGradU();

  for (auto & type : mesh.elementTypes(Model::spatial_dimension, ghost_type)) {
    auto & strain_on_qpoints_vect = strain_on_qpoints(type, ghost_type);
    auto & gradu_on_qpoints_vect = gradu_on_qpoints(type, ghost_type);
    for (auto && values :
         zip(make_view(strain_on_qpoints_vect, Model::spatial_dimension,
                       Model::spatial_dimension),
             make_view(gradu_on_qpoints_vect, Model::spatial_dimension,
                       Model::spatial_dimension))) {
      auto & strain = std::get<0>(values);
      auto & grad_u = std::get<1>(values);
      this->gradUToEpsilon(grad_u, strain);
    }
  }

  AKANTU_DEBUG_OUT();
}

/* ------------------------------------------------------------------------- */
void CouplerSolidPhaseField::solve() {

  solid->solveStep();
  this->computeStrainOnQuadPoints(_not_ghost);

  phasefield->solveStep();
  this->computeDamageOnQuadPoints(_not_ghost);

  solid->assembleInternalForces();
}

/* ------------------------------------------------------------------------- */
void CouplerSolidPhaseField::gradUToEpsilon(const Matrix<Real> & grad_u,
                                            Matrix<Real> & epsilon) {
  for (UInt i = 0; i < Model::spatial_dimension; ++i) {
    for (UInt j = 0; j < Model::spatial_dimension; ++j)
      epsilon(i, j) = 0.5 * (grad_u(i, j) + grad_u(j, i));
  }
}

/* ------------------------------------------------------------------------- */
bool CouplerSolidPhaseField::checkConvergence(Array<Real> & u_new,
                                              Array<Real> & u_old,
                                              Array<Real> & d_new,
                                              Array<Real> & d_old) {

  const Array<bool> & blocked_dofs = solid->getBlockedDOFs();
  UInt nb_degree_of_freedom = u_new.size();

  auto u_n_it = u_new.begin();
  auto u_o_it = u_old.begin();
  auto bld_it = blocked_dofs.begin();

  Real norm = 0;
  for (UInt n = 0; n < nb_degree_of_freedom;
       ++n, ++u_n_it, ++u_o_it, ++bld_it) {
    if ((!*bld_it)) {
      norm += (*u_n_it - *u_o_it) * (*u_n_it - *u_o_it);
    }
  }

  norm = std::sqrt(norm);

  auto d_n_it = d_new.begin();
  auto d_o_it = d_old.begin();
  nb_degree_of_freedom = d_new.size();

  Real norm2 = 0;
  for (UInt i = 0; i < nb_degree_of_freedom; ++i) {
    norm2 += (*d_n_it - *d_o_it);
  }

  norm2 = std::sqrt(norm2);

  Real error = std::max(norm, norm2);

  Real tolerance = 1e-8;
  if (error < tolerance) {

    return true;
  }

  return false;
}

/* -------------------------------------------------------------------------- */
#ifdef AKANTU_USE_IOHELPER

/* -------------------------------------------------------------------------- */
std::shared_ptr<dumper::Field> CouplerSolidPhaseField::createElementalField(
    const std::string & field_name, const std::string & group_name,
    bool padding_flag, const UInt & spatial_dimension,
    const ElementKind & kind) {

  return solid->createElementalField(field_name, group_name, padding_flag,
                                     spatial_dimension, kind);

  std::shared_ptr<dumper::Field> field;
  return field;
}

/* -------------------------------------------------------------------------- */
std::shared_ptr<dumper::Field>
CouplerSolidPhaseField::createNodalFieldReal(const std::string & field_name,
                                             const std::string & group_name,
                                             bool padding_flag) {

  return solid->createNodalFieldReal(field_name, group_name, padding_flag);

  std::shared_ptr<dumper::Field> field;
  return field;
}

/* -------------------------------------------------------------------------- */
std::shared_ptr<dumper::Field>
CouplerSolidPhaseField::createNodalFieldBool(const std::string & field_name,
                                             const std::string & group_name,
                                             bool padding_flag) {

  return solid->createNodalFieldBool(field_name, group_name, padding_flag);

  std::shared_ptr<dumper::Field> field;
  return field;
}

#else

/* -------------------------------------------------------------------------- */
std::shared_ptr<dumper::Field> CouplerSolidPhaseField::createElementalField(
    const std::string &, const std::string &, bool, const UInt &,
    const ElementKind &) {
  return nullptr;
}

/* ----------------------------------------------------------------------- */
std::shared_ptr<dumper::Field>
CouplerSolidPhaseField::createNodalFieldReal(const std::string &,
                                             const std::string &, bool) {
  return nullptr;
}

/*-------------------------------------------------------------------*/
std::shared_ptr<dumper::Field>
CouplerSolidPhaseField::createNodalFieldBool(const std::string &,
                                             const std::string &, bool) {
  return nullptr;
}

#endif

/* -----------------------------------------------------------------------*/
void CouplerSolidPhaseField::dump(const std::string & dumper_name) {
  solid->onDump();
  mesh.dump(dumper_name);
}

/* ------------------------------------------------------------------------*/
void CouplerSolidPhaseField::dump(const std::string & dumper_name, UInt step) {
  solid->onDump();
  mesh.dump(dumper_name, step);
}

/* ----------------------------------------------------------------------- */
void CouplerSolidPhaseField::dump(const std::string & dumper_name, Real time,
                                  UInt step) {
  solid->onDump();
  mesh.dump(dumper_name, time, step);
}

/* -------------------------------------------------------------------------- */
void CouplerSolidPhaseField::dump() {
  solid->onDump();
  mesh.dump();
}

/* -------------------------------------------------------------------------- */
void CouplerSolidPhaseField::dump(UInt step) {
  solid->onDump();
  mesh.dump(step);
}

/* -------------------------------------------------------------------------- */
void CouplerSolidPhaseField::dump(Real time, UInt step) {
  solid->onDump();
  mesh.dump(time, step);
}

} // namespace akantu

/* -------------------------------------------------------------------------- */
#include "non_linear_solver.hh"
#include "phase_field_model.hh"
#include "solid_mechanics_model.hh"
#include "coupler_solid_phasefield.hh"
/* -------------------------------------------------------------------------- */
#include <iostream>
/* -------------------------------------------------------------------------- */

using namespace akantu;
const UInt spatial_dimension = 2;
void computeStrainOnQuadPoints(SolidMechanicsModel &, PhaseFieldModel &,
                               const GhostType &);
void computeDamageOnQuadPoints(SolidMechanicsModel &, PhaseFieldModel &,
                               const GhostType &);
void gradUToEpsilon(const Matrix<Real> &, Matrix<Real> &);
bool testConvergence(SolidMechanicsModel &, PhaseFieldModel &, Array<Real> &,
                     Array<Real> &, Array<Real> &, Array<Real> &);

/* -------------------------------------------------------------------------- */

int main(int argc, char *argv[])
{
  initialize("material_notch.dat", argc, argv);

  // create mesh
  Mesh mesh(spatial_dimension);
  mesh.read("square_notch.msh");

  PhaseFieldModel phasefield(mesh);
  phasefield.initFull(_analysis_method = _static);

  SolidMechanicsModel solid(mesh);
  solid.initFull(_analysis_method = _static);
 
  solid.applyBC(BC::Dirichlet::FixedValue(0., _y), "bottom");
  solid.applyBC(BC::Dirichlet::FixedValue(0., _x), "left");

  solid.setBaseName(        "square_notch_modified");
  solid.addDumpFieldVector( "displacement");
  solid.addDumpFieldVector( "internal_force");
  solid.addDumpField(       "stress");
  solid.addDumpField(       "grad_u");
  solid.addDumpField(       "damage");
  solid.addDumpField(       "blocked_dofs");
  solid.dump();  

  UInt nbSteps   = 1000;
  Real increment = 1.e-5;

  solid.getNewSolver("solid_linear", TimeStepSolverType::_static,
                     NonLinearSolverType::_linear);
  solid.setIntegrationScheme("solid_linear", "displacement",
			     IntegrationSchemeType::_pseudo_time);

  phasefield.getNewSolver("phase_linear", TimeStepSolverType::_static,
			  NonLinearSolverType::_linear);
  phasefield.setIntegrationScheme("phase_linear", "damage",
			     IntegrationSchemeType::_pseudo_time);

  solid.setDefaultSolver("solid_linear");
  phasefield.setDefaultSolver("phase_linear");
  
  for (UInt s = 1; s < nbSteps; ++s) {
    solid.applyBC(BC::Dirichlet::IncrementValue(increment, _y), "top");
    solid.solveStep();
    computeStrainOnQuadPoints(solid, phasefield, _not_ghost);

    phasefield.solveStep();
    computeDamageOnQuadPoints(solid, phasefield, _not_ghost);

    if (s % 50 == 0) {
      solid.dump();
    }
    
    std::cout << "Step " << s << "/" << nbSteps << std::endl;
    
  }

  finalize();
  return EXIT_SUCCESS;
}

/* -------------------------------------------------------------------------- */
bool testConvergence(SolidMechanicsModel &solid, PhaseFieldModel &phase,
                     Array<Real> &u_new, Array<Real> &u_old, Array<Real> &d_new,
                     Array<Real> &d_old) {
  Real tolerance = 1e-8;
  Real norm_u = 0;
  Real norm_d = 0;

  UInt nb_degree_of_freedom = u_new.size();

  const Array<bool> &blocked_dofs = solid.getBlockedDOFs();

  for (auto &&values : zip(make_view(blocked_dofs, 1), make_view(u_new, 1),
                           make_view(u_old, 1))) {
    auto &bld = std::get<0>(values);
    auto &u_n = std::get<1>(values);
    auto &u_o = std::get<2>(values);
    if (!bld[0]) {
      norm_u += (u_n[0] - u_o[0]) * (u_n[0] - u_o[0]);
    }
  }

  auto d_n_it = d_new.begin();
  auto d_o_it = d_old.begin();

  nb_degree_of_freedom = d_new.size();

  for (UInt i = 0; i < nb_degree_of_freedom; ++i) {
    norm_d += (*d_n_it - *d_o_it);
  }

  norm_u = std::sqrt(norm_u);
  norm_d = std::sqrt(norm_d);

  std::cerr << norm_u << "--------" << norm_d << std::endl;

  Real error = std::max(norm_u, norm_d);

  if (error < tolerance) {
    return true;
  }

  return false;
}
/* -------------------------------------------------------------------------- */
void computeStrainOnQuadPoints(SolidMechanicsModel &solid,
                               PhaseFieldModel &phase,
                               const GhostType &ghost_type) {
  AKANTU_DEBUG_IN();

  auto &mesh = solid.getMesh();

  auto &strain_on_qpoints = phase.getStrain();
  auto &gradu_on_qpoints = solid.getMaterial(0).getGradU();

  for (auto &type : mesh.elementTypes(spatial_dimension, ghost_type)) {
    auto &strain_on_qpoints_vect = strain_on_qpoints(type, ghost_type);
    auto &gradu_on_qpoints_vect = gradu_on_qpoints(type, ghost_type);
    for (auto &&values : zip(make_view(strain_on_qpoints_vect,
                                       spatial_dimension, spatial_dimension),
                             make_view(gradu_on_qpoints_vect, spatial_dimension,
                                       spatial_dimension))) {
      auto &strain = std::get<0>(values);
      auto &grad_u = std::get<1>(values);
      gradUToEpsilon(grad_u, strain);
    }
  }
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void computeDamageOnQuadPoints(SolidMechanicsModel &solid,
                               PhaseFieldModel &phase,
                               const GhostType &ghost_type) {
  AKANTU_DEBUG_IN();

  auto &fem = phase.getFEEngine();
  auto &mesh = phase.getMesh();

  switch (spatial_dimension) {
  case 1: {
    auto &mat = static_cast<MaterialPhaseField<1> &>(solid.getMaterial(0));
    auto &damage = mat.getDamage();
    for (auto &type : mesh.elementTypes(spatial_dimension, ghost_type)) {
      auto &damage_on_qpoints_vect = damage(type, ghost_type);
      fem.interpolateOnIntegrationPoints(
          phase.getDamage(), damage_on_qpoints_vect, 1, type, ghost_type);
    }

    break;
  }
  case 2: {
    auto &mat = static_cast<MaterialPhaseField<2> &>(solid.getMaterial(0));
    auto &damage = mat.getDamage();

    for (auto &type : mesh.elementTypes(spatial_dimension, ghost_type)) {
      auto &damage_on_qpoints_vect = damage(type, ghost_type);
      fem.interpolateOnIntegrationPoints(
          phase.getDamage(), damage_on_qpoints_vect, 1, type, ghost_type);
    }
    break;
  }
  default:
    auto &mat = static_cast<MaterialPhaseField<3> &>(solid.getMaterial(0));
    break;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void gradUToEpsilon(const Matrix<Real> &grad_u, Matrix<Real> &epsilon) {
  for (UInt i = 0; i < spatial_dimension; ++i) {
    for (UInt j = 0; j < spatial_dimension; ++j)
      epsilon(i, j) = 0.5 * (grad_u(i, j) + grad_u(j, i));
  }
}

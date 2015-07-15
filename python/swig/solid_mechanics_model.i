%{
  #include "solid_mechanics_model.hh"
  #include "sparse_matrix.hh"
%}

namespace akantu {
  %ignore SolidMechanicsModel::initFEEngineBoundary;
  %ignore SolidMechanicsModel::initParallel;
  %ignore SolidMechanicsModel::initArrays;
  %ignore SolidMechanicsModel::initMaterials;
  %ignore SolidMechanicsModel::initModel;
  %ignore SolidMechanicsModel::initPBC;


  %ignore SolidMechanicsModel::initExplicit;
  %ignore SolidMechanicsModel::isExplicit;
  %ignore SolidMechanicsModel::updateCurrentPosition;
  %ignore SolidMechanicsModel::updateAcceleration;
  %ignore SolidMechanicsModel::updateIncrement;
  %ignore SolidMechanicsModel::updatePreviousDisplacement;
  %ignore SolidMechanicsModel::saveStressAndStrainBeforeDamage;
  %ignore SolidMechanicsModel::updateEnergiesAfterDamage;
  %ignore SolidMechanicsModel::solveLumped;
  %ignore SolidMechanicsModel::explicitPred;
  %ignore SolidMechanicsModel::explicitCorr;

  %ignore SolidMechanicsModel::initSolver;
  %ignore SolidMechanicsModel::initImplicit;
  %ignore SolidMechanicsModel::initialAcceleration;


  %ignore SolidMechanicsModel::testConvergence;
  %ignore SolidMechanicsModel::testConvergenceIncrement;
  %ignore SolidMechanicsModel::testConvergenceResidual;
  %ignore SolidMechanicsModel::initVelocityDampingMatrix;

  %ignore SolidMechanicsModel::getNbDataForElements;
  %ignore SolidMechanicsModel::packElementData;
  %ignore SolidMechanicsModel::unpackElementData;
  %ignore SolidMechanicsModel::getNbDataToPack;
  %ignore SolidMechanicsModel::getNbDataToUnpack;
  %ignore SolidMechanicsModel::packData;
  %ignore SolidMechanicsModel::unpackData;

  %ignore SolidMechanicsModel::setMaterialSelector;
  %ignore SolidMechanicsModel::getSolver;
  %ignore SolidMechanicsModel::getSynchronizer;

  %ignore Dumpable::registerExternalDumper;

}

%template(SolidMechanicsBoundaryCondition) akantu::BoundaryCondition<akantu::SolidMechanicsModel>;

%include "dumpable.hh"

print_self(SolidMechanicsModel)

%include "solid_mechanics_model.hh"

%template(testConvergenceResidual) akantu::SolidMechanicsModel::testConvergence<akantu::SolveConvergenceCriteria::_scc_residual>;

%extend akantu::SolidMechanicsModel {
  
  void solveStaticDisplacement(Real tolerance, UInt max_iteration) {

    $self->solveStatic<akantu::SolveConvergenceMethod::_scm_newton_raphson_tangent_not_computed, akantu::SolveConvergenceCriteria::_scc_residual>(tolerance, max_iteration);

  }

  void solveDisplCorr(bool need_factorize, bool has_profile_changed) {

    akantu::Array<akantu::Real> & increment = $self->getIncrement();

    $self->solve<akantu::IntegrationScheme2ndOrder::_displacement_corrector>(increment, 1., need_factorize, has_profile_changed);

  }
  
  void clearDispl() {

    akantu::Array<akantu::Real> & displ = $self->getDisplacement();
    displ.clear();
  }

}

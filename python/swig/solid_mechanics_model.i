%{
  #include "solid_mechanics_model.hh"
  #include "sparse_matrix.hh"
  #include "boundary_condition.hh"
  #include "boundary_condition_functor.hh"
  #include "boundary_condition_python_functor.hh"
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

%extend akantu::SolidMechanicsModel {
  

  bool testConvergenceSccRes(Real tolerance) {

    Real error = 0;
    bool res = self->testConvergence<akantu::SolveConvergenceCriteria::_scc_residual>(tolerance,error);
    return res;
  }

  void solveStaticDisplacement(Real tolerance, UInt max_iteration) {

    $self->solveStatic<akantu::SolveConvergenceMethod::_scm_newton_raphson_tangent_not_computed, akantu::SolveConvergenceCriteria::_scc_residual>(tolerance, max_iteration);

  }

  void applyDirichletBC(PyObject * func_obj,
			const std::string & group_name) {
    
    akantu::BC::PythonFunctorDirichlet functor(func_obj);
    $self->applyBC(functor,group_name);
  }

  void applyNeumannBC(PyObject * func_obj,
		      const std::string & group_name) {
    
    akantu::BC::PythonFunctorNeumann functor(func_obj);
    $self->applyBC(functor,group_name);
  }

  
  void solveDisplCorr(bool need_factorize, bool has_profile_changed) {

    akantu::Array<akantu::Real> & increment = $self->getIncrement();

    $self->solve<akantu::IntegrationScheme2ndOrder::_displacement_corrector>(increment, 1., need_factorize, has_profile_changed);

  }
  
  void clearDispl() {

    akantu::Array<akantu::Real> & displ = $self->getDisplacement();
    displ.clear();
  }

  void solveStep_TgModifIncr(Real tolerance, UInt max_iteration) {
      
  $self->solveStep<akantu::SolveConvergenceMethod::_scm_newton_raphson_tangent_modified, akantu::SolveConvergenceCriteria::_scc_residual>(tolerance, max_iteration);
    
 }
  
  void clearDisplVeloAcc() {

  akantu::Array<akantu::Real> & displ = $self->getDisplacement();
  akantu::Array<akantu::Real> & velo = $self->getVelocity();
  akantu::Array<akantu::Real> & acc = $self->getAcceleration();

  displ.clear();
  velo.clear();
  acc.clear();
 }

  void applyUniformPressure(Real pressure, const std::string surface_name){
  
  UInt spatial_dimension = $self->getSpatialDimension();
  akantu::Matrix<akantu::Real> surface_stress(spatial_dimension, spatial_dimension, 0.0);

  for(UInt i = 0; i < spatial_dimension; ++i) {
  surface_stress(i,i) = -pressure;
 }  
  $self->applyBC(akantu::BC::Neumann::FromStress(surface_stress), surface_name);
 }

  void blockDOF(const std::string surface_name, SpacialDirection direction){

    $self->applyBC(akantu::BC::Dirichlet::FixedValue(0.0, direction), surface_name);
  }
}

%{
  #include "structural_mechanics_model.hh"
  #include "sparse_matrix.hh"
%}

namespace akantu {
  %ignore StructuralMechanicsModel::initArrays;
  %ignore StructuralMechanicsModel::initModel;
  
  %ignore StructuralMechanicsModel::initSolver;
  %ignore StructuralMechanicsModel::initImplicit;

  static void lin_load(double * position, double * load, 
		       __attribute__ ((unused)) Real * normal,
		       __attribute__ ((unused)) UInt surface_id) {
      
    memset(load,0,sizeof(Real)*3);
    if (position[0]<=10){
      load[1]= -6000;
    }
  }  
}

%include "dumpable.hh"
%include "structural_mechanics_model.hh"

%extend akantu::StructuralMechanicsModel {
  
  bool solveStep(Real tolerance, UInt max_iteration) {

    return $self->solveStep<akantu::SolveConvergenceMethod::_scm_newton_raphson_tangent, akantu::SolveConvergenceCriteria::_scc_residual>(tolerance, max_iteration);

  }

  void computeForcesFromFunctionBeam2d(BoundaryFunctionType function_type) {

    $self->computeForcesFromFunction<akantu::ElementType::_bernoulli_beam_2>(akantu::lin_load, function_type);

  }
 }

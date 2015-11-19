%{
  #include "boundary_condition_python_functor.hh"
%}


namespace akantu {
  %ignore Model::createSynchronizerRegistry;
  %ignore Model::createParallelSynch;
  %ignore Model::getDOFSynchronizer;
  //%ignore Model::getSynchronizerRegistry;
  %ignore Model::registerFEEngineObject;
  %ignore Model::unregisterFEEngineObject;
  %ignore Model::getFEEngineBoundary;
  //  %ignore Model::getFEEngine;
  %ignore Model::getFEEngineClass;
  %ignore Model::getFEEngineClassBoundary;
  %ignore Model::setParser;

  %ignore IntegrationPoint::operator=;

  %ignore FEEngine::getNbIntegrationPoints;
  %ignore FEEngine::getShapes;
  %ignore FEEngine::getShapesDerivatives;
  %ignore FEEngine::getIntegrationPoints;
  %ignore FEEngine::getIGFEMElementTypes;
  %ignore FEEngine::interpolateOnIntegrationPoints(const Array<Real> &,ElementTypeMapArray<Real> &,const ElementTypeMapArray<UInt> *) const;
  %ignore FEEngine::interpolateOnIntegrationPoints(const Array<Real> &,ElementTypeMapArray<Real> &) const;
  %ignore FEEngine::interpolateOnIntegrationPoints(const Array<Real> &,Array<Real> &,UInt,const ElementType&,const GhostType &,const Array< UInt > &) const;
  %ignore FEEngine::interpolateOnIntegrationPoints(const Array<Real> &,Array<Real> &,UInt,const ElementType&,const GhostType &) const;

}

%include "sparse_matrix.hh"
%include "fe_engine.hh"

%rename(FreeBoundaryDirichlet) akantu::BC::Dirichlet::FreeBoundary;
%rename(FreeBoundaryNeumann) akantu::BC::Neumann::FreeBoundary;
%rename(PythonBoundary) akantu::BC::Dirichlet::PythonFunctor;

%include "boundary_condition_functor.hh"
%include "boundary_condition.hh"
%include "boundary_condition_python_functor.hh"
%include "communication_buffer.hh"
%include "data_accessor.hh"
%include "synchronizer.hh"
%include "synchronizer_registry.hh"
%include "model.hh"


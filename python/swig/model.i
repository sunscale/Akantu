namespace akantu {
  %ignore Model::createSynchronizerRegistry;
  %ignore Model::createParallelSynch;
  %ignore Model::getDOFSynchronizer;
  %ignore Model::getSynchronizerRegistry;
  %ignore Model::registerFEEngineObject;
  %ignore Model::unregisterFEEngineObject;
  %ignore Model::getFEEngineBoundary;
  //  %ignore Model::getFEEngine;
  %ignore Model::getFEEngineClass;
  %ignore Model::getFEEngineClassBoundary;
  %ignore Model::setParser;

  %ignore QuadraturePoint::operator=;

  %ignore FEEngine::getNbQuadraturePoints;
  %ignore FEEngine::getShapes;
  %ignore FEEngine::getShapesDerivatives;
  %ignore FEEngine::getQuadraturePoints;
  %ignore FEEngine::getIGFEMElementTypes;
  %ignore FEEngine::interpolateOnQuadraturePoints(const Array<Real> &,ElementTypeMapArray<Real> &,const ElementTypeMapArray<UInt> *) const;
  %ignore FEEngine::interpolateOnQuadraturePoints(const Array<Real> &,ElementTypeMapArray<Real> &) const;
  %ignore FEEngine::interpolateOnQuadraturePoints(const Array<Real> &,Array<Real> &,UInt,const ElementType&,const GhostType &,const Array< UInt > &) const;
  %ignore FEEngine::interpolateOnQuadraturePoints(const Array<Real> &,Array<Real> &,UInt,const ElementType&,const GhostType &) const;
}

%include "sparse_matrix.hh"
%include "fe_engine.hh"

%rename(FreeBoundaryDirichlet) akantu::BC::Dirichlet::FreeBoundary;
%rename(FreeBoundaryNeumann) akantu::BC::Neumann::FreeBoundary;



%include "boundary_condition_functor.hh"

%include "boundary_condition.hh"

%include "model.hh"

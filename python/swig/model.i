/**
 * @file   model.i
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Dec 12 2014
 * @date last modification: Wed Nov 11 2015
 *
 * @brief  model wrapper
 *
 * @section LICENSE
 *
 * Copyright (©) 2015 EPFL (Ecole Polytechnique Fédérale de Lausanne) Laboratory
 * (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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

%{
  #include "boundary_condition_python_functor.hh"
  #include "model_solver.hh"
  #include "non_linear_solver.hh"
  %}


namespace akantu {
  %ignore Model::createSynchronizerRegistry;
  %ignore Model::getSynchronizerRegistry;
  %ignore Model::createParallelSynch;
  %ignore Model::getDOFSynchronizer;
  %ignore Model::registerFEEngineObject;
  %ignore Model::unregisterFEEngineObject;
  %ignore Model::getFEEngineBoundary;
  //  %ignore Model::getFEEngine;
  %ignore Model::getFEEngineClass;
  %ignore Model::getFEEngineClassBoundary;
  %ignore Model::setParser;
  %ignore Model::updateDataForNonLocalCriterion;
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
  %ignore FEEngine::onNodesAdded;
  %ignore FEEngine::onNodesRemoved;
  %ignore FEEngine::onElementsAdded;
  %ignore FEEngine::onElementsChanged;
  %ignore FEEngine::onElementsRemoved;
  %ignore FEEngine::elementTypes;
}

%include "sparse_matrix.i"
%include "fe_engine.hh"

%rename(FreeBoundaryDirichlet) akantu::BC::Dirichlet::FreeBoundary;
%rename(FreeBoundaryNeumann) akantu::BC::Neumann::FreeBoundary;
%rename(PythonBoundary) akantu::BC::Dirichlet::PythonFunctor;

%include "boundary_condition_functor.hh"
%include "boundary_condition.hh"
%include "boundary_condition_python_functor.hh"
%include "communication_buffer.hh"
%include "data_accessor.hh"
//%include "synchronizer.hh"
//%include "synchronizer_registry.hh"
%include "model.hh"
%include "non_linear_solver.hh"

%extend akantu::Model {

  void solveStep(){
    $self->solveStep();
  }

  akantu::NonLinearSolver & getNonLinearSolver(){
   return $self->getNonLinearSolver();
  }
 }

%extend akantu::NonLinearSolver {

  void set(const std::string & id, akantu::Real val){
    if (id == "max_iterations")
      $self->set(id, int(val));
    else if (id == "convergence_type")
      $self->set(id, akantu::SolveConvergenceCriteria(UInt(val)));
    else $self->set(id, val);
  }
 }

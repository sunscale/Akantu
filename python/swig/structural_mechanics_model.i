/**
 * @file   structural_mechanics_model.i
 *
 * @author Fabian Barras <fabian.barras@epfl.ch>
 *
 * @date creation: Wed Apr 01 2015
 *
 * @brief  structural mechanics model wrapper
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

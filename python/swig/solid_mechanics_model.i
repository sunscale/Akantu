/**
 * @file   solid_mechanics_model.i
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Dec 12 2014
 * @date last modification: Wed Jan 06 2016
 *
 * @brief  solid mechanics model wrapper
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
  #include "solid_mechanics_model.hh"
  #include "sparse_matrix.hh"
  #include "boundary_condition.hh"
  #include "boundary_condition_functor.hh"
  #include "boundary_condition_python_functor.hh"
  #include "material_selector.hh"
  #include "material_python.hh"
%}

namespace akantu {
  %ignore SolidMechanicsModel::initFEEngineBoundary;
  %ignore SolidMechanicsModel::initParallel;
  %ignore SolidMechanicsModel::initArrays;
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
  %ignore Material::onNodesAdded;
  %ignore Material::onNodesRemoved;
  %ignore Material::onElementsAdded;
  %ignore Material::onElementsRemoved;
  %ignore Material::onElementsChanged;
}

%template(SolidMechanicsBoundaryCondition) akantu::BoundaryCondition<akantu::SolidMechanicsModel>;

%include "dumpable.hh"

print_self(SolidMechanicsModel)

%include "material.i"
%include "solid_mechanics_model.hh"


%extend akantu::SolidMechanicsModel {

  /* ------------------------------------------------------------------------ */
  void setPhysicalNamesMaterialSelector(){
    akantu::MeshDataMaterialSelector<std::string> * selector = new
      akantu::MeshDataMaterialSelector<std::string>("physical_names", *self);
    self->setMaterialSelector(*selector);
  }
  
  /* ------------------------------------------------------------------------ */
  bool testConvergenceSccRes(Real tolerance) {
    Real error = 0;
    bool res = self->testConvergence<akantu::_scc_residual>(tolerance, error);
    return res;
  }

  /* ------------------------------------------------------------------------ */
  void solveStaticDisplacement(Real tolerance, UInt max_iteration) {
    $self->solveStatic<akantu::_scm_newton_raphson_tangent_not_computed,
                       akantu::_scc_residual>(tolerance, max_iteration);
  }

  /* ------------------------------------------------------------------------ */
  /// register an empty material of a given type
  void registerNewPythonMaterial(PyObject * obj, const akantu::ID & mat_type) {
    std::pair<akantu::Parser::const_section_iterator,
              akantu::Parser::const_section_iterator>
        sub_sect = akantu::getStaticParser().getSubSections(akantu::_st_material);

    akantu::Parser::const_section_iterator it = sub_sect.first;
    for (; it != sub_sect.second; ++it) {
      if (it->getName() == mat_type) {

        AKANTU_DEBUG_ASSERT($self->materials_names_to_id.find(mat_type) ==
                            $self->materials_names_to_id.end(),
                            "A material with this name '"
                            << mat_type << "' has already been registered. "
                            << "Please use unique names for materials");

        UInt mat_count = $self->materials.size();
        $self->materials_names_to_id[mat_type] = mat_count;

        std::stringstream sstr_mat;
        sstr_mat << $self->getID() << ":" << mat_count << ":" << mat_type;
        akantu::ID mat_id = sstr_mat.str();

        akantu::Material * material = new akantu::MaterialPython(*$self, obj, mat_id);
        $self->materials.push_back(material);

        material->parseSection(*it);
      }
    }
  }

  /* ------------------------------------------------------------------------ */
  void applyDirichletBC(PyObject * func_obj, const std::string & group_name) {
    akantu::BC::PythonFunctorDirichlet functor(func_obj);
    $self->applyBC(functor, group_name);
  }

  /* ------------------------------------------------------------------------ */
  void applyNeumannBC(PyObject * func_obj, const std::string & group_name) {
    akantu::BC::PythonFunctorNeumann functor(func_obj);
    $self->applyBC(functor, group_name);
  }

  /* ------------------------------------------------------------------------ */
  void solveDisplCorr(bool need_factorize, bool has_profile_changed) {
    akantu::Array<akantu::Real> & increment = $self->getIncrement();

    $self->solve<akantu::IntegrationScheme2ndOrder::_displacement_corrector>(
        increment, 1., need_factorize, has_profile_changed);
  }

  /* ------------------------------------------------------------------------ */
  void clearDispl() {
    akantu::Array<akantu::Real> & displ = $self->getDisplacement();
    displ.clear();
  }

  /* ------------------------------------------------------------------------ */
  void solveStep_TgModifIncr(Real tolerance, UInt max_iteration) {
    $self->solveStep<akantu::_scm_newton_raphson_tangent_modified,
                     akantu::_scc_increment>(tolerance, max_iteration);
  }

  /* ------------------------------------------------------------------------ */
  void solveStep_TgIncr(Real tolerance, Real & error, UInt max_iteration) {
    
    $self->solveStep<akantu::_scm_newton_raphson_tangent,
      akantu::_scc_increment>(tolerance, error, max_iteration);
  }

  /* ------------------------------------------------------------------------ */
  void clearDisplVeloAcc() {
    akantu::Array<akantu::Real> & displ = $self->getDisplacement();
    akantu::Array<akantu::Real> & velo = $self->getVelocity();
    akantu::Array<akantu::Real> & acc = $self->getAcceleration();

    displ.clear();
    velo.clear();
    acc.clear();
  }

  /* ------------------------------------------------------------------------ */
  void applyUniformPressure(Real pressure, const std::string surface_name) {
    UInt spatial_dimension = $self->getSpatialDimension();
    akantu::Matrix<akantu::Real> surface_stress(spatial_dimension,
                                                spatial_dimension, 0.0);

    for (UInt i = 0; i < spatial_dimension; ++i) {
      surface_stress(i, i) = -pressure;
    }
    $self->applyBC(akantu::BC::Neumann::FromStress(surface_stress),
                   surface_name);
  }

  /* ------------------------------------------------------------------------ */
  void blockDOF(const std::string surface_name, SpacialDirection direction) {
    $self->applyBC(akantu::BC::Dirichlet::FixedValue(0.0, direction),
                   surface_name);
  }
}

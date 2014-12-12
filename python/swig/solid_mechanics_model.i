%{
  #include "solid_mechanics_model.hh"
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
  %ignore SolidMechanicsModel::initializeUpdateResidualData;
  %ignore SolidMechanicsModel::updateCurrentPosition;
  %ignore SolidMechanicsModel::updateResidual;
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
  %ignore SolidMechanicsModel::assembleStiffnessMatrix;


  %ignore SolidMechanicsModel::testConvergence;
  %ignore SolidMechanicsModel::testConvergenceIncrement;
  %ignore SolidMechanicsModel::testConvergenceResidual;
  %ignore SolidMechanicsModel::initVelocityDampingMatrix;
  %ignore SolidMechanicsModel::implicitPred;
  %ignore SolidMechanicsModel::implicitCorr;

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


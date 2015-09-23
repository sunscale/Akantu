%{
  #include "heat_transfer_model.hh"
  #include "data_accessor.hh"
%}

namespace akantu {
  %ignore HeatTransferModel::initFEEngineBoundary;
  %ignore HeatTransferModel::initParallel;
  %ignore HeatTransferModel::initArrays;
  %ignore HeatTransferModel::initMaterials;
  %ignore HeatTransferModel::initModel;
  %ignore HeatTransferModel::initPBC;

  %ignore HeatTransferModel::initSolver;
  %ignore HeatTransferModel::initImplicit;

  %ignore HeatTransferModel::getNbDataForElements;
  %ignore HeatTransferModel::packElementData;
  %ignore HeatTransferModel::unpackElementData;
  %ignore HeatTransferModel::getNbDataToPack;
  %ignore HeatTransferModel::getNbDataToUnpack;
  %ignore HeatTransferModel::packData;
  %ignore HeatTransferModel::unpackData;

}

%include "heat_transfer_model.hh"

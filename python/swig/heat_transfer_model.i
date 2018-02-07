/**
 * @file   heat_transfer_model.i
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 *
 * @date creation: Wed Jul 15 2015
 *
 * @brief  heat transfer model wrapper
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

  %ignore HeatTransferModel::getNbDataToPack;
  %ignore HeatTransferModel::getNbData;
  %ignore HeatTransferModel::packData;
  %ignore HeatTransferModel::unpackData;

}

%include "heat_transfer_model.hh"


%extend akantu::HeatTransferModel {

  Real getParamReal(const ID & param){
     Real res = $self->get(param);
     return res;
   }
   UInt getParamUInt(const ID & param){
     UInt res = $self->get(param);
     return res;
   }
   int getParamInt(const ID & param){
     int res = $self->get(param);
     return res;
   }
}


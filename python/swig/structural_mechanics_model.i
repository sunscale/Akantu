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
%ignore StructuralMechanicsModel::onNodesAdded;
%ignore StructuralMechanicsModel::onNodesRemoved;
%ignore StructuralMechanicsModel::onElementsAdded;
%ignore StructuralMechanicsModel::onElementsRemoved;
%ignore StructuralMechanicsModel::onElementsChanged;
%ignore StructuralMechanicsModel::getNbData;
%ignore StructuralMechanicsModel::packData;
%ignore StructuralMechanicsModel::unpackData;
}

%include "dumpable.hh"
%include "structural_mechanics_model.hh"

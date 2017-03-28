/**
 * @file   material.i
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @brief  material wrapper
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
  #include "material_python.hh"
%}

namespace akantu {
  %ignore Material::onNodesAdded;
  %ignore Material::onNodesRemoved;
  %ignore Material::onElementsAdded;
  %ignore Material::onElementsRemoved;
  %ignore Material::onElementsChanged;
}

%include "material.hh"


%extend akantu::Material {
   Array<Real> & getArrayReal(const ID & id, const ElementType & type,
                              const GhostType & ghost_type = _not_ghost) {
      return $self->getArray<Real>(id, type, ghost_type);
   }
}
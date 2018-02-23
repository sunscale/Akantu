/**
 * @file   solid_mechanics_model_cohesive.i
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Mon Nov 23 2015
 * @date last modification: Wed Jan 13 2016
 *
 * @brief  
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
#include "cohesive_element_inserter.hh"
#include "solid_mechanics_model_cohesive.hh"
#include "material_cohesive.hh"
%}

namespace akantu {
%ignore SolidMechanicsModelCohesive::initFacetFilter;
%ignore SolidMechanicsModelCohesive::initParallel;
%ignore CohesiveElementInserter::initParallel;
}

%extend akantu::SolidMechanicsModelCohesive {
  Array<Real> & getRealInternalCohesiveField(const std::string field_name) {
    akantu::Mesh & mesh = $self->getMesh();
    akantu::ElementType type = *(mesh.firstType(mesh.getSpatialDimension(), akantu::_not_ghost, akantu::_ek_cohesive));
    return ($self->flattenInternal(field_name,akantu::_ek_cohesive, akantu::_not_ghost))(type);
  }
}

%include "cohesive_element_inserter.hh"
%include "solid_mechanics_model_cohesive.hh"

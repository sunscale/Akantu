/**
 * @file   fe_engine.i
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Dec 12 2014
 * @date last modification: Wed Nov 11 2015
 *
 * @brief  FEEngine wrapper
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
  #include "mesh.hh"
  %}

namespace akantu {
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

%extend akantu::FEEngine {
  void interpolateField(const Array<Real> & in, Array<Real> & out, ElementType type) {
    $self->interpolateOnIntegrationPoints(in, out, in.getNbComponent(), type);
  }
}

%include "sparse_matrix.i"
%include "fe_engine.hh"

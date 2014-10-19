/**
 * @file   structural_mechanics_model_mass.cc
 *
 * @author Sébastien Hartmann <sebastien.hartmann@epfl.ch>
 *
 * @date creation: Mon Jul 07 2014
 * @date last modification: Mon Jul 07 2014
 *
 * @brief  function handling mass computation
 *
 * @section LICENSE
 *
 * Copyright (©) 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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

/* -------------------------------------------------------------------------- */

#include "structural_mechanics_model.hh"
#include "material.hh"
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */

void StructuralMechanicsModel::assembleMass(){
  AKANTU_DEBUG_IN();

  if(!mass_matrix){
    std::stringstream sstr; sstr << id << ":mass_matrix";
    mass_matrix = new SparseMatrix(*jacobian_matrix, sstr.str(), memory_id);
  }

  assembleMass(_not_ghost);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void StructuralMechanicsModel::assembleMass(GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Array<Real> rho_1(0,1);
  mass_matrix->clear();

  Mesh::type_iterator it = mesh.firstType(spatial_dimension, ghost_type, _ek_structural);
  Mesh::type_iterator end = mesh.lastType(spatial_dimension, ghost_type, _ek_structural);
  for(; it != end; ++it){
    ElementType type = *it;

#define ASSEMBLE_MASS(type)		\
    assembleMass<type>();

    AKANTU_BOOST_STRUCTURAL_ELEMENT_SWITCH(ASSEMBLE_MASS);
#undef ASSEMBLE_MASS

  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

void StructuralMechanicsModel::computeRho(Array<Real> & rho,
					  ElementType type,
					  GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  UInt nb_element = getFEEngine().getMesh().getNbElement(type);
  UInt nb_quadrature_points = getFEEngine().getNbQuadraturePoints(type);

  Array<UInt> & el_mat = element_material(type, ghost_type);

  for (UInt e = 0; e < nb_element; ++e){
    UInt mat = el_mat(e);
    Real rho_el = materials[mat].rho;
    for (UInt q = e*nb_quadrature_points; q < e*nb_quadrature_points + nb_quadrature_points; ++q){
      rho(q) = rho_el;
    }
  }

  AKANTU_DEBUG_OUT();
}

__END_AKANTU__

/**
 * @file   solid_mechanics_model_mass.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Tue Oct 05 2010
 * @date last modification: Thu Jun 05 2014
 *
 * @brief  function handling mass computation
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "solid_mechanics_model.hh"
#include "material.hh"
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::assembleMassLumped() {
  AKANTU_DEBUG_IN();

  UInt nb_nodes = mesh.getNbNodes();

  if (!mass) {
    std::stringstream sstr_mass; sstr_mass << id << ":mass";
    mass         = &(alloc<Real>(sstr_mass.str(), nb_nodes, spatial_dimension, 0));
  } else
    mass->clear();

  assembleMassLumped(_not_ghost);
  assembleMassLumped(_ghost);

  /// for not connected nodes put mass to one in order to avoid
  /// wrong range in paraview
  Real * mass_values = mass->storage();
  for (UInt i = 0; i < nb_nodes; ++i) {
    if (fabs(mass_values[i]) < std::numeric_limits<Real>::epsilon() || Math::isnan(mass_values[i]))
      mass_values[i] = 1.;
  }

  synch_registry->synchronize(_gst_smm_mass);
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::assembleMassLumped(GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  FEEngine & fem = getFEEngine();

  Array<Real> rho_1(0,1);

  Mesh::type_iterator it  = mesh.firstType(spatial_dimension, ghost_type);
  Mesh::type_iterator end = mesh.lastType(spatial_dimension, ghost_type);
  for(; it != end; ++it) {
    ElementType type = *it;

    computeRho(rho_1, type, ghost_type);

    AKANTU_DEBUG_ASSERT(dof_synchronizer,
			"DOFSynchronizer number must not be initialized");
    fem.assembleFieldLumped(rho_1, spatial_dimension,*mass,
			    dof_synchronizer->getLocalDOFEquationNumbers(),
			    type, ghost_type);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::assembleMass() {
  AKANTU_DEBUG_IN();

  if(!mass_matrix) {
    std::stringstream sstr; sstr << id << ":mass_matrix";
    mass_matrix = new SparseMatrix(*jacobian_matrix, sstr.str(), memory_id);
  }

  assembleMass(_not_ghost);
  //  assembleMass(_ghost);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::assembleMass(GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  MyFEEngineType & fem = getFEEngineClass<MyFEEngineType>();

  Array<Real> rho_1(0,1);
  //UInt nb_element;
  mass_matrix->clear();

  Mesh::type_iterator it  = mesh.firstType(spatial_dimension, ghost_type);
  Mesh::type_iterator end = mesh.lastType(spatial_dimension, ghost_type);
  for(; it != end; ++it) {
    ElementType type = *it;

    computeRho(rho_1, type, ghost_type);
    fem.assembleFieldMatrix(rho_1, spatial_dimension, *mass_matrix, type, ghost_type);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::computeRho(Array<Real> & rho,
				     ElementType type,
				     GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Material ** mat_val = &(materials.at(0));

  FEEngine & fem = getFEEngine();
  UInt nb_element = fem.getMesh().getNbElement(type,ghost_type);

  Array<UInt> & elem_mat_val = element_index_by_material(type, ghost_type);

  UInt nb_quadrature_points = fem.getNbQuadraturePoints(type, ghost_type);

  rho.resize(nb_element * nb_quadrature_points);
  Real * rho_1_val = rho.storage();

  /// compute @f$ rho @f$ for each nodes of each element
  for (UInt el = 0; el < nb_element; ++el) {
    Real mat_rho = mat_val[elem_mat_val(el, 0)]->getParam<Real>("rho"); /// here rho is constant in an element

    for (UInt n = 0; n < nb_quadrature_points; ++n) {
      *rho_1_val++ = mat_rho;
    }
  }

  AKANTU_DEBUG_OUT();
}


__END_AKANTU__

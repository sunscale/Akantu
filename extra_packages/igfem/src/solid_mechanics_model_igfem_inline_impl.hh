/**
 * @file   solid_mechanics_model_igfem_inline_impl.hh
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @date   Wed Nov  4 15:53:52 2015
 *
 * @brief  Implementation on inline functions for SMMIGFEM
 *
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#ifndef __AKANTU_SOLID_MECHANICS_MODEL_IGFEM_INLINE_IMPL_HH__
#define __AKANTU_SOLID_MECHANICS_MODEL_IGFEM_INLINE_IMPL_HH__

namespace akantu {
/* -------------------------------------------------------------------------- */
inline void SolidMechanicsModelIGFEM::getSubElementBarycenter(
    UInt element, UInt sub_element, const ElementType & type,
    Vector<Real> & barycenter, GhostType ghost_type) const {
  UInt * conn_val = this->mesh.getConnectivity(type, ghost_type).storage();
  UInt nb_sub_element_nodes =
      IGFEMHelper::getNbNodesPerSubElement(type, sub_element);
  UInt * sub_el_conn =
      IGFEMHelper::getSubElementConnectivity(type, sub_element);
  UInt nb_nodes_per_element = this->mesh.getNbNodesPerElement(type);
  const Array<Real> & node_coords = this->mesh.getNodes();

  Real local_coord[spatial_dimension * nb_sub_element_nodes];

  UInt offset = element * nb_nodes_per_element;
  for (UInt n = 0; n < nb_sub_element_nodes; ++n) {
    UInt index = conn_val[offset + sub_el_conn[n]];
    memcpy(local_coord + n * spatial_dimension,
           node_coords.storage() + index * spatial_dimension,
           spatial_dimension * sizeof(Real));
  }
  Math::barycenter(local_coord, nb_sub_element_nodes, spatial_dimension,
                   barycenter.storage());
}

} // namespace akantu

#endif /* __AKANTU_SOLID_MECHANICS_MODEL_IGFEM_INLINE_IMPL_HH__ */

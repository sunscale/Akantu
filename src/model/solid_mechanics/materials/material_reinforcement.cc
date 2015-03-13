/**
 * @file   material_reinforcement.cc
 *
 * @author Lucas Frérot <lucas.frerot@epfl.ch>
 *
 * @date creation: Thu Mar 12 2015
 * @date last modification: Thu Mar 12 2015
 *
 * @brief  Reinforcement material
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

#include "material_reinforcement.hh"

__BEGIN_AKANTU__

template<UInt d>
MaterialReinforcement<d>::MaterialReinforcement(SolidMechanicsModel & model, const ID & id):
  Material(model, id),
  model(NULL),
  directing_cosines("directing_cosines", *this),
  reinforcement_stiffness("reinforcement_stiffness", *this),
  area("area", *this),
  shape_derivatives()
{
  this->model = dynamic_cast<EmbeddedInterfaceModel *>(&model);

  const UInt steel_dof = Mesh::getNbNodesPerElement(_segment_2) * d;
  const UInt voigt_size = static_cast<UInt>(d * (d - (d - 1) / 2.));

  directing_cosines.initialize(steel_dof * voigt_size);
}

/* -------------------------------------------------------------------------- */

template<UInt d>
MaterialReinforcement<d>::~MaterialReinforcement()
{}

/* -------------------------------------------------------------------------- */

template<UInt d>
void MaterialReinforcement<d>::initBackgroundShapeDerivatives() {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

template<UInt d>
void MaterialReinforcement<d>::assembleStiffnessMatrix(GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

template<UInt d>
void MaterialReinforcement<d>::computeStiffness(const ElementType & type, GhostType ghost_type) {
  // This methods should be virtual pure
}

/* -------------------------------------------------------------------------- */

template<UInt d>
void MaterialReinforcement<d>::filterInterfaceBackgroundElements(Array<UInt> & filter,
                                                              const ElementType & type,
                                                              const ElementType & interface_type,
                                                              GhostType ghost_type,
                                                              GhostType interface_ghost_type) {
  AKANTU_DEBUG_IN();

  filter.clear();

  Array<Element> & elements = model->getInterfaceAssociatedElements(interface_type, interface_ghost_type);

  for (Array<Element>::scalar_iterator it = elements.begin() ;
      it != elements.end() ;
      ++it) {
    if (it->type == type) filter.push_back(it->element);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

template<UInt d>
void MaterialReinforcement<d>::computeDirectingCosines(const ElementType & type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Mesh & interface_mesh = this->model->getInterfaceMesh();
  
  const UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
  const UInt steel_dof = nb_nodes_per_element * d;
  const UInt voigt_size = static_cast<UInt>(d * (d - (d - 1) / 2.));
  const UInt nb_quad_points = model->getFEEngine().getNbQuadraturePoints(type, ghost_type);

  Array<Real> node_coordinates(this->element_filter(type, ghost_type).getSize(), steel_dof);

  this->model->getFEEngine().template extractNodalToElementField<Real>(interface_mesh,
                                                                       interface_mesh.getNodes(),
                                                                       node_coordinates,
                                                                       type,
                                                                       ghost_type,
                                                                       this->element_filter(type, ghost_type));

  Array<Real>::matrix_iterator
    directing_cosines_it = directing_cosines(type, ghost_type).begin(steel_dof, voigt_size);

  Array<Real>::matrix_iterator node_coordinates_it = node_coordinates.begin(d, nb_nodes_per_element);
  Array<Real>::matrix_iterator node_coordinates_end = node_coordinates.end(d, nb_nodes_per_element);

  for (; node_coordinates_it != node_coordinates_end ; ++node_coordinates_it) {
    for (UInt i = 0 ; i < nb_quad_points ; i++, ++directing_cosines_it) {
      Matrix<Real> & nodes = *node_coordinates_it;
      Matrix<Real> & cosines = *directing_cosines_it;

      computeDirectingCosinesOnQuad(nodes, cosines);
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

template<UInt d>
void MaterialReinforcement<d>::assembleStiffnessMatrix(const ElementType & type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

/// In this function, type and ghost type refer to background elements
template<UInt d>
void MaterialReinforcement<d>::computeBackgroundShapeDerivatives(const ElementType & type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  Mesh & mesh = model->getMesh();
  Mesh & interface_mesh = model->getInterfaceMesh();

  FEEngine & engine = model->getFEEngine();
  FEEngine & interface_engine = model->getFEEngine("EmbeddedInterfaceFEEngine");

  ElementTypeMapArray<Real> quad_pos_map;
  interface_engine.interpolateOnQuadraturePoints(interface_mesh.getNodes(), quad_pos_map);

  Mesh::type_iterator interface_type = interface_mesh.firstType();
  Mesh::type_iterator interface_last = interface_mesh.lastType();

  for (; interface_type != interface_last ; ++interface_type) {
    Array<UInt> filter;

    filterInterfaceBackgroundElements(filter, type, *interface_type, ghost_type, _not_ghost);

    const UInt nb_elements = filter.getSize();
    const UInt nb_nodes = Mesh::getNbNodesPerElement(type);
    const UInt nb_quad_per_element = interface_engine.getNbQuadraturePoints(*interface_type);

    shape_derivatives(*interface_type).alloc(nb_elements * nb_quad_per_element, d * nb_nodes, type, ghost_type);

    Array<Real> & background_shapesd = shape_derivatives(*interface_type)(type, ghost_type);
    Array<Real> & quad_pos = quad_pos_map(*interface_type);

    Array<UInt>::scalar_iterator filter_it = filter.begin();
    Array<UInt>::scalar_iterator filter_end = filter.end();

    Array<Real>::matrix_iterator shapesd_it = background_shapesd.begin(d, nb_nodes);

    Array<Real>::vector_iterator quad_pos_it = quad_pos.begin(d);

    for (; filter_it != filter_end ; ++filter_it) {
      for (UInt i = 0 ; i < nb_quad_per_element ; i++, ++shapesd_it, ++quad_pos_it);
        //engine.computeShapeDerivatives(*quad_pos_it, *filter_it, type, *shapesd_it, ghost_type);
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

INSTANSIATE_MATERIAL(MaterialReinforcement);

/* -------------------------------------------------------------------------- */

__END_AKANTU__

/**
 * @file   boundary_condition_tmpl.hh
 *
 * @author Dana Christen <dana.christen@gmail.com>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri May 03 2013
 * @date last modification: Tue Feb 20 2018
 *
 * @brief  implementation of the applyBC
 *
 *
 * Copyright (©) 2014-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#include "boundary_condition.hh"
#include "element_group.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_BOUNDARY_CONDITION_TMPL_HH__
#define __AKANTU_BOUNDARY_CONDITION_TMPL_HH__

namespace akantu {

/* -------------------------------------------------------------------------- */
template <typename ModelType>
void BoundaryCondition<ModelType>::initBC(ModelType & model,
                                          Array<Real> & primal,
                                          Array<Real> & dual) {
  this->model = &model;
  this->primal = &primal;
  this->dual = &dual;
}

/* -------------------------------------------------------------------------- */
template <typename ModelType>
void BoundaryCondition<ModelType>::initBC(ModelType & model,
                                          Array<Real> & primal,
                                          Array<Real> & primal_increment,
                                          Array<Real> & dual) {
  this->initBC(model, primal, dual);
  this->primal_increment = &primal_increment;
}

/* -------------------------------------------------------------------------- */
/* Partial specialization for DIRICHLET functors */
template <typename ModelType>
template <typename FunctorType>
struct BoundaryCondition<ModelType>::TemplateFunctionWrapper<
    FunctorType, BC::Functor::_dirichlet> {
  static inline void applyBC(const FunctorType & func,
                             const ElementGroup & group,
                             BoundaryCondition<ModelType> & bc_instance) {
    auto & model = bc_instance.getModel();
    auto & primal = bc_instance.getPrimal();

    const auto & coords = model.getMesh().getNodes();
    auto & boundary_flags = model.getBlockedDOFs();
    UInt dim = model.getMesh().getSpatialDimension();

    auto primal_iter = primal.begin(primal.getNbComponent());
    auto coords_iter = coords.begin(dim);
    auto flags_iter = boundary_flags.begin(boundary_flags.getNbComponent());

    for (auto n : group.getNodeGroup()) {
      Vector<bool> flag(flags_iter[n]);
      Vector<Real> primal(primal_iter[n]);
      Vector<Real> coords(coords_iter[n]);
      func(n, flag, primal, coords);
    }
  }
};

/* -------------------------------------------------------------------------- */
/* Partial specialization for NEUMANN functors */
template <typename ModelType>
template <typename FunctorType>
struct BoundaryCondition<ModelType>::TemplateFunctionWrapper<
    FunctorType, BC::Functor::_neumann> {
  static inline void applyBC(const FunctorType & func,
                             const ElementGroup & group,
                             BoundaryCondition<ModelType> & bc_instance) {
    UInt dim = bc_instance.getModel().getSpatialDimension();
    switch (dim) {
    case 1: {
      AKANTU_TO_IMPLEMENT();
      break;
    }
    case 2:
    case 3: {
      applyBC(func, group, bc_instance, _not_ghost);
      applyBC(func, group, bc_instance, _ghost);
      break;
    }
    }
  }

  static inline void applyBC(const FunctorType & func,
                             const ElementGroup & group,
                             BoundaryCondition<ModelType> & bc_instance,
                             GhostType ghost_type) {
    auto & model = bc_instance.getModel();
    auto & dual = bc_instance.getDual();
    const auto & mesh = model.getMesh();
    const auto & nodes_coords = mesh.getNodes();
    const auto & fem_boundary = model.getFEEngineBoundary();

    UInt dim = model.getSpatialDimension();
    UInt nb_degree_of_freedom = dual.getNbComponent();

    IntegrationPoint quad_point;
    quad_point.ghost_type = ghost_type;

    // Loop over the boundary element types
    for (auto && type : group.elementTypes(dim - 1, ghost_type)) {
      const auto & element_ids = group.getElements(type, ghost_type);

      UInt nb_quad_points =
          fem_boundary.getNbIntegrationPoints(type, ghost_type);
      UInt nb_elements = element_ids.size();
      UInt nb_nodes_per_element = mesh.getNbNodesPerElement(type);

      Array<Real> dual_before_integ(nb_elements * nb_quad_points,
                                    nb_degree_of_freedom, 0.);
      Array<Real> quad_coords(nb_elements * nb_quad_points, dim);

      const auto & normals_on_quad =
          fem_boundary.getNormalsOnIntegrationPoints(type, ghost_type);

      fem_boundary.interpolateOnIntegrationPoints(
          nodes_coords, quad_coords, dim, type, ghost_type, element_ids);
      auto normals_begin = normals_on_quad.begin(dim);
      decltype(normals_begin) normals_iter;
      auto quad_coords_iter = quad_coords.begin(dim);
      auto dual_iter = dual_before_integ.begin(nb_degree_of_freedom);

      quad_point.type = type;
      for (auto el : element_ids) {
        quad_point.element = el;
        normals_iter = normals_begin + el * nb_quad_points;
        for (auto q : arange(nb_quad_points)) {
          quad_point.num_point = q;
          func(quad_point, *dual_iter, *quad_coords_iter, *normals_iter);
          ++dual_iter;
          ++quad_coords_iter;
          ++normals_iter;
        }
      }

      Array<Real> dual_by_shapes(nb_elements * nb_quad_points,
                                 nb_degree_of_freedom * nb_nodes_per_element);

      fem_boundary.computeNtb(dual_before_integ, dual_by_shapes, type,
                              ghost_type, element_ids);

      Array<Real> dual_by_shapes_integ(nb_elements, nb_degree_of_freedom *
                                                        nb_nodes_per_element);
      fem_boundary.integrate(dual_by_shapes, dual_by_shapes_integ,
                             nb_degree_of_freedom * nb_nodes_per_element, type,
                             ghost_type, element_ids);

      // assemble the result into force vector
      model.getDOFManager().assembleElementalArrayLocalArray(
          dual_by_shapes_integ, dual, type, ghost_type, 1., element_ids);
    }
  }
};

/* -------------------------------------------------------------------------- */
template <typename ModelType>
template <typename FunctorType>
inline void BoundaryCondition<ModelType>::applyBC(const FunctorType & func) {
  auto bit = model->getMesh().getGroupManager().element_group_begin();
  auto bend = model->getMesh().getGroupManager().element_group_end();
  for (; bit != bend; ++bit)
    applyBC(func, *bit);
}

/* -------------------------------------------------------------------------- */
template <typename ModelType>
template <typename FunctorType>
inline void
BoundaryCondition<ModelType>::applyBC(const FunctorType & func,
                                      const std::string & group_name) {
  try {
    const ElementGroup & element_group =
        model->getMesh().getElementGroup(group_name);
    applyBC(func, element_group);
  } catch (akantu::debug::Exception & e) {
    AKANTU_EXCEPTION("Error applying a boundary condition onto \""
                     << group_name << "\"! [" << e.what() << "]");
  }
}

/* -------------------------------------------------------------------------- */
template <typename ModelType>
template <typename FunctorType>
inline void
BoundaryCondition<ModelType>::applyBC(const FunctorType & func,
                                      const ElementGroup & element_group) {
#if !defined(AKANTU_NDEBUG)
  if (element_group.getDimension() != model->getSpatialDimension() - 1)
    AKANTU_DEBUG_WARNING("The group "
                         << element_group.getName()
                         << " does not contain only boundaries elements");
#endif

  TemplateFunctionWrapper<FunctorType>::applyBC(func, element_group, *this);
}

#endif /* __AKANTU_BOUNDARY_CONDITION_TMPL_HH__ */

} // namespace akantu

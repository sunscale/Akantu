/**
 * @file   boundary_condition_tmpl.hh
 *
 * @author Dana Christen <dana.christen@gmail.com>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri May 03 2013
 * @date last modification: Mon Jun 23 2014
 *
 * @brief  XXX
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
#include "element_group.hh"
/* -------------------------------------------------------------------------- */


__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template<typename ModelType>
void BoundaryCondition<ModelType>::initBC(ModelType & model, Array<Real> & primal, Array<Real> & dual)
{
  this->model  = &model;
  this->primal = &primal;
  this->dual   = &dual;
  if(this->model->getSpatialDimension() > 1)
    this->model->initFEEngineBoundary();
}

/* -------------------------------------------------------------------------- */
template<typename ModelType>
void BoundaryCondition<ModelType>::initBC(ModelType & model,
					  Array<Real> & primal,
					  Array<Real> & primal_increment,
					  Array<Real> & dual)
{
  this->initBC(model, primal, dual);
  this->primal_increment = &primal_increment;
}
/* -------------------------------------------------------------------------- */
/* Partial specialization for DIRICHLET functors */
template<typename ModelType>
template<typename FunctorType>
struct BoundaryCondition<ModelType>::TemplateFunctionWrapper<FunctorType, BC::Functor::_dirichlet> {
  static inline void applyBC(const FunctorType & func,
                             const ElementGroup & group,
                             BoundaryCondition<ModelType> & bc_instance) {
    ModelType         & model          = bc_instance.getModel();
    Array<Real>       & primal         = bc_instance.getPrimal();
    const Array<Real> & coords         = model.getMesh().getNodes();
    Array<bool>       & boundary_flags = model.getBlockedDOFs();
    UInt dim = model.getMesh().getSpatialDimension();

    Array<Real>::vector_iterator primal_iter = primal.begin(primal.getNbComponent());
    Array<Real>::const_vector_iterator coords_iter = coords.begin(dim);
    Array<bool>::vector_iterator flags_iter = boundary_flags.begin(boundary_flags.getNbComponent());

    for(ElementGroup::const_node_iterator nodes_it(group.node_begin());
        nodes_it!= group.node_end();
        ++nodes_it) {
      UInt n = *nodes_it;
      func(n, flags_iter[n], primal_iter[n], coords_iter[n]);
    }
  }
};

/* -------------------------------------------------------------------------- */
/* Partial specialization for NEUMANN functors */
template<typename ModelType>
template<typename FunctorType>
struct BoundaryCondition<ModelType>::TemplateFunctionWrapper<FunctorType, BC::Functor::_neumann> {
  static inline void applyBC(const FunctorType & func,
                             const ElementGroup & group,
                             BoundaryCondition<ModelType> & bc_instance) {
    UInt dim = bc_instance.getModel().getSpatialDimension();
    switch(dim) {
    case 1: { AKANTU_DEBUG_TO_IMPLEMENT(); break; }
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
    ModelType &         model        = bc_instance.getModel();
    Array<Real>       & dual         = bc_instance.getDual();
    const Mesh        & mesh         = model.getMesh();
    const Array<Real> & nodes_coords = mesh.getNodes();
    const FEEngine         & fem_boundary = model.getFEEngineBoundary();

    UInt dim = model.getSpatialDimension();
    UInt nb_degree_of_freedom = dual.getNbComponent();

    QuadraturePoint quad_point;
    quad_point.ghost_type = ghost_type;

    ElementGroup::type_iterator type_it  = group.firstType(dim - 1, ghost_type);
    ElementGroup::type_iterator type_end = group.lastType (dim - 1, ghost_type);

    // Loop over the boundary element types
    for(; type_it != type_end; ++type_it) {
      const Array<UInt> & element_ids = group.getElements(*type_it, ghost_type);

      Array<UInt>::const_scalar_iterator elem_iter = element_ids.begin();
      Array<UInt>::const_scalar_iterator elem_iter_end = element_ids.end();

      UInt nb_quad_points = fem_boundary.getNbQuadraturePoints(*type_it, ghost_type);
      UInt nb_elements = element_ids.getSize();
      UInt nb_nodes_per_element = mesh.getNbNodesPerElement(*type_it);

      Array<Real> * dual_before_integ = new Array<Real>(nb_elements * nb_quad_points,
							nb_degree_of_freedom,
							0.);
      Array<Real> * quad_coords = new Array<Real>(nb_elements * nb_quad_points, dim);

      const Array<Real> & normals_on_quad =
	fem_boundary.getNormalsOnQuadPoints(*type_it, ghost_type);

      fem_boundary.interpolateOnQuadraturePoints(nodes_coords, *quad_coords,
						 dim,
						 *type_it, ghost_type,
						 element_ids);
      Array<Real>::const_vector_iterator normals_begin = normals_on_quad.begin(dim);
      Array<Real>::const_vector_iterator normals_iter;
      Array<Real>::const_vector_iterator quad_coords_iter  = quad_coords->begin(dim);
      Array<Real>::vector_iterator dual_iter = dual_before_integ->begin(nb_degree_of_freedom);

      quad_point.type = *type_it;
      for(; elem_iter != elem_iter_end; ++elem_iter) {
        UInt el = *elem_iter;
        quad_point.element = el;
	normals_iter = normals_begin + el * nb_quad_points;
        for(UInt q(0); q < nb_quad_points; ++q) {
          quad_point.num_point = q;
          func(quad_point,
               *dual_iter,
               *quad_coords_iter,
               *normals_iter);
          ++dual_iter;
          ++quad_coords_iter;
          ++normals_iter;
        }
      }

      delete quad_coords;

      /* -------------------------------------------------------------------- */
      // Initialization of iterators
      Array<Real>::matrix_iterator dual_iter_mat =
	dual_before_integ->begin(nb_degree_of_freedom,1);
      elem_iter = element_ids.begin();
      Array<Real>::const_matrix_iterator shapes_iter_begin =
	fem_boundary.getShapes(*type_it, ghost_type).begin(1, nb_nodes_per_element);

      Array<Real> * dual_by_shapes =
	new Array<Real>(nb_elements*nb_quad_points, nb_degree_of_freedom*nb_nodes_per_element);

      Array<Real>::matrix_iterator dual_by_shapes_iter =
	dual_by_shapes->begin(nb_degree_of_freedom, nb_nodes_per_element);

      Array<Real>::const_matrix_iterator shapes_iter;

      /* -------------------------------------------------------------------- */
      // Loop computing dual x shapes
      for(; elem_iter != elem_iter_end; ++elem_iter) {
	shapes_iter = shapes_iter_begin + *elem_iter*nb_quad_points;

        for(UInt q(0); q < nb_quad_points; ++q,
	      ++dual_iter_mat, ++dual_by_shapes_iter, ++shapes_iter) {
          dual_by_shapes_iter->mul<false, false>(*dual_iter_mat, *shapes_iter);
        }
      }

      delete dual_before_integ;

      Array<Real> * dual_by_shapes_integ =
        new Array<Real>(nb_elements, nb_degree_of_freedom*nb_nodes_per_element);
      fem_boundary.integrate(*dual_by_shapes,
                             *dual_by_shapes_integ,
                             nb_degree_of_freedom*nb_nodes_per_element,
                             *type_it,
                             ghost_type,
                             element_ids);
      delete dual_by_shapes;

      // assemble the result into force vector
      fem_boundary.assembleArray(*dual_by_shapes_integ,
                                 dual,
                                 model.getDOFSynchronizer().getLocalDOFEquationNumbers(),
                                 nb_degree_of_freedom,
                                 *type_it,
                                 ghost_type,
                                 element_ids);
      delete dual_by_shapes_integ;
    }
  }
};

/* -------------------------------------------------------------------------- */
template<typename ModelType>
template<typename FunctorType>
inline void BoundaryCondition<ModelType>::applyBC(const FunctorType & func) {
  GroupManager::const_element_group_iterator bit = model->getMesh().getGroupManager().element_group_begin();
  GroupManager::const_element_group_iterator bend = model->getMesh().getGroupManager().element_group_end();
  for(; bit != bend; ++bit) applyBC(func, *bit);
}

/* -------------------------------------------------------------------------- */
template<typename ModelType>
template<typename FunctorType>
inline void BoundaryCondition<ModelType>::applyBC(const FunctorType & func,
						  const std::string & group_name) {
  try {
    const ElementGroup & element_group = model->getMesh().getElementGroup(group_name);
    applyBC(func, element_group);
  } catch(akantu::debug::Exception e) {
    AKANTU_EXCEPTION("Error applying a boundary condition onto \""
		     << group_name << "\"! [" << e.what() <<"]");
  }
}

/* -------------------------------------------------------------------------- */
template<typename ModelType>
template<typename FunctorType>
inline void BoundaryCondition<ModelType>::applyBC(const FunctorType & func,
						  const ElementGroup & element_group) {
#if !defined(AKANTU_NDEBUG)
  if(element_group.getDimension() != model->getSpatialDimension() - 1)
    AKANTU_DEBUG_WARNING("The group " << element_group.getName()
                         << " does not contain only boundaries elements");
#endif

  TemplateFunctionWrapper<FunctorType>::applyBC(func, element_group, *this);
}

__END_AKANTU__

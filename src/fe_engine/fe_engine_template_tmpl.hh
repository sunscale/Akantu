/**
 * @file   fe_engine_template_tmpl.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Mauro Corrado <mauro.corrado@epfl.ch>
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Tue Feb 15 2011
 * @date last modification: Thu Nov 19 2015
 *
 * @brief  Template implementation of FEEngineTemplate
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2012, 2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de
 * Lausanne)  Laboratory (LSMS  -  Laboratoire de  Simulation  en Mécanique  des
 * Solides)
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
#include "aka_common.hh"
#include "dof_manager.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::FEEngineTemplate(
    Mesh & mesh, UInt spatial_dimension, ID id, MemoryID memory_id)
    : FEEngine(mesh, spatial_dimension, id, memory_id),
      integrator(mesh, id, memory_id), shape_functions(mesh, id, memory_id) {}

/* -------------------------------------------------------------------------- */
template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::~FEEngineTemplate() {}

/* -------------------------------------------------------------------------- */
/**
 * Helper class to be able to write a partial specialization on the element kind
 */
template <ElementKind kind> struct GradientOnIntegrationPointsHelper {
  template <class S>
  static void call(__attribute__((unused)) const S & shape_functions,
                   __attribute__((unused)) Mesh & mesh,
                   __attribute__((unused)) const Array<Real> & u,
                   __attribute__((unused)) Array<Real> & nablauq,
                   __attribute__((unused)) const UInt nb_degree_of_freedom,
                   __attribute__((unused)) const ElementType & type,
                   __attribute__((unused)) const GhostType & ghost_type,
                   __attribute__((unused))
                   const Array<UInt> & filter_elements) {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }
};

#define COMPUTE_GRADIENT(type)                                                 \
  if (element_dimension == ElementClass<type>::getSpatialDimension())          \
    shape_functions.template gradientOnIntegrationPoints<type>(                \
        u, nablauq, nb_degree_of_freedom, ghost_type, filter_elements);

#define AKANTU_SPECIALIZE_GRADIENT_ON_INTEGRATION_POINTS_HELPER(kind)          \
  template <> struct GradientOnIntegrationPointsHelper<kind> {                 \
    template <class S>                                                         \
    static void call(const S & shape_functions, Mesh & mesh,                   \
                     const Array<Real> & u, Array<Real> & nablauq,             \
                     const UInt nb_degree_of_freedom,                          \
                     const ElementType & type, const GhostType & ghost_type,   \
                     const Array<UInt> & filter_elements) {                    \
      UInt element_dimension = mesh.getSpatialDimension(type);                 \
      AKANTU_BOOST_KIND_ELEMENT_SWITCH(COMPUTE_GRADIENT, kind);                \
    }                                                                          \
  };
AKANTU_BOOST_ALL_KIND_LIST(
    AKANTU_SPECIALIZE_GRADIENT_ON_INTEGRATION_POINTS_HELPER,
    AKANTU_FE_ENGINE_LIST_GRADIENT_ON_INTEGRATION_POINTS)

#undef AKANTU_SPECIALIZE_GRADIENT_ON_INTEGRATION_POINTS_HELPER
#undef COMPUTE_GRADIENT

template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
void FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::
    gradientOnIntegrationPoints(const Array<Real> & u, Array<Real> & nablauq,
                                const UInt nb_degree_of_freedom,
                                const ElementType & type,
                                const GhostType & ghost_type,
                                const Array<UInt> & filter_elements) const {
  AKANTU_DEBUG_IN();

  UInt nb_element = mesh.getNbElement(type, ghost_type);
  if (filter_elements != empty_filter)
    nb_element = filter_elements.getSize();
  UInt nb_points =
      shape_functions.getIntegrationPoints(type, ghost_type).cols();

#ifndef AKANTU_NDEBUG

  UInt element_dimension = mesh.getSpatialDimension(type);

  AKANTU_DEBUG_ASSERT(u.getSize() == mesh.getNbNodes(),
                      "The vector u(" << u.getID()
                                      << ") has not the good size.");
  AKANTU_DEBUG_ASSERT(u.getNbComponent() == nb_degree_of_freedom,
                      "The vector u("
                          << u.getID()
                          << ") has not the good number of component.");

  AKANTU_DEBUG_ASSERT(
      nablauq.getNbComponent() == nb_degree_of_freedom * element_dimension,
      "The vector nablauq(" << nablauq.getID()
                            << ") has not the good number of component.");

// AKANTU_DEBUG_ASSERT(nablauq.getSize() == nb_element * nb_points,
//                  "The vector nablauq(" << nablauq.getID()
//                  << ") has not the good size.");
#endif

  nablauq.resize(nb_element * nb_points);

  GradientOnIntegrationPointsHelper<kind>::call(
      shape_functions, mesh, u, nablauq, nb_degree_of_freedom, type, ghost_type,
      filter_elements);
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
void FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::initShapeFunctions(
    const GhostType & ghost_type) {
  initShapeFunctions(mesh.getNodes(), ghost_type);
}

/* -------------------------------------------------------------------------- */
template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
void FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::initShapeFunctions(
    const Array<Real> & nodes, const GhostType & ghost_type) {
  AKANTU_DEBUG_IN();

  Mesh::type_iterator it = mesh.firstType(element_dimension, ghost_type, kind);
  Mesh::type_iterator end = mesh.lastType(element_dimension, ghost_type, kind);
  for (; it != end; ++it) {
    ElementType type = *it;
    integrator.initIntegrator(nodes, type, ghost_type);
    const Matrix<Real> & control_points =
        getIntegrationPoints(type, ghost_type);
    shape_functions.initShapeFunctions(nodes, control_points, type, ghost_type);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/**
 * Helper class to be able to write a partial specialization on the element kind
 */
template <ElementKind kind> struct IntegrateHelper {};

#define INTEGRATE(type)                                                        \
  integrator.template integrate<type>(f, intf, nb_degree_of_freedom,           \
                                      ghost_type, filter_elements);

#define AKANTU_SPECIALIZE_INTEGRATE_HELPER(kind)                               \
  template <> struct IntegrateHelper<kind> {                                   \
    template <class I>                                                         \
    static void call(const I & integrator, const Array<Real> & f,              \
                     Array<Real> & intf, UInt nb_degree_of_freedom,            \
                     const ElementType & type, const GhostType & ghost_type,   \
                     const Array<UInt> & filter_elements) {                    \
      AKANTU_BOOST_KIND_ELEMENT_SWITCH(INTEGRATE, kind);                       \
    }                                                                          \
  };

AKANTU_BOOST_ALL_KIND(AKANTU_SPECIALIZE_INTEGRATE_HELPER)

#undef AKANTU_SPECIALIZE_INTEGRATE_HELPER
#undef INTEGRATE

template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
void FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::integrate(
    const Array<Real> & f, Array<Real> & intf, UInt nb_degree_of_freedom,
    const ElementType & type, const GhostType & ghost_type,
    const Array<UInt> & filter_elements) const {

  UInt nb_element = mesh.getNbElement(type, ghost_type);
  if (filter_elements != empty_filter)
    nb_element = filter_elements.getSize();
#ifndef AKANTU_NDEBUG

  UInt nb_quadrature_points = getNbIntegrationPoints(type);

  AKANTU_DEBUG_ASSERT(f.getSize() == nb_element * nb_quadrature_points,
                      "The vector f(" << f.getID() << " size " << f.getSize()
                                      << ") has not the good size ("
                                      << nb_element << ").");
  AKANTU_DEBUG_ASSERT(f.getNbComponent() == nb_degree_of_freedom,
                      "The vector f("
                          << f.getID()
                          << ") has not the good number of component.");
  AKANTU_DEBUG_ASSERT(intf.getNbComponent() == nb_degree_of_freedom,
                      "The vector intf("
                          << intf.getID()
                          << ") has not the good number of component.");
#endif

  intf.resize(nb_element);

  IntegrateHelper<kind>::call(integrator, f, intf, nb_degree_of_freedom, type,
                              ghost_type, filter_elements);
}

/* -------------------------------------------------------------------------- */
/**
 * Helper class to be able to write a partial specialization on the element kind
 */
template <ElementKind kind> struct IntegrateScalarHelper {};

#define INTEGRATE(type)                                                        \
  integral =                                                                   \
      integrator.template integrate<type>(f, ghost_type, filter_elements);

#define AKANTU_SPECIALIZE_INTEGRATE_SCALAR_HELPER(kind)                        \
  template <> struct IntegrateScalarHelper<kind> {                             \
    template <class I>                                                         \
    static Real call(const I & integrator, const Array<Real> & f,              \
                     const ElementType & type, const GhostType & ghost_type,   \
                     const Array<UInt> & filter_elements) {                    \
      Real integral = 0.;                                                      \
      AKANTU_BOOST_KIND_ELEMENT_SWITCH(INTEGRATE, kind);                       \
      return integral;                                                         \
    }                                                                          \
  };

AKANTU_BOOST_ALL_KIND(AKANTU_SPECIALIZE_INTEGRATE_SCALAR_HELPER)

#undef AKANTU_SPECIALIZE_INTEGRATE_SCALAR_HELPER
#undef INTEGRATE

template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
Real FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::integrate(
    const Array<Real> & f, const ElementType & type,
    const GhostType & ghost_type, const Array<UInt> & filter_elements) const {
  AKANTU_DEBUG_IN();

#ifndef AKANTU_NDEBUG
  //   std::stringstream sstr; sstr << ghost_type;
  //   AKANTU_DEBUG_ASSERT(sstr.str() == nablauq.getTag(),
  //                  "The vector " << nablauq.getID() << " is not taged " <<
  //                  ghost_type);
  UInt nb_element = mesh.getNbElement(type, ghost_type);
  if (filter_elements != empty_filter)
    nb_element = filter_elements.getSize();

  UInt nb_quadrature_points = getNbIntegrationPoints(type, ghost_type);

  AKANTU_DEBUG_ASSERT(f.getSize() == nb_element * nb_quadrature_points,
                      "The vector f("
                          << f.getID() << ") has not the good size. ("
                          << f.getSize()
                          << "!=" << nb_quadrature_points * nb_element << ")");
  AKANTU_DEBUG_ASSERT(f.getNbComponent() == 1,
                      "The vector f("
                          << f.getID()
                          << ") has not the good number of component.");
#endif

  Real integral = IntegrateScalarHelper<kind>::call(
      integrator, f, type, ghost_type, filter_elements);
  AKANTU_DEBUG_OUT();
  return integral;
}

/* -------------------------------------------------------------------------- */
/**
 * Helper class to be able to write a partial specialization on the element kind
 */
template <ElementKind kind> struct IntegrateScalarOnOneElementHelper {};

#define INTEGRATE(type)                                                        \
  res = integrator.template integrate<type>(f, index, ghost_type);

#define AKANTU_SPECIALIZE_INTEGRATE_SCALAR_ON_ONE_ELEMENT_HELPER(kind)         \
  template <> struct IntegrateScalarOnOneElementHelper<kind> {                 \
    template <class I>                                                         \
    static Real call(const I & integrator, const Vector<Real> & f,             \
                     const ElementType & type, UInt index,                     \
                     const GhostType & ghost_type) {                           \
      Real res = 0.;                                                           \
      AKANTU_BOOST_KIND_ELEMENT_SWITCH(INTEGRATE, kind);                       \
      return res;                                                              \
    }                                                                          \
  };

AKANTU_BOOST_ALL_KIND(AKANTU_SPECIALIZE_INTEGRATE_SCALAR_ON_ONE_ELEMENT_HELPER)

#undef AKANTU_SPECIALIZE_INTEGRATE_SCALAR_ON_ONE_ELEMENT_HELPER
#undef INTEGRATE

template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
Real FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::integrate(
    const Vector<Real> & f, const ElementType & type, UInt index,
    const GhostType & ghost_type) const {

  Real res = IntegrateScalarOnOneElementHelper<kind>::call(integrator, f, type,
                                                           index, ghost_type);
  return res;
}

/* -------------------------------------------------------------------------- */
/**
 * Helper class to be able to write a partial specialization on the element kind
 */
template <ElementKind kind> struct IntegrateOnIntegrationPointsHelper {};

#define INTEGRATE(type)                                                        \
  integrator.template integrateOnIntegrationPoints<type>(                      \
      f, intf, nb_degree_of_freedom, ghost_type, filter_elements);

#define AKANTU_SPECIALIZE_INTEGRATE_ON_INTEGRATION_POINTS_HELPER(kind)         \
  template <> struct IntegrateOnIntegrationPointsHelper<kind> {                \
    template <class I>                                                         \
    static void call(const I & integrator, const Array<Real> & f,              \
                     Array<Real> & intf, UInt nb_degree_of_freedom,            \
                     const ElementType & type, const GhostType & ghost_type,   \
                     const Array<UInt> & filter_elements) {                    \
      AKANTU_BOOST_KIND_ELEMENT_SWITCH(INTEGRATE, kind);                       \
    }                                                                          \
  };

AKANTU_BOOST_ALL_KIND(AKANTU_SPECIALIZE_INTEGRATE_ON_INTEGRATION_POINTS_HELPER)

#undef AKANTU_SPECIALIZE_INTEGRATE_ON_INTEGRATION_POINTS_HELPER
#undef INTEGRATE

template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
void FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::
    integrateOnIntegrationPoints(const Array<Real> & f, Array<Real> & intf,
                                 UInt nb_degree_of_freedom,
                                 const ElementType & type,
                                 const GhostType & ghost_type,
                                 const Array<UInt> & filter_elements) const {

  UInt nb_element = mesh.getNbElement(type, ghost_type);
  if (filter_elements != empty_filter)
    nb_element = filter_elements.getSize();
  UInt nb_quadrature_points = getNbIntegrationPoints(type);
#ifndef AKANTU_NDEBUG
  //   std::stringstream sstr; sstr << ghost_type;
  //   AKANTU_DEBUG_ASSERT(sstr.str() == nablauq.getTag(),
  //                  "The vector " << nablauq.getID() << " is not taged " <<
  //                  ghost_type);

  AKANTU_DEBUG_ASSERT(f.getSize() == nb_element * nb_quadrature_points,
                      "The vector f(" << f.getID() << " size " << f.getSize()
                                      << ") has not the good size ("
                                      << nb_element << ").");
  AKANTU_DEBUG_ASSERT(f.getNbComponent() == nb_degree_of_freedom,
                      "The vector f("
                          << f.getID()
                          << ") has not the good number of component.");
  AKANTU_DEBUG_ASSERT(intf.getNbComponent() == nb_degree_of_freedom,
                      "The vector intf("
                          << intf.getID()
                          << ") has not the good number of component.");
#endif

  intf.resize(nb_element * nb_quadrature_points);
  IntegrateOnIntegrationPointsHelper<kind>::call(integrator, f, intf,
                                                 nb_degree_of_freedom, type,
                                                 ghost_type, filter_elements);
}

/* -------------------------------------------------------------------------- */
/**
 * Helper class to be able to write a partial specialization on the element kind
 */
template <ElementKind kind> struct InterpolateOnIntegrationPointsHelper {
  template <class S>
  static void call(__attribute__((unused)) const S & shape_functions,
                   __attribute__((unused)) const Array<Real> & u,
                   __attribute__((unused)) Array<Real> & uq,
                   __attribute__((unused)) const UInt nb_degree_of_freedom,
                   __attribute__((unused)) const ElementType & type,
                   __attribute__((unused)) const GhostType & ghost_type,
                   __attribute__((unused))
                   const Array<UInt> & filter_elements) {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }
};

#define INTERPOLATE(type)                                                      \
  shape_functions.template interpolateOnIntegrationPoints<type>(               \
      u, uq, nb_degree_of_freedom, ghost_type, filter_elements);

#define AKANTU_SPECIALIZE_INTERPOLATE_ON_INTEGRATION_POINTS_HELPER(kind)       \
  template <> struct InterpolateOnIntegrationPointsHelper<kind> {              \
    template <class S>                                                         \
    static void call(const S & shape_functions, const Array<Real> & u,         \
                     Array<Real> & uq, const UInt nb_degree_of_freedom,        \
                     const ElementType & type, const GhostType & ghost_type,   \
                     const Array<UInt> & filter_elements) {                    \
      AKANTU_BOOST_KIND_ELEMENT_SWITCH(INTERPOLATE, kind);                     \
    }                                                                          \
  };
AKANTU_BOOST_ALL_KIND_LIST(
    AKANTU_SPECIALIZE_INTERPOLATE_ON_INTEGRATION_POINTS_HELPER,
    AKANTU_FE_ENGINE_LIST_INTERPOLATE_ON_INTEGRATION_POINTS)

#undef AKANTU_SPECIALIZE_INTERPOLATE_ON_INTEGRATION_POINTS_HELPER
#undef INTERPOLATE

template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
void FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::
    interpolateOnIntegrationPoints(const Array<Real> & u, Array<Real> & uq,
                                   const UInt nb_degree_of_freedom,
                                   const ElementType & type,
                                   const GhostType & ghost_type,
                                   const Array<UInt> & filter_elements) const {
  AKANTU_DEBUG_IN();

  UInt nb_points =
      shape_functions.getIntegrationPoints(type, ghost_type).cols();
  UInt nb_element = mesh.getNbElement(type, ghost_type);
  if (filter_elements != empty_filter)
    nb_element = filter_elements.getSize();
#ifndef AKANTU_NDEBUG

  AKANTU_DEBUG_ASSERT(u.getSize() == mesh.getNbNodes(),
                      "The vector u(" << u.getID()
                                      << ") has not the good size.");
  AKANTU_DEBUG_ASSERT(u.getNbComponent() == nb_degree_of_freedom,
                      "The vector u("
                          << u.getID()
                          << ") has not the good number of component.");
  AKANTU_DEBUG_ASSERT(uq.getNbComponent() == nb_degree_of_freedom,
                      "The vector uq("
                          << uq.getID()
                          << ") has not the good number of component.");
#endif

  uq.resize(nb_element * nb_points);

  InterpolateOnIntegrationPointsHelper<kind>::call(shape_functions, u, uq,
                                                   nb_degree_of_freedom, type,
                                                   ghost_type, filter_elements);
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
void FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::
    interpolateOnIntegrationPoints(
        const Array<Real> & u, ElementTypeMapArray<Real> & uq,
        const ElementTypeMapArray<UInt> * filter_elements) const {
  AKANTU_DEBUG_IN();

  const Array<UInt> * filter = NULL;

  for (ghost_type_t::iterator gt = ghost_type_t::begin();
       gt != ghost_type_t::end(); ++gt) {

    GhostType ghost_type = *gt;
    ElementTypeMapArray<Real>::type_iterator it =
        uq.firstType(_all_dimensions, ghost_type, kind);
    ElementTypeMapArray<Real>::type_iterator last =
        uq.lastType(_all_dimensions, ghost_type, kind);

    for (; it != last; ++it) {
      ElementType type = *it;

      UInt nb_quad_per_element = getNbIntegrationPoints(type, ghost_type);

      UInt nb_element = 0;

      if (filter_elements) {
        filter = &((*filter_elements)(type, ghost_type));
        nb_element = filter->getSize();
      } else {
        filter = &empty_filter;
        nb_element = mesh.getNbElement(type, ghost_type);
      }

      UInt nb_tot_quad = nb_quad_per_element * nb_element;

      Array<Real> & quad = uq(type, ghost_type);
      quad.resize(nb_tot_quad);

      interpolateOnIntegrationPoints(u, quad, quad.getNbComponent(), type,
                                     ghost_type, *filter);
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
inline void FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::
    computeIntegrationPointsCoordinates(
        ElementTypeMapArray<Real> & quadrature_points_coordinates,
        const ElementTypeMapArray<UInt> * filter_elements) const {

  const Array<Real> & nodes_coordinates = mesh.getNodes();

  interpolateOnIntegrationPoints(
      nodes_coordinates, quadrature_points_coordinates, filter_elements);
}

/* -------------------------------------------------------------------------- */
template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
inline void FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::
    computeIntegrationPointsCoordinates(
        Array<Real> & quadrature_points_coordinates, const ElementType & type,
        const GhostType & ghost_type,
        const Array<UInt> & filter_elements) const {

  const Array<Real> & nodes_coordinates = mesh.getNodes();

  UInt spatial_dimension = mesh.getSpatialDimension();

  interpolateOnIntegrationPoints(
      nodes_coordinates, quadrature_points_coordinates, spatial_dimension, type,
      ghost_type, filter_elements);
}

/* -------------------------------------------------------------------------- */
template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
inline void FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::
    initElementalFieldInterpolationFromIntegrationPoints(
        const ElementTypeMapArray<Real> & interpolation_points_coordinates,
        ElementTypeMapArray<Real> & interpolation_points_coordinates_matrices,
        ElementTypeMapArray<Real> & quad_points_coordinates_inv_matrices,
        const ElementTypeMapArray<UInt> * element_filter) const {

  AKANTU_DEBUG_IN();

  UInt spatial_dimension = this->mesh.getSpatialDimension();

  ElementTypeMapArray<Real> quadrature_points_coordinates(
      "quadrature_points_coordinates_for_interpolation", getID(),
      getMemoryID());

  quadrature_points_coordinates.initialize(
      *this, _nb_component = spatial_dimension,
      _spatial_dimension = spatial_dimension);

  computeIntegrationPointsCoordinates(quadrature_points_coordinates,
                                      element_filter);
  shape_functions.initElementalFieldInterpolationFromIntegrationPoints(
      interpolation_points_coordinates,
      interpolation_points_coordinates_matrices,
      quad_points_coordinates_inv_matrices, quadrature_points_coordinates,
      element_filter);
}

/* -------------------------------------------------------------------------- */
template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
inline void FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::
    interpolateElementalFieldFromIntegrationPoints(
        const ElementTypeMapArray<Real> & field,
        const ElementTypeMapArray<Real> & interpolation_points_coordinates,
        ElementTypeMapArray<Real> & result, const GhostType ghost_type,
        const ElementTypeMapArray<UInt> * element_filter) const {

  ElementTypeMapArray<Real> interpolation_points_coordinates_matrices(
      "interpolation_points_coordinates_matrices", id, memory_id);
  ElementTypeMapArray<Real> quad_points_coordinates_inv_matrices(
      "quad_points_coordinates_inv_matrices", id, memory_id);

  initElementalFieldInterpolationFromIntegrationPoints(
      interpolation_points_coordinates,
      interpolation_points_coordinates_matrices,
      quad_points_coordinates_inv_matrices, element_filter);

  interpolateElementalFieldFromIntegrationPoints(
      field, interpolation_points_coordinates_matrices,
      quad_points_coordinates_inv_matrices, result, ghost_type, element_filter);
}

/* -------------------------------------------------------------------------- */
template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
inline void FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::
    interpolateElementalFieldFromIntegrationPoints(
        const ElementTypeMapArray<Real> & field,
        const ElementTypeMapArray<Real> &
            interpolation_points_coordinates_matrices,
        const ElementTypeMapArray<Real> & quad_points_coordinates_inv_matrices,
        ElementTypeMapArray<Real> & result, const GhostType ghost_type,
        const ElementTypeMapArray<UInt> * element_filter) const {

  shape_functions.interpolateElementalFieldFromIntegrationPoints(
      field, interpolation_points_coordinates_matrices,
      quad_points_coordinates_inv_matrices, result, ghost_type, element_filter);
}

/* -------------------------------------------------------------------------- */
/**
 * Helper class to be able to write a partial specialization on the element kind
 */
template <ElementKind kind> struct InterpolateHelper {
  template <class S>
  static void call(__attribute__((unused)) const S & shape_functions,
                   __attribute__((unused)) const Vector<Real> & real_coords,
                   __attribute__((unused)) UInt elem,
                   __attribute__((unused)) const Matrix<Real> & nodal_values,
                   __attribute__((unused)) Vector<Real> & interpolated,
                   __attribute__((unused)) const ElementType & type,
                   __attribute__((unused)) const GhostType & ghost_type) {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }
};

#define INTERPOLATE(type)                                                      \
  shape_functions.template interpolate<type>(                                  \
      real_coords, element, nodal_values, interpolated, ghost_type);

#define AKANTU_SPECIALIZE_INTERPOLATE_HELPER(kind)                             \
  template <> struct InterpolateHelper<kind> {                                 \
    template <class S>                                                         \
    static void call(const S & shape_functions,                                \
                     const Vector<Real> & real_coords, UInt element,           \
                     const Matrix<Real> & nodal_values,                        \
                     Vector<Real> & interpolated, const ElementType & type,    \
                     const GhostType & ghost_type) {                           \
      AKANTU_BOOST_KIND_ELEMENT_SWITCH(INTERPOLATE, kind);                     \
    }                                                                          \
  };

AKANTU_BOOST_ALL_KIND_LIST(AKANTU_SPECIALIZE_INTERPOLATE_HELPER,
                           AKANTU_FE_ENGINE_LIST_INTERPOLATE)

#undef AKANTU_SPECIALIZE_INTERPOLATE_HELPER
#undef INTERPOLATE

template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
inline void FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::interpolate(
    const Vector<Real> & real_coords, const Matrix<Real> & nodal_values,
    Vector<Real> & interpolated, const Element & element) const {

  AKANTU_DEBUG_IN();

  InterpolateHelper<kind>::call(shape_functions, real_coords, element.element,
                                nodal_values, interpolated, element.type,
                                element.ghost_type);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
void FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::
    computeNormalsOnIntegrationPoints(const GhostType & ghost_type) {
  AKANTU_DEBUG_IN();

  computeNormalsOnIntegrationPoints(mesh.getNodes(), ghost_type);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
void FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::
    computeNormalsOnIntegrationPoints(const Array<Real> & field,
                                      const GhostType & ghost_type) {
  AKANTU_DEBUG_IN();

  //  Real * coord = mesh.getNodes().storage();
  UInt spatial_dimension = mesh.getSpatialDimension();

  // allocate the normal arrays
  Mesh::type_iterator it = mesh.firstType(element_dimension, ghost_type, kind);
  Mesh::type_iterator end = mesh.lastType(element_dimension, ghost_type, kind);
  for (; it != end; ++it) {
    ElementType type = *it;
    UInt size = mesh.getNbElement(type, ghost_type);
    if (normals_on_integration_points.exists(type, ghost_type)) {
      normals_on_integration_points(type, ghost_type).resize(size);
    } else {
      normals_on_integration_points.alloc(size, spatial_dimension, type,
                                          ghost_type);
    }
  }

  // loop over the type to build the normals
  it = mesh.firstType(element_dimension, ghost_type, kind);
  for (; it != end; ++it) {
    Array<Real> & normals_on_quad =
        normals_on_integration_points(*it, ghost_type);
    computeNormalsOnIntegrationPoints(field, normals_on_quad, *it, ghost_type);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/**
 * Helper class to be able to write a partial specialization on the element kind
 */
template <ElementKind kind> struct ComputeNormalsOnIntegrationPoints {
  template <template <ElementKind, class> class I,
            template <ElementKind> class S, ElementKind k, class IOF>
  static void call(__attribute__((unused))
                   const FEEngineTemplate<I, S, k, IOF> & fem,
                   __attribute__((unused)) const Array<Real> & field,
                   __attribute__((unused)) Array<Real> & normal,
                   __attribute__((unused)) const ElementType & type,
                   __attribute__((unused)) const GhostType & ghost_type) {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }
};

#define COMPUTE_NORMALS_ON_INTEGRATION_POINTS(type)                            \
  fem.template computeNormalsOnIntegrationPoints<type>(field, normal,          \
                                                       ghost_type);

#define AKANTU_SPECIALIZE_COMPUTE_NORMALS_ON_INTEGRATION_POINTS(kind)          \
  template <> struct ComputeNormalsOnIntegrationPoints<kind> {                 \
    template <template <ElementKind, class> class I,                           \
              template <ElementKind> class S, ElementKind k, class IOF>        \
    static void call(const FEEngineTemplate<I, S, k, IOF> & fem,               \
                     const Array<Real> & field, Array<Real> & normal,          \
                     const ElementType & type, const GhostType & ghost_type) { \
      AKANTU_BOOST_KIND_ELEMENT_SWITCH(COMPUTE_NORMALS_ON_INTEGRATION_POINTS,  \
                                       kind);                                  \
    }                                                                          \
  };

AKANTU_BOOST_ALL_KIND_LIST(
    AKANTU_SPECIALIZE_COMPUTE_NORMALS_ON_INTEGRATION_POINTS,
    AKANTU_FE_ENGINE_LIST_COMPUTE_NORMALS_ON_INTEGRATION_POINTS)

#undef AKANTU_SPECIALIZE_COMPUTE_NORMALS_ON_INTEGRATION_POINTS
#undef COMPUTE_NORMALS_ON_INTEGRATION_POINTS

template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
void FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::
    computeNormalsOnIntegrationPoints(const Array<Real> & field,
                                      Array<Real> & normal,
                                      const ElementType & type,
                                      const GhostType & ghost_type) const {
  ComputeNormalsOnIntegrationPoints<kind>::call(*this, field, normal, type,
                                                ghost_type);
}

/* -------------------------------------------------------------------------- */
template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
template <ElementType type>
void FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::
    computeNormalsOnIntegrationPoints(const Array<Real> & field,
                                      Array<Real> & normal,
                                      const GhostType & ghost_type) const {
  AKANTU_DEBUG_IN();

  UInt spatial_dimension = mesh.getSpatialDimension();
  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
  UInt nb_points = getNbIntegrationPoints(type, ghost_type);

  UInt nb_element = mesh.getConnectivity(type, ghost_type).getSize();
  normal.resize(nb_element * nb_points);
  Array<Real>::matrix_iterator normals_on_quad =
      normal.begin_reinterpret(spatial_dimension, nb_points, nb_element);
  Array<Real> f_el(0, spatial_dimension * nb_nodes_per_element);
  FEEngine::extractNodalToElementField(mesh, field, f_el, type, ghost_type);

  const Matrix<Real> & quads =
      integrator.template getIntegrationPoints<type>(ghost_type);

  Array<Real>::matrix_iterator f_it =
      f_el.begin(spatial_dimension, nb_nodes_per_element);

  for (UInt elem = 0; elem < nb_element; ++elem) {
    ElementClass<type>::computeNormalsOnNaturalCoordinates(quads, *f_it,
                                                           *normals_on_quad);
    ++normals_on_quad;
    ++f_it;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/* Matrix lumping functions                                                   */
/* -------------------------------------------------------------------------- */

/**
 * Helper class to be able to write a partial specialization on the element kind
 */
template <ElementKind kind> struct AssembleLumpedTemplateHelper {};

#define ASSEMBLE_LUMPED(type)                                                  \
  fem.template assembleLumpedTemplate<type>(field_1, lumped, dof_id,           \
                                            dof_manager, ghost_type)

#define AKANTU_SPECIALIZE_ASSEMBLE_HELPER(kind)                                \
  template <> struct AssembleLumpedTemplateHelper<kind> {                      \
    template <template <ElementKind, class> class I,                           \
              template <ElementKind> class S, ElementKind k, class IOF>        \
    static void call(const FEEngineTemplate<I, S, k, IOF> & fem,               \
                     const Array<Real> & field_1, const ID & lumped,           \
                     const ID & dof_id, DOFManager & dof_manager,              \
                     ElementType type, const GhostType & ghost_type) {         \
      AKANTU_BOOST_KIND_ELEMENT_SWITCH(ASSEMBLE_LUMPED, kind);                 \
    }                                                                          \
  };

AKANTU_BOOST_ALL_KIND(AKANTU_SPECIALIZE_ASSEMBLE_HELPER)

#undef AKANTU_SPECIALIZE_ASSEMBLE_HELPER
#undef ASSEMBLE_LUMPED

/* -------------------------------------------------------------------------- */
template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IOF>
void FEEngineTemplate<I, S, kind, IOF>::assembleFieldLumped(
    const Array<Real> & field, const ID & lumped, const ID & dof_id,
    DOFManager & dof_manager, ElementType type,
    const GhostType & ghost_type) const {
  AKANTU_DEBUG_IN();

  AssembleLumpedTemplateHelper<kind>::call(*this, field, lumped, dof_id,
                                           dof_manager, type, ghost_type);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/**
 * Helper class to be able to write a partial specialization on the element kind
 */
template <ElementKind kind> struct AssembleFieldMatrixHelper {};

#define ASSEMBLE_MATRIX(type)                                                  \
  fem.template assembleFieldMatrix<Functor, type>(                             \
      field_funct, matrix_id, dof_id, dof_manager, ghost_type)

#define AKANTU_SPECIALIZE_ASSEMBLE_FIELD_MATRIX_HELPER(kind)                   \
  template <> struct AssembleFieldMatrixHelper<kind> {                         \
    template <template <ElementKind, class> class I,                           \
              template <ElementKind> class S, ElementKind k, class IOF,        \
              class Functor>                                                   \
    static void call(const FEEngineTemplate<I, S, k, IOF> & fem,               \
                     Functor field_funct, const ID & matrix_id,                \
                     const ID & dof_id, DOFManager & dof_manager,              \
                     ElementType type, const GhostType & ghost_type) {         \
      AKANTU_BOOST_KIND_ELEMENT_SWITCH(ASSEMBLE_MATRIX, kind);                 \
    }                                                                          \
  };

AKANTU_BOOST_ALL_KIND(AKANTU_SPECIALIZE_ASSEMBLE_FIELD_MATRIX_HELPER)

#undef AKANTU_SPECIALIZE_ASSEMBLE_FIELD_MATRIX_HELPER
#undef ASSEMBLE_MATRIX

template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IOF>
template <class Functor>
void FEEngineTemplate<I, S, kind, IOF>::assembleFieldMatrix(
    Functor field_funct, const ID & matrix_id, const ID & dof_id,
    DOFManager & dof_manager, ElementType type,
    const GhostType & ghost_type) const {
  AKANTU_DEBUG_IN();
  AssembleFieldMatrixHelper<kind>::template call(
      *this, field_funct, matrix_id, dof_id, dof_manager, type, ghost_type);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
template <ElementType type>
void FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::
    assembleLumpedTemplate(const Array<Real> & field, const ID & lumped,
                           const ID & dof_id, DOFManager & dof_manager,
                           const GhostType & ghost_type) const {
  AKANTU_DEBUG_IN();
  this->template assembleLumpedRowSum<type>(field, lumped, dof_id, dof_manager,
                                            ghost_type);
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/**
 * @f$ \tilde{M}_{i} = \sum_j M_{ij} = \sum_j \int \rho \varphi_i \varphi_j dV =
 * \int \rho \varphi_i dV @f$
 */
template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
template <ElementType type>
void FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::
    assembleLumpedRowSum(const Array<Real> & field, const ID & lumped,
                         const ID & dof_id, DOFManager & dof_manager,
                         const GhostType & ghost_type) const {
  AKANTU_DEBUG_IN();

  UInt shapes_size = ElementClass<type>::getShapeSize();
  UInt nb_degree_of_freedom = field.getNbComponent();

  Array<Real> * field_times_shapes =
      new Array<Real>(0, shapes_size * nb_degree_of_freedom);

  shape_functions.template fieldTimesShapes<type>(field, *field_times_shapes,
                                                  ghost_type);

  UInt nb_element = mesh.getNbElement(type, ghost_type);
  Array<Real> * int_field_times_shapes = new Array<Real>(
      nb_element, shapes_size * nb_degree_of_freedom, "inte_rho_x_shapes");

  integrator.template integrate<type>(
      *field_times_shapes, *int_field_times_shapes,
      nb_degree_of_freedom * shapes_size, ghost_type, empty_filter);

  delete field_times_shapes;

  dof_manager.assembleElementalArrayToLumpedMatrix(
      dof_id, *int_field_times_shapes, lumped, type, ghost_type);

  delete int_field_times_shapes;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/**
 * @f$ \tilde{M}_{i} = c * M_{ii} = \int_{V_e} \rho dV @f$
 */
template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
template <ElementType type>
void FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::
    assembleLumpedDiagonalScaling(const Array<Real> & field, const ID & lumped,
                                  const ID & dof_id, DOFManager & dof_manager,
                                  const GhostType & ghost_type) const {
  AKANTU_DEBUG_IN();

  const ElementType & type_p1 = ElementClass<type>::getP1ElementType();
  UInt nb_nodes_per_element_p1 = Mesh::getNbNodesPerElement(type_p1);
  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
  UInt nb_quadrature_points =
      integrator.template getIntegrationPoints<type>(ghost_type).cols();
  UInt nb_degree_of_freedom = field.getNbComponent();

  UInt nb_element = field.getSize() / nb_quadrature_points;
  Vector<Real> nodal_factor(nb_nodes_per_element);

#define ASSIGN_WEIGHT_TO_NODES(corner, mid)                                    \
  {                                                                            \
    for (UInt n = 0; n < nb_nodes_per_element_p1; n++)                         \
      nodal_factor(n) = corner;                                                \
    for (UInt n = nb_nodes_per_element_p1; n < nb_nodes_per_element; n++)      \
      nodal_factor(n) = mid;                                                   \
  }

  if (type == _triangle_6)
    ASSIGN_WEIGHT_TO_NODES(1. / 12., 1. / 4.);
  if (type == _tetrahedron_10)
    ASSIGN_WEIGHT_TO_NODES(1. / 32., 7. / 48.);
  if (type == _quadrangle_8)
    ASSIGN_WEIGHT_TO_NODES(
        3. / 76.,
        16. / 76.); /** coeff. derived by scaling
                     * the diagonal terms of the corresponding
                     * consistent mass computed with 3x3 gauss points;
                     * coeff. are (1./36., 8./36.) for 2x2 gauss points */
  if (type == _hexahedron_20)
    ASSIGN_WEIGHT_TO_NODES(
        7. / 248.,
        16. / 248.); /** coeff. derived by scaling
                      * the diagonal terms of the corresponding
                      * consistent mass computed with 3x3x3 gauss points;
                      * coeff. are (1./40.,
                      * 1./15.) for 2x2x2 gauss points */
  if (type == _pentahedron_15) {
    // coefficients derived by scaling the diagonal terms of the corresponding
    // consistent mass computed with 8 gauss points;
    for (UInt n = 0; n < nb_nodes_per_element_p1; n++)
      nodal_factor(n) = 51. / 2358.;

    Real mid_triangle = 192. / 2358.;
    Real mid_quadrangle = 300. / 2358.;

    nodal_factor(6) = mid_triangle;
    nodal_factor(7) = mid_triangle;
    nodal_factor(8) = mid_triangle;
    nodal_factor(9) = mid_quadrangle;
    nodal_factor(10) = mid_quadrangle;
    nodal_factor(11) = mid_quadrangle;
    nodal_factor(12) = mid_triangle;
    nodal_factor(13) = mid_triangle;
    nodal_factor(14) = mid_triangle;
  }

  if (nb_element == 0) {
    AKANTU_DEBUG_OUT();
    return;
  }

#undef ASSIGN_WEIGHT_TO_NODES

  /// compute @f$ \int \rho dV = \rho V @f$ for each element
  Array<Real> * int_field =
      new Array<Real>(field.getSize(), nb_degree_of_freedom, "inte_rho_x");
  integrator.template integrate<type>(field, *int_field, nb_degree_of_freedom,
                                      ghost_type, empty_filter);

  /// distribute the mass of the element to the nodes
  Array<Real> * lumped_per_node = new Array<Real>(
      nb_element, nb_degree_of_freedom * nb_nodes_per_element, "mass_per_node");
  Array<Real>::const_vector_iterator int_field_it =
      int_field->begin(nb_degree_of_freedom);
  Array<Real>::matrix_iterator lumped_per_node_it =
      lumped_per_node->begin(nb_degree_of_freedom, nb_nodes_per_element);

  for (UInt e = 0; e < nb_element; ++e) {
    for (UInt n = 0; n < nb_nodes_per_element; ++n) {
      Vector<Real> l = (*lumped_per_node_it)(n);
      l = *int_field_it;
      l *= nodal_factor(n);
    }
    ++int_field_it;
    ++lumped_per_node_it;
  }
  delete int_field;

  dof_manager.assembleElementalArrayToLumpedMatrix(dof_id, *lumped_per_node,
                                                   lumped, type, ghost_type);
  delete lumped_per_node;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/**
 * @f$ \tilde{M}_{i} = \sum_j M_{ij} = \sum_j \int \rho \varphi_i \varphi_j dV =
 * \int \rho \varphi_i dV @f$
 */
template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
template <class Functor, ElementType type>
void FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::assembleFieldMatrix(
    Functor field_funct, const ID & matrix_id, const ID & dof_id,
    DOFManager & dof_manager, const GhostType & ghost_type) const {
  AKANTU_DEBUG_IN();

  UInt shapes_size = ElementClass<type>::getShapeSize();
  UInt nb_degree_of_freedom = dof_manager.getDOFs(dof_id).getNbComponent();
  UInt lmat_size = nb_degree_of_freedom * shapes_size;
  UInt nb_element = mesh.getNbElement(type, ghost_type);

  // \int N * N  so degree 2 * degree of N
  const UInt polynomial_degree =
      2 * ElementClassProperty<type>::polynomial_degree;

  Matrix<Real> integration_points =
      integrator.template getIntegrationPoints<type, polynomial_degree>();

  UInt nb_integration_points = integration_points.cols();
  UInt vect_size = nb_integration_points * nb_element;

  Array<Real> shapes(0, shapes_size);
  shape_functions.template computeShapesOnIntegrationPoints<type>(
      mesh.getNodes(), integration_points, shapes, ghost_type);

  Array<Real> integration_points_pos(vect_size, mesh.getSpatialDimension());
  shape_functions.template interpolateOnIntegrationPoints<type>(
      mesh.getNodes(), integration_points_pos, mesh.getSpatialDimension(),
      shapes, ghost_type, empty_filter);

  Array<Real> * modified_shapes =
      new Array<Real>(vect_size, lmat_size * nb_degree_of_freedom);
  modified_shapes->clear();
  Array<Real> * local_mat = new Array<Real>(vect_size, lmat_size * lmat_size);

  Array<Real>::matrix_iterator mshapes_it =
      modified_shapes->begin(nb_degree_of_freedom, lmat_size);
  Array<Real>::const_vector_iterator shapes_it = shapes.begin(shapes_size);

  for (UInt q = 0; q < vect_size; ++q, ++mshapes_it, ++shapes_it) {
    for (UInt d = 0; d < nb_degree_of_freedom; ++d) {
      for (UInt s = 0; s < shapes_size; ++s) {
        (*mshapes_it)(d, s * nb_degree_of_freedom + d) = (*shapes_it)(s);
      }
    }
  }

  Array<Real> field(vect_size, nb_degree_of_freedom);
  Array<Real>::matrix_iterator field_c_it = field.begin_reinterpret(
      nb_degree_of_freedom, nb_integration_points, nb_element);
  Array<Real>::const_matrix_iterator pos_it =
      integration_points_pos.begin_reinterpret(
          mesh.getSpatialDimension(), nb_integration_points, nb_element);
  Element el;
  el.type = type, el.ghost_type = ghost_type;

  for (el.element = 0; el.element < nb_element;
       ++el.element, ++field_c_it, ++pos_it) {
    field_funct(*field_c_it, el, *pos_it);
  }

  mshapes_it = modified_shapes->begin(nb_degree_of_freedom, lmat_size);
  Array<Real>::matrix_iterator lmat = local_mat->begin(lmat_size, lmat_size);
  Array<Real>::const_vector_iterator field_it =
      field.begin_reinterpret(nb_degree_of_freedom, field.getSize());

  for (UInt q = 0; q < vect_size; ++q, ++lmat, ++mshapes_it, ++field_it) {
    const Vector<Real> & rho = *field_it;
    const Matrix<Real> & N = *mshapes_it;
    Matrix<Real> & mat = *lmat;

    Matrix<Real> Nt = N.transpose();
    for (UInt d = 0; d < Nt.cols(); ++d) {
      Nt(d) *= rho(d);
    }

    mat.mul<false, false>(Nt, N);
  }

  delete modified_shapes;

  Array<Real> * int_field_times_shapes =
      new Array<Real>(nb_element, lmat_size * lmat_size, "inte_rho_x_shapes");
  this->integrator.template integrate<type, polynomial_degree>(
      *local_mat, *int_field_times_shapes, lmat_size * lmat_size, ghost_type);
  delete local_mat;

  dof_manager.assembleElementalMatricesToMatrix(
      matrix_id, dof_id, *int_field_times_shapes, type, ghost_type);
  delete int_field_times_shapes;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/**
 * Helper class to be able to write a partial specialization on the element kind
 */
template <ElementKind kind> struct InverseMapHelper {
  template <class S>
  static void call(__attribute__((unused)) const S & shape_functions,
                   __attribute__((unused)) const Vector<Real> & real_coords,
                   __attribute__((unused)) UInt element,
                   __attribute__((unused)) const ElementType & type,
                   __attribute__((unused)) Vector<Real> & natural_coords,
                   __attribute__((unused)) const GhostType & ghost_type) {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }
};

#define INVERSE_MAP(type)                                                      \
  shape_functions.template inverseMap<type>(real_coords, element,              \
                                            natural_coords, ghost_type);

#define AKANTU_SPECIALIZE_INVERSE_MAP_HELPER(kind)                             \
  template <> struct InverseMapHelper<kind> {                                  \
    template <class S>                                                         \
    static void call(const S & shape_functions,                                \
                     const Vector<Real> & real_coords, UInt element,           \
                     const ElementType & type, Vector<Real> & natural_coords,  \
                     const GhostType & ghost_type) {                           \
      AKANTU_BOOST_KIND_ELEMENT_SWITCH(INVERSE_MAP, kind);                     \
    }                                                                          \
  };

AKANTU_BOOST_ALL_KIND_LIST(AKANTU_SPECIALIZE_INVERSE_MAP_HELPER,
                           AKANTU_FE_ENGINE_LIST_INVERSE_MAP)

#undef AKANTU_SPECIALIZE_INVERSE_MAP_HELPER
#undef INVERSE_MAP

template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
inline void FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::inverseMap(
    const Vector<Real> & real_coords, UInt element, const ElementType & type,
    Vector<Real> & natural_coords, const GhostType & ghost_type) const {

  AKANTU_DEBUG_IN();

  InverseMapHelper<kind>::call(shape_functions, real_coords, element, type,
                               natural_coords, ghost_type);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/**
 * Helper class to be able to write a partial specialization on the element kind
 */
template <ElementKind kind> struct ContainsHelper {
  template <class S>
  static void call(__attribute__((unused)) const S & shape_functions,
                   __attribute__((unused)) const Vector<Real> & real_coords,
                   __attribute__((unused)) UInt element,
                   __attribute__((unused)) const ElementType & type,
                   __attribute__((unused)) const GhostType & ghost_type) {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }
};

#define CONTAINS(type)                                                         \
  contain = shape_functions.template contains<type>(real_coords, element,      \
                                                    ghost_type);

#define AKANTU_SPECIALIZE_CONTAINS_HELPER(kind)                                \
  template <> struct ContainsHelper<kind> {                                    \
    template <template <ElementKind> class S, ElementKind k>                   \
    static bool call(const S<k> & shape_functions,                             \
                     const Vector<Real> & real_coords, UInt element,           \
                     const ElementType & type, const GhostType & ghost_type) { \
      bool contain = false;                                                    \
      AKANTU_BOOST_KIND_ELEMENT_SWITCH(CONTAINS, kind);                        \
      return contain;                                                          \
    }                                                                          \
  };

AKANTU_BOOST_ALL_KIND_LIST(AKANTU_SPECIALIZE_CONTAINS_HELPER,
                           AKANTU_FE_ENGINE_LIST_CONTAINS)

#undef AKANTU_SPECIALIZE_CONTAINS_HELPER
#undef CONTAINS

template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
inline bool FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::contains(
    const Vector<Real> & real_coords, UInt element, const ElementType & type,
    const GhostType & ghost_type) const {
  return ContainsHelper<kind>::call(shape_functions, real_coords, element, type,
                                    ghost_type);
}

/* -------------------------------------------------------------------------- */
/**
 * Helper class to be able to write a partial specialization on the element kind
 */
template <ElementKind kind> struct ComputeShapesHelper {
  template <class S>
  static void call(__attribute__((unused)) const S & shape_functions,
                   __attribute__((unused)) const Vector<Real> & real_coords,
                   __attribute__((unused)) UInt element,
                   __attribute__((unused)) const ElementType type,
                   __attribute__((unused)) Vector<Real> & shapes,
                   __attribute__((unused)) const GhostType & ghost_type) {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }
};

#define COMPUTE_SHAPES(type)                                                   \
  shape_functions.template computeShapes<type>(real_coords, element, shapes,   \
                                               ghost_type);

#define AKANTU_SPECIALIZE_COMPUTE_SHAPES_HELPER(kind)                          \
  template <> struct ComputeShapesHelper<kind> {                               \
    template <class S>                                                         \
    static void call(const S & shape_functions,                                \
                     const Vector<Real> & real_coords, UInt element,           \
                     const ElementType type, Vector<Real> & shapes,            \
                     const GhostType & ghost_type) {                           \
      AKANTU_BOOST_KIND_ELEMENT_SWITCH(COMPUTE_SHAPES, kind);                  \
    }                                                                          \
  };

AKANTU_BOOST_ALL_KIND_LIST(AKANTU_SPECIALIZE_COMPUTE_SHAPES_HELPER,
                           AKANTU_FE_ENGINE_LIST_COMPUTE_SHAPES)

#undef AKANTU_SPECIALIZE_COMPUTE_SHAPES_HELPER
#undef COMPUTE_SHAPES

template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
inline void
FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::computeShapes(
    const Vector<Real> & real_coords, UInt element, const ElementType & type,
    Vector<Real> & shapes, const GhostType & ghost_type) const {

  AKANTU_DEBUG_IN();

  ComputeShapesHelper<kind>::call(shape_functions, real_coords, element, type,
                                  shapes, ghost_type);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/**
 * Helper class to be able to write a partial specialization on the element kind
 */
template <ElementKind kind> struct ComputeShapeDerivativesHelper {
  template <class S>
  static void call(__attribute__((unused)) const S & shape_functions,
                   __attribute__((unused)) const Vector<Real> & real_coords,
                   __attribute__((unused)) UInt element,
                   __attribute__((unused)) const ElementType type,
                   __attribute__((unused)) Matrix<Real> & shape_derivatives,
                   __attribute__((unused)) const GhostType & ghost_type) {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }
};

#define COMPUTE_SHAPE_DERIVATIVES(type)                                        \
  Matrix<Real> coords_mat(real_coords.storage(), shape_derivatives.rows(), 1); \
  Tensor3<Real> shapesd_tensor(shape_derivatives.storage(),                    \
                               shape_derivatives.rows(),                       \
                               shape_derivatives.cols(), 1);                   \
  shape_functions.template computeShapeDerivatives<type>(                      \
      coords_mat, element, shapesd_tensor, ghost_type);

#define AKANTU_SPECIALIZE_COMPUTE_SHAPE_DERIVATIVES_HELPER(kind)               \
  template <> struct ComputeShapeDerivativesHelper<kind> {                     \
    template <class S>                                                         \
    static void call(const S & shape_functions,                                \
                     const Vector<Real> & real_coords, UInt element,           \
                     const ElementType type, Matrix<Real> & shape_derivatives, \
                     const GhostType & ghost_type) {                           \
      AKANTU_BOOST_KIND_ELEMENT_SWITCH(COMPUTE_SHAPE_DERIVATIVES, kind);       \
    }                                                                          \
  };

AKANTU_BOOST_ALL_KIND_LIST(AKANTU_SPECIALIZE_COMPUTE_SHAPE_DERIVATIVES_HELPER,
                           AKANTU_FE_ENGINE_LIST_COMPUTE_SHAPES_DERIVATIVES)

#undef AKANTU_SPECIALIZE_COMPUTE_SHAPE_DERIVATIVES_HELPER
#undef COMPUTE_SHAPE_DERIVATIVES

template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
inline void
FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::computeShapeDerivatives(
    const Vector<Real> & real_coords, UInt element, const ElementType & type,
    Matrix<Real> & shape_derivatives, const GhostType & ghost_type) const {
  AKANTU_DEBUG_IN();

  ComputeShapeDerivativesHelper<kind>::call(shape_functions, real_coords,
                                            element, type, shape_derivatives,
                                            ghost_type);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/**
 * Helper class to be able to write a partial specialization on the element kind
 */
template <ElementKind kind> struct GetNbIntegrationPointsHelper {};

#define GET_NB_INTEGRATION_POINTS(type)                                        \
  nb_quad_points = integrator.template getNbIntegrationPoints<type>(ghost_type);

#define AKANTU_SPECIALIZE_GET_NB_INTEGRATION_POINTS_HELPER(kind)               \
  template <> struct GetNbIntegrationPointsHelper<kind> {                      \
    template <template <ElementKind, class> class I, ElementKind k, class IOF> \
    static UInt call(const I<k, IOF> & integrator, const ElementType type,     \
                     const GhostType & ghost_type) {                           \
      UInt nb_quad_points = 0;                                                 \
      AKANTU_BOOST_KIND_ELEMENT_SWITCH(GET_NB_INTEGRATION_POINTS, kind);       \
      return nb_quad_points;                                                   \
    }                                                                          \
  };

AKANTU_BOOST_ALL_KIND(AKANTU_SPECIALIZE_GET_NB_INTEGRATION_POINTS_HELPER)

#undef AKANTU_SPECIALIZE_GET_NB_INTEGRATION_POINTS_HELPER
#undef GET_NB_INTEGRATION

template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
inline UInt
FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::getNbIntegrationPoints(
    const ElementType & type, const GhostType & ghost_type) const {
  return GetNbIntegrationPointsHelper<kind>::call(integrator, type, ghost_type);
}

/* -------------------------------------------------------------------------- */
/**
 * Helper class to be able to write a partial specialization on the element kind
 */
template <ElementKind kind> struct GetShapesHelper {};

#define GET_SHAPES(type) ret = &(shape_functions.getShapes(type, ghost_type));

#define AKANTU_SPECIALIZE_GET_SHAPES_HELPER(kind)                              \
  template <> struct GetShapesHelper<kind> {                                   \
    template <class S>                                                         \
    static const Array<Real> & call(const S & shape_functions,                 \
                                    const ElementType type,                    \
                                    const GhostType & ghost_type) {            \
      const Array<Real> * ret = NULL;                                          \
      AKANTU_BOOST_KIND_ELEMENT_SWITCH(GET_SHAPES, kind);                      \
      return *ret;                                                             \
    }                                                                          \
  };

AKANTU_BOOST_ALL_KIND(AKANTU_SPECIALIZE_GET_SHAPES_HELPER)

#undef AKANTU_SPECIALIZE_GET_SHAPES_HELPER
#undef GET_SHAPES

template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
inline const Array<Real> &
FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::getShapes(
    const ElementType & type, const GhostType & ghost_type,
    __attribute__((unused)) UInt id) const {
  return GetShapesHelper<kind>::call(shape_functions, type, ghost_type);
}

/* -------------------------------------------------------------------------- */
/**
 * Helper class to be able to write a partial specialization on the element kind
 */
template <ElementKind kind> struct GetShapesDerivativesHelper {
  template <template <ElementKind> class S, ElementKind k>
  static const Array<Real> &
  call(__attribute__((unused)) const S<k> & shape_functions,
       __attribute__((unused)) const ElementType & type,
       __attribute__((unused)) const GhostType & ghost_type,
       __attribute__((unused)) UInt id) {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }
};

#define GET_SHAPES_DERIVATIVES(type)                                           \
  ret = &(shape_functions.getShapesDerivatives(type, ghost_type));

#define AKANTU_SPECIALIZE_GET_SHAPES_DERIVATIVES_HELPER(kind)                  \
  template <> struct GetShapesDerivativesHelper<kind> {                        \
    template <template <ElementKind> class S, ElementKind k>                   \
    static const Array<Real> &                                                 \
    call(const S<k> & shape_functions, const ElementType type,                 \
         const GhostType & ghost_type, __attribute__((unused)) UInt id) {      \
      const Array<Real> * ret = NULL;                                          \
      AKANTU_BOOST_KIND_ELEMENT_SWITCH(GET_SHAPES_DERIVATIVES, kind);          \
      return *ret;                                                             \
    }                                                                          \
  };

AKANTU_BOOST_ALL_KIND_LIST(AKANTU_SPECIALIZE_GET_SHAPES_DERIVATIVES_HELPER,
                           AKANTU_FE_ENGINE_LIST_GET_SHAPES_DERIVATIVES)

#undef AKANTU_SPECIALIZE_GET_SHAPE_DERIVATIVES_HELPER
#undef GET_SHAPES_DERIVATIVES

template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
inline const Array<Real> &
FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::getShapesDerivatives(
    const ElementType & type, const GhostType & ghost_type,
    __attribute__((unused)) UInt id) const {
  return GetShapesDerivativesHelper<kind>::call(shape_functions, type,
                                                ghost_type, id);
}

/* -------------------------------------------------------------------------- */
/**
 * Helper class to be able to write a partial specialization on the element kind
 */
template <ElementKind kind> struct GetIntegrationPointsHelper {};

#define GET_INTEGRATION_POINTS(type)                                           \
  ret = &(integrator.template getIntegrationPoints<type>(ghost_type));

#define AKANTU_SPECIALIZE_GET_INTEGRATION_POINTS_HELPER(kind)                  \
  template <> struct GetIntegrationPointsHelper<kind> {                        \
    template <template <ElementKind, class> class I, ElementKind k, class IOF> \
    static const Matrix<Real> & call(const I<k, IOF> & integrator,             \
                                     const ElementType type,                   \
                                     const GhostType & ghost_type) {           \
      const Matrix<Real> * ret = NULL;                                         \
      AKANTU_BOOST_KIND_ELEMENT_SWITCH(GET_INTEGRATION_POINTS, kind);          \
      return *ret;                                                             \
    }                                                                          \
  };

AKANTU_BOOST_ALL_KIND(AKANTU_SPECIALIZE_GET_INTEGRATION_POINTS_HELPER)

#undef AKANTU_SPECIALIZE_GET_INTEGRATION_POINTS_HELPER
#undef GET_INTEGRATION_POINTS

template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
inline const Matrix<Real> &
FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::getIntegrationPoints(
    const ElementType & type, const GhostType & ghost_type) const {
  return GetIntegrationPointsHelper<kind>::call(integrator, type, ghost_type);
}

/* -------------------------------------------------------------------------- */
template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
void FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::printself(
    std::ostream & stream, int indent) const {
  std::string space;
  for (Int i = 0; i < indent; i++, space += AKANTU_INDENT)
    ;

  stream << space << "FEEngineTemplate [" << std::endl;
  stream << space << " + parent [" << std::endl;
  FEEngine::printself(stream, indent + 3);
  stream << space << "   ]" << std::endl;
  stream << space << " + shape functions [" << std::endl;
  shape_functions.printself(stream, indent + 3);
  stream << space << "   ]" << std::endl;
  stream << space << " + integrator [" << std::endl;
  integrator.printself(stream, indent + 3);
  stream << space << "   ]" << std::endl;
  stream << space << "]" << std::endl;
}

/* -------------------------------------------------------------------------- */

} // akantu

#include "integrator_gauss.hh"
#include "shape_lagrange.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
template <>
template <>
inline void FEEngineTemplate<IntegratorGauss, ShapeLagrange, _ek_regular>::
    assembleLumpedTemplate<_triangle_6>(const Array<Real> & field,
                                        const ID & lumped, const ID & dof_id,
                                        DOFManager & dof_manager,
                                        const GhostType & ghost_type) const {
  assembleLumpedDiagonalScaling<_triangle_6>(field, lumped, dof_id, dof_manager,
                                             ghost_type);
}

/* -------------------------------------------------------------------------- */
template <>
template <>
inline void FEEngineTemplate<IntegratorGauss, ShapeLagrange, _ek_regular,
                             DefaultIntegrationOrderFunctor>::
    assembleLumpedTemplate<_tetrahedron_10>(
        const Array<Real> & field, const ID & lumped, const ID & dof_id,
        DOFManager & dof_manager, const GhostType & ghost_type) const {
  assembleLumpedDiagonalScaling<_tetrahedron_10>(field, lumped, dof_id,
                                                 dof_manager, ghost_type);
}

/* -------------------------------------------------------------------------- */
template <>
template <>
inline void FEEngineTemplate<IntegratorGauss, ShapeLagrange, _ek_regular>::
    assembleLumpedTemplate<_quadrangle_8>(const Array<Real> & field,
                                          const ID & lumped, const ID & dof_id,
                                          DOFManager & dof_manager,
                                          const GhostType & ghost_type) const {
  assembleLumpedDiagonalScaling<_quadrangle_8>(field, lumped, dof_id,
                                               dof_manager, ghost_type);
}

/* -------------------------------------------------------------------------- */
template <>
template <>
inline void FEEngineTemplate<IntegratorGauss, ShapeLagrange, _ek_regular>::
    assembleLumpedTemplate<_hexahedron_20>(const Array<Real> & field,
                                           const ID & lumped, const ID & dof_id,
                                           DOFManager & dof_manager,
                                           const GhostType & ghost_type) const {
  assembleLumpedDiagonalScaling<_hexahedron_20>(field, lumped, dof_id,
                                                dof_manager, ghost_type);
}

/* -------------------------------------------------------------------------- */
template <>
template <>
inline void FEEngineTemplate<IntegratorGauss, ShapeLagrange, _ek_regular,
                             DefaultIntegrationOrderFunctor>::
    assembleLumpedTemplate<_pentahedron_15>(
        const Array<Real> & field, const ID & lumped, const ID & dof_id,
        DOFManager & dof_manager, const GhostType & ghost_type) const {
  assembleLumpedDiagonalScaling<_pentahedron_15>(field, lumped, dof_id,
                                                 dof_manager, ghost_type);
}

/* -------------------------------------------------------------------------- */
template <>
template <>
inline void FEEngineTemplate<IntegratorGauss, ShapeLagrange, _ek_regular,
                             DefaultIntegrationOrderFunctor>::
    computeNormalsOnIntegrationPoints<_point_1>(
        const Array<Real> &, Array<Real> & normal,
        const GhostType & ghost_type) const {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_ASSERT(mesh.getSpatialDimension() == 1,
                      "Mesh dimension must be 1 to compute normals on points!");
  const ElementType type = _point_1;
  UInt spatial_dimension = mesh.getSpatialDimension();
  // UInt nb_nodes_per_element  = Mesh::getNbNodesPerElement(type);
  UInt nb_points = getNbIntegrationPoints(type, ghost_type);

  UInt nb_element = mesh.getConnectivity(type, ghost_type).getSize();
  normal.resize(nb_element * nb_points);
  Array<Real>::matrix_iterator normals_on_quad =
      normal.begin_reinterpret(spatial_dimension, nb_points, nb_element);
  const Array<std::vector<Element>> & segments =
      mesh.getElementToSubelement(type, ghost_type);
  const Array<Real> & coords = mesh.getNodes();

  const Mesh * mesh_segment;
  if (mesh.isMeshFacets())
    mesh_segment = &(mesh.getMeshParent());
  else
    mesh_segment = &mesh;

  for (UInt elem = 0; elem < nb_element; ++elem) {
    UInt nb_segment = segments(elem).size();

    AKANTU_DEBUG_ASSERT(
        nb_segment > 0,
        "Impossible to compute a normal on a point connected to 0 segments");

    Real normal_value = 1;
    if (nb_segment == 1) {
      const Element & segment = segments(elem)[0];
      const Array<UInt> & segment_connectivity =
          mesh_segment->getConnectivity(segment.type, segment.ghost_type);
      // const Vector<UInt> & segment_points =
      // segment_connectivity.begin(Mesh::getNbNodesPerElement(segment.type))[segment.element];
      Real difference;
      if (segment_connectivity(0) == elem) {
        difference = coords(elem) - coords(segment_connectivity(1));
      } else {
        difference = coords(elem) - coords(segment_connectivity(0));
      }
      normal_value = difference / std::abs(difference);
    }

    for (UInt n(0); n < nb_points; ++n) {
      (*normals_on_quad)(0, n) = normal_value;
    }
    ++normals_on_quad;
  }

  AKANTU_DEBUG_OUT();
}

} // akantu

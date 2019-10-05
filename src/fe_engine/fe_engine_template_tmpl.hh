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
 * @date last modification: Tue Feb 20 2018
 *
 * @brief  Template implementation of FEEngineTemplate
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "aka_common.hh"
#include "dof_manager.hh"
#include "fe_engine_template.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::FEEngineTemplate(
    Mesh & mesh, UInt spatial_dimension, ID id, MemoryID memory_id)
    : FEEngine(mesh, spatial_dimension, id, memory_id),
      integrator(mesh, spatial_dimension, id, memory_id),
      shape_functions(mesh, spatial_dimension, id, memory_id) {}

/* -------------------------------------------------------------------------- */
template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::~FEEngineTemplate() =
    default;

/* -------------------------------------------------------------------------- */
/**
 * Helper class to be able to write a partial specialization on the element kind
 */
namespace fe_engine {
  namespace details {
    template <ElementKind kind> struct GradientOnIntegrationPointsHelper {
      template <class S>
      static void call(const S &, Mesh &, const Array<Real> &, Array<Real> &,
                       const UInt, const ElementType &, const GhostType &,
                       const Array<UInt> &) {
        AKANTU_TO_IMPLEMENT();
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
  } // namespace details
} // namespace fe_engine

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
    nb_element = filter_elements.size();
  UInt nb_points =
      shape_functions.getIntegrationPoints(type, ghost_type).cols();

#ifndef AKANTU_NDEBUG

  UInt element_dimension = mesh.getSpatialDimension(type);

  AKANTU_DEBUG_ASSERT(u.size() == mesh.getNbNodes(),
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

// AKANTU_DEBUG_ASSERT(nablauq.size() == nb_element * nb_points,
//                  "The vector nablauq(" << nablauq.getID()
//                  << ") has not the good size.");
#endif

  nablauq.resize(nb_element * nb_points);

  fe_engine::details::GradientOnIntegrationPointsHelper<kind>::call(
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

  for (auto & type : mesh.elementTypes(element_dimension, ghost_type, kind)) {
    integrator.initIntegrator(nodes, type, ghost_type);
    const auto & control_points = getIntegrationPoints(type, ghost_type);
    shape_functions.initShapeFunctions(nodes, control_points, type, ghost_type);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/**
 * Helper class to be able to write a partial specialization on the element kind
 */
namespace fe_engine {
  namespace details {
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
  } // namespace details
} // namespace fe_engine

template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
void FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::integrate(
    const Array<Real> & f, Array<Real> & intf, UInt nb_degree_of_freedom,
    const ElementType & type, const GhostType & ghost_type,
    const Array<UInt> & filter_elements) const {

  UInt nb_element = mesh.getNbElement(type, ghost_type);
  if (filter_elements != empty_filter)
    nb_element = filter_elements.size();
#ifndef AKANTU_NDEBUG

  UInt nb_quadrature_points = getNbIntegrationPoints(type);

  AKANTU_DEBUG_ASSERT(f.size() == nb_element * nb_quadrature_points,
                      "The vector f(" << f.getID() << " size " << f.size()
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

  fe_engine::details::IntegrateHelper<kind>::call(integrator, f, intf,
                                                  nb_degree_of_freedom, type,
                                                  ghost_type, filter_elements);
}

/* -------------------------------------------------------------------------- */
/**
 * Helper class to be able to write a partial specialization on the element kind
 */
namespace fe_engine {
  namespace details {
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
  } // namespace details
} // namespace fe_engine

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
    nb_element = filter_elements.size();

  UInt nb_quadrature_points = getNbIntegrationPoints(type, ghost_type);

  AKANTU_DEBUG_ASSERT(
      f.size() == nb_element * nb_quadrature_points,
      "The vector f(" << f.getID() << ") has not the good size. (" << f.size()
                      << "!=" << nb_quadrature_points * nb_element << ")");
  AKANTU_DEBUG_ASSERT(f.getNbComponent() == 1,
                      "The vector f("
                          << f.getID()
                          << ") has not the good number of component.");
#endif

  Real integral = fe_engine::details::IntegrateScalarHelper<kind>::call(
      integrator, f, type, ghost_type, filter_elements);
  AKANTU_DEBUG_OUT();
  return integral;
}

/* -------------------------------------------------------------------------- */
/**
 * Helper class to be able to write a partial specialization on the element kind
 */
namespace fe_engine {
  namespace details {
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

    AKANTU_BOOST_ALL_KIND(
        AKANTU_SPECIALIZE_INTEGRATE_SCALAR_ON_ONE_ELEMENT_HELPER)

#undef AKANTU_SPECIALIZE_INTEGRATE_SCALAR_ON_ONE_ELEMENT_HELPER
#undef INTEGRATE
  } // namespace details
} // namespace fe_engine

template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
Real FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::integrate(
    const Vector<Real> & f, const ElementType & type, UInt index,
    const GhostType & ghost_type) const {

  Real res = fe_engine::details::IntegrateScalarOnOneElementHelper<kind>::call(
      integrator, f, type, index, ghost_type);
  return res;
}

/* -------------------------------------------------------------------------- */
/**
 * Helper class to be able to write a partial specialization on the element kind
 */
namespace fe_engine {
  namespace details {
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

    AKANTU_BOOST_ALL_KIND(
        AKANTU_SPECIALIZE_INTEGRATE_ON_INTEGRATION_POINTS_HELPER)

#undef AKANTU_SPECIALIZE_INTEGRATE_ON_INTEGRATION_POINTS_HELPER
#undef INTEGRATE
  } // namespace details
} // namespace fe_engine

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
    nb_element = filter_elements.size();
  UInt nb_quadrature_points = getNbIntegrationPoints(type);
#ifndef AKANTU_NDEBUG
  //   std::stringstream sstr; sstr << ghost_type;
  //   AKANTU_DEBUG_ASSERT(sstr.str() == nablauq.getTag(),
  //                  "The vector " << nablauq.getID() << " is not taged " <<
  //                  ghost_type);

  AKANTU_DEBUG_ASSERT(f.size() == nb_element * nb_quadrature_points,
                      "The vector f(" << f.getID() << " size " << f.size()
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
  fe_engine::details::IntegrateOnIntegrationPointsHelper<kind>::call(
      integrator, f, intf, nb_degree_of_freedom, type, ghost_type,
      filter_elements);
}

/* -------------------------------------------------------------------------- */
/**
 * Helper class to be able to write a partial specialization on the element kind
 */
namespace fe_engine {
  namespace details {
    template <ElementKind kind> struct InterpolateOnIntegrationPointsHelper {
      template <class S>
      static void call(const S &, const Array<Real> &, Array<Real> &,
                       const UInt, const ElementType &, const GhostType &,
                       const Array<UInt> &) {
        AKANTU_TO_IMPLEMENT();
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
  } // namespace details
} // namespace fe_engine

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
    nb_element = filter_elements.size();
#ifndef AKANTU_NDEBUG

  AKANTU_DEBUG_ASSERT(u.size() == mesh.getNbNodes(),
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

  fe_engine::details::InterpolateOnIntegrationPointsHelper<kind>::call(
      shape_functions, u, uq, nb_degree_of_freedom, type, ghost_type,
      filter_elements);

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

  const Array<UInt> * filter = nullptr;

  for (auto ghost_type : ghost_types) {
    for (auto & type : uq.elementTypes(_all_dimensions, ghost_type, kind)) {
      UInt nb_quad_per_element = getNbIntegrationPoints(type, ghost_type);

      UInt nb_element = 0;

      if (filter_elements) {
        filter = &((*filter_elements)(type, ghost_type));
        nb_element = filter->size();
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
namespace fe_engine {
  namespace details {
    template <ElementKind kind> struct ComputeBtDHelper {};

#define COMPUTE_BTD(type)                                                      \
  shape_functions.template computeBtD<type>(Ds, BtDs, ghost_type,              \
                                            filter_elements);

#define AKANTU_SPECIALIZE_COMPUTE_BtD_HELPER(kind)                             \
  template <> struct ComputeBtDHelper<kind> {                                  \
    template <class S>                                                         \
    static void call(const S & shape_functions, const Array<Real> & Ds,        \
                     Array<Real> & BtDs, const ElementType & type,             \
                     const GhostType & ghost_type,                             \
                     const Array<UInt> & filter_elements) {                    \
      AKANTU_BOOST_KIND_ELEMENT_SWITCH(COMPUTE_BTD, kind);                     \
    }                                                                          \
  };

    AKANTU_BOOST_ALL_KIND(AKANTU_SPECIALIZE_COMPUTE_BtD_HELPER)

#undef AKANTU_SPECIALIZE_COMPUTE_BtD_HELPER
#undef COMPUTE_BTD
  } // namespace details
} // namespace fe_engine

template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
inline void FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::computeBtD(
    const Array<Real> & Ds, Array<Real> & BtDs, const ElementType & type,
    const GhostType & ghost_type, const Array<UInt> & filter_elements) const {
  fe_engine::details::ComputeBtDHelper<kind>::call(
      shape_functions, Ds, BtDs, type, ghost_type, filter_elements);
}

/* -------------------------------------------------------------------------- */
namespace fe_engine {
  namespace details {
    template <ElementKind kind> struct ComputeBtDBHelper {};

#define COMPUTE_BTDB(type)                                                     \
  shape_functions.template computeBtDB<type>(Ds, BtDBs, order_d, ghost_type,   \
                                             filter_elements);

#define AKANTU_SPECIALIZE_COMPUTE_BtDB_HELPER(kind)                            \
  template <> struct ComputeBtDBHelper<kind> {                                 \
    template <class S>                                                         \
    static void call(const S & shape_functions, const Array<Real> & Ds,        \
                     Array<Real> & BtDBs, UInt order_d,                        \
                     const ElementType & type, const GhostType & ghost_type,   \
                     const Array<UInt> & filter_elements) {                    \
      AKANTU_BOOST_KIND_ELEMENT_SWITCH(COMPUTE_BTDB, kind);                    \
    }                                                                          \
  };

    AKANTU_BOOST_ALL_KIND(AKANTU_SPECIALIZE_COMPUTE_BtDB_HELPER)

#undef AKANTU_SPECIALIZE_COMPUTE_BtDB_HELPER
#undef COMPUTE_BTDB
  } // namespace details
} // namespace fe_engine

template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
inline void FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::computeBtDB(
    const Array<Real> & Ds, Array<Real> & BtDBs, UInt order_d,
    const ElementType & type, const GhostType & ghost_type,
    const Array<UInt> & filter_elements) const {
  fe_engine::details::ComputeBtDBHelper<kind>::call(
      shape_functions, Ds, BtDBs, order_d, type, ghost_type, filter_elements);
}

/* -------------------------------------------------------------------------- */
namespace fe_engine {
  namespace details {
    template <ElementKind kind> struct ComputeNtbNHelper {};

#define COMPUTE_NTbN(type)                                                     \
  shape_functions.template computeNtbN<type>(bs, NtbNs, order_d, ghost_type,   \
                                             filter_elements);

#define AKANTU_SPECIALIZE_COMPUTE_NtbN_HELPER(kind)                            \
  template <> struct ComputeNtbNHelper<kind> {                                 \
    template <class S>                                                         \
    static void call(const S & shape_functions, const Array<Real> & bs,        \
                     Array<Real> & NtbNs, UInt order_d,                        \
                     const ElementType & type, const GhostType & ghost_type,   \
                     const Array<UInt> & filter_elements) {                    \
      AKANTU_BOOST_KIND_ELEMENT_SWITCH(COMPUTE_NTbN, kind);                    \
    }                                                                          \
  };

    AKANTU_BOOST_ALL_KIND(AKANTU_SPECIALIZE_COMPUTE_NtbN_HELPER)

#undef AKANTU_SPECIALIZE_COMPUTE_NtbN_HELPER
#undef COMPUTE_NTbN
  } // namespace details
} // namespace fe_engine

template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
inline void FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::computeNtbN(
    const Array<Real> & bs, Array<Real> & NtbNs, UInt order_d,
    const ElementType & type, const GhostType & ghost_type,
    const Array<UInt> & filter_elements) const {
  fe_engine::details::ComputeNtbNHelper<kind>::call(
      shape_functions, bs, NtbNs, order_d, type, ghost_type, filter_elements);
}

  
/* -------------------------------------------------------------------------- */
namespace fe_engine {
  namespace details {
    template <ElementKind kind> struct ComputeNtbHelper {};

#define COMPUTE_Ntb(type)                                                      \
  shape_functions.template computeNtb<type>(bs, Ntbs, ghost_type,              \
                                            filter_elements);

#define AKANTU_SPECIALIZE_COMPUTE_Ntb_HELPER(kind)                             \
  template <> struct ComputeNtbHelper<kind> {                                  \
    template <class S>                                                         \
    static void call(const S & shape_functions, const Array<Real> & bs,        \
                     Array<Real> & Ntbs, const ElementType & type,             \
                     const GhostType & ghost_type,                             \
                     const Array<UInt> & filter_elements) {                    \
      AKANTU_BOOST_KIND_ELEMENT_SWITCH(COMPUTE_Ntb, kind);                     \
    }                                                                          \
  };

    AKANTU_BOOST_ALL_KIND(AKANTU_SPECIALIZE_COMPUTE_Ntb_HELPER)

#undef AKANTU_SPECIALIZE_COMPUTE_Ntb_HELPER
#undef COMPUTE_Ntb
  } // namespace details
} // namespace fe_engine

template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
inline void FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::computeNtb(
    const Array<Real> & bs, Array<Real> & Ntbs, const ElementType & type,
    const GhostType & ghost_type, const Array<UInt> & filter_elements) const {
  fe_engine::details::ComputeNtbHelper<kind>::call(
      shape_functions, bs, Ntbs, type, ghost_type, filter_elements);
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

  quadrature_points_coordinates.initialize(*this,
                                           _nb_component = spatial_dimension);

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
namespace fe_engine {
  namespace details {
    template <ElementKind kind> struct InterpolateHelper {
      template <class S>
      static void call(const S &, const Vector<Real> &, UInt,
                       const Matrix<Real> &, Vector<Real> &,
                       const ElementType &, const GhostType &) {
        AKANTU_TO_IMPLEMENT();
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
  } // namespace details
} // namespace fe_engine

template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
inline void FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::interpolate(
    const Vector<Real> & real_coords, const Matrix<Real> & nodal_values,
    Vector<Real> & interpolated, const Element & element) const {

  AKANTU_DEBUG_IN();

  fe_engine::details::InterpolateHelper<kind>::call(
      shape_functions, real_coords, element.element, nodal_values, interpolated,
      element.type, element.ghost_type);

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
  normals_on_integration_points.initialize(
      *this, _nb_component = spatial_dimension,
      _spatial_dimension = element_dimension, _ghost_type = ghost_type,
      _element_kind = kind);

  // loop over the type to build the normals
  for (auto & type : mesh.elementTypes(element_dimension, ghost_type, kind)) {
    auto & normals_on_quad = normals_on_integration_points(type, ghost_type);
    computeNormalsOnIntegrationPoints(field, normals_on_quad, type, ghost_type);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/**
 * Helper class to be able to write a partial specialization on the element kind
 */
namespace fe_engine {
  namespace details {
    template <ElementKind kind> struct ComputeNormalsOnIntegrationPoints {
      template <template <ElementKind, class> class I,
                template <ElementKind> class S, ElementKind k, class IOF>
      static void call(const FEEngineTemplate<I, S, k, IOF> &,
                       const Array<Real> &, Array<Real> &, const ElementType &,
                       const GhostType &) {
        AKANTU_TO_IMPLEMENT();
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
  } // namespace details
} // namespace fe_engine

template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
void FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::
    computeNormalsOnIntegrationPoints(const Array<Real> & field,
                                      Array<Real> & normal,
                                      const ElementType & type,
                                      const GhostType & ghost_type) const {
  fe_engine::details::ComputeNormalsOnIntegrationPoints<kind>::call(
      *this, field, normal, type, ghost_type);
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

  if (type == _point_1) {
    computeNormalsOnIntegrationPointsPoint1(field, normal, ghost_type);
    return;
  }

  UInt spatial_dimension = mesh.getSpatialDimension();
  UInt nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
  UInt nb_points = getNbIntegrationPoints(type, ghost_type);

  UInt nb_element = mesh.getConnectivity(type, ghost_type).size();
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
/**
 * Helper class to be able to write a partial specialization on the element kind
 */
template <ElementKind kind> struct InverseMapHelper {
  template <class S>
  static void
  call(const S & /*shape_functions*/, const Vector<Real> & /*real_coords*/,
       UInt /*element*/, const ElementType & /*type*/,
       Vector<Real> & /*natural_coords*/, const GhostType & /*ghost_type*/) {
    AKANTU_TO_IMPLEMENT();
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
namespace fe_engine {
  namespace details {
    template <ElementKind kind> struct ContainsHelper {
      template <class S>
      static void call(const S &, const Vector<Real> &, UInt,
                       const ElementType &, const GhostType &) {
        AKANTU_TO_IMPLEMENT();
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
  } // namespace details
} // namespace fe_engine

template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
inline bool FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::contains(
    const Vector<Real> & real_coords, UInt element, const ElementType & type,
    const GhostType & ghost_type) const {
  return fe_engine::details::ContainsHelper<kind>::call(
      shape_functions, real_coords, element, type, ghost_type);
}

/* -------------------------------------------------------------------------- */
/**
 * Helper class to be able to write a partial specialization on the element kind
 */
namespace fe_engine {
  namespace details {
    template <ElementKind kind> struct ComputeShapesHelper {
      template <class S>
      static void call(const S &, const Vector<Real> &, UInt, const ElementType,
                       Vector<Real> &, const GhostType &) {
        AKANTU_TO_IMPLEMENT();
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
  } // namespace details
} // namespace fe_engine
template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
inline void
FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::computeShapes(
    const Vector<Real> & real_coords, UInt element, const ElementType & type,
    Vector<Real> & shapes, const GhostType & ghost_type) const {

  AKANTU_DEBUG_IN();

  fe_engine::details::ComputeShapesHelper<kind>::call(
      shape_functions, real_coords, element, type, shapes, ghost_type);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/**
 * Helper class to be able to write a partial specialization on the element kind
 */
namespace fe_engine {
  namespace details {
    template <ElementKind kind> struct ComputeShapeDerivativesHelper {
      template <class S>
      static void call(__attribute__((unused)) const S & shape_functions,
                       __attribute__((unused)) const Vector<Real> & real_coords,
                       __attribute__((unused)) UInt element,
                       __attribute__((unused)) const ElementType type,
                       __attribute__((unused)) Matrix<Real> & shape_derivatives,
                       __attribute__((unused)) const GhostType & ghost_type) {
        AKANTU_TO_IMPLEMENT();
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

    AKANTU_BOOST_ALL_KIND_LIST(
        AKANTU_SPECIALIZE_COMPUTE_SHAPE_DERIVATIVES_HELPER,
        AKANTU_FE_ENGINE_LIST_COMPUTE_SHAPES_DERIVATIVES)

#undef AKANTU_SPECIALIZE_COMPUTE_SHAPE_DERIVATIVES_HELPER
#undef COMPUTE_SHAPE_DERIVATIVES
  } // namespace details
} // namespace fe_engine

template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
inline void
FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::computeShapeDerivatives(
    const Vector<Real> & real_coords, UInt element, const ElementType & type,
    Matrix<Real> & shape_derivatives, const GhostType & ghost_type) const {
  AKANTU_DEBUG_IN();

  fe_engine::details::ComputeShapeDerivativesHelper<kind>::call(
      shape_functions, real_coords, element, type, shape_derivatives,
      ghost_type);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/**
 * Helper class to be able to write a partial specialization on the element kind
 */
namespace fe_engine {
  namespace details {
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
  } // namespace details
} // namespace fe_engine

template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
inline UInt
FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::getNbIntegrationPoints(
    const ElementType & type, const GhostType & ghost_type) const {
  return fe_engine::details::GetNbIntegrationPointsHelper<kind>::call(
      integrator, type, ghost_type);
}

/* -------------------------------------------------------------------------- */
/**
 * Helper class to be able to write a partial specialization on the element kind
 */
namespace fe_engine {
  namespace details {
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
  } // namespace details
} // namespace fe_engine

template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
inline const Array<Real> &
FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::getShapes(
    const ElementType & type, const GhostType & ghost_type,
    __attribute__((unused)) UInt id) const {
  return fe_engine::details::GetShapesHelper<kind>::call(shape_functions, type,
                                                         ghost_type);
}

/* -------------------------------------------------------------------------- */
/**
 * Helper class to be able to write a partial specialization on the element kind
 */
namespace fe_engine {
  namespace details {
    template <ElementKind kind> struct GetShapesDerivativesHelper {
      template <template <ElementKind> class S, ElementKind k>
      static const Array<Real> & call(const S<k> &, const ElementType &,
                                      const GhostType &, UInt) {
        AKANTU_TO_IMPLEMENT();
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
  } // namespace details
} // namespace fe_engine

template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
inline const Array<Real> &
FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::getShapesDerivatives(
    const ElementType & type, const GhostType & ghost_type,
    __attribute__((unused)) UInt id) const {
  return fe_engine::details::GetShapesDerivativesHelper<kind>::call(
      shape_functions, type, ghost_type, id);
}

/* -------------------------------------------------------------------------- */
/**
 * Helper class to be able to write a partial specialization on the element kind
 */
namespace fe_engine {
  namespace details {
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
  } // namespace details
} // namespace fe_engine

template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
inline const Matrix<Real> &
FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::getIntegrationPoints(
    const ElementType & type, const GhostType & ghost_type) const {
  return fe_engine::details::GetIntegrationPointsHelper<kind>::call(
      integrator, type, ghost_type);
}

/* -------------------------------------------------------------------------- */
template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
void FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::printself(
    std::ostream & stream, int indent) const {
  std::string space(indent, AKANTU_INDENT);

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
template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
void FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::onElementsAdded(
    const Array<Element> & new_elements, const NewElementsEvent &) {
  integrator.onElementsAdded(new_elements);
  shape_functions.onElementsAdded(new_elements);
}

/* -------------------------------------------------------------------------- */
template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
void FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::onElementsRemoved(
    const Array<Element> &, const ElementTypeMapArray<UInt> &,
    const RemovedElementsEvent &) {}

/* -------------------------------------------------------------------------- */
template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
void FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::onElementsChanged(
    const Array<Element> &, const Array<Element> &,
    const ElementTypeMapArray<UInt> &, const ChangedElementsEvent &) {}

/* -------------------------------------------------------------------------- */
template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind, class IntegrationOrderFunctor>
inline void FEEngineTemplate<I, S, kind, IntegrationOrderFunctor>::
    computeNormalsOnIntegrationPointsPoint1(
        const Array<Real> &, Array<Real> & normal,
        const GhostType & ghost_type) const {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_ASSERT(mesh.getSpatialDimension() == 1,
                      "Mesh dimension must be 1 to compute normals on points!");
  const auto type = _point_1;
  auto spatial_dimension = mesh.getSpatialDimension();
  // UInt nb_nodes_per_element  = Mesh::getNbNodesPerElement(type);
  auto nb_points = getNbIntegrationPoints(type, ghost_type);
  const auto & connectivity = mesh.getConnectivity(type, ghost_type);
  auto nb_element = connectivity.size();

  normal.resize(nb_element * nb_points);
  auto normals_on_quad =
      normal.begin_reinterpret(spatial_dimension, nb_points, nb_element);
  const auto & segments = mesh.getElementToSubelement(type, ghost_type);
  const auto & coords = mesh.getNodes();

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
      auto point = connectivity(elem);
      const auto segment = segments(elem)[0];
      const auto & segment_connectivity =
          mesh_segment->getConnectivity(segment.type, segment.ghost_type);
      Vector<UInt> segment_points = segment_connectivity.begin(
          Mesh::getNbNodesPerElement(segment.type))[segment.element];
      Real difference;
      if (segment_points(0) == point) {
        difference = coords(elem) - coords(segment_points(1));
      } else {
        difference = coords(elem) - coords(segment_points(0));
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

} // namespace akantu

/**
 * @file   fe_engine_template.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Mon Jan 29 2018
 *
 * @brief  templated class that calls integration and shape objects
 *
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
#include "fe_engine.hh"
#include "integrator.hh"
#include "shape_functions.hh"
/* -------------------------------------------------------------------------- */
#include <type_traits>
/* -------------------------------------------------------------------------- */


/* -------------------------------------------------------------------------- */
#ifndef __AKANTU_FE_ENGINE_TEMPLATE_HH__
#define __AKANTU_FE_ENGINE_TEMPLATE_HH__

namespace akantu {
class DOFManager;
namespace fe_engine {
  namespace details {
    template <ElementKind> struct AssembleLumpedTemplateHelper;
    template <ElementKind> struct AssembleFieldMatrixHelper;
  } // namespace details
} // namespace fe_engine

template <ElementKind, typename> struct AssembleFieldMatrixStructHelper;

struct DefaultIntegrationOrderFunctor {
  template <ElementType type> static inline constexpr int getOrder() {
    return ElementClassProperty<type>::polynomial_degree;
  }
};

/* -------------------------------------------------------------------------- */
template <template <ElementKind, class> class I, template <ElementKind> class S,
          ElementKind kind = _ek_regular,
          class IntegrationOrderFunctor = DefaultIntegrationOrderFunctor>
class FEEngineTemplate : public FEEngine {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  using Integ = I<kind, IntegrationOrderFunctor>;
  using Shape = S<kind>;

  FEEngineTemplate(Mesh & mesh, UInt spatial_dimension = _all_dimensions,
                   ID id = "fem", MemoryID memory_id = 0);

  ~FEEngineTemplate() override;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// pre-compute all the shape functions, their derivatives and the jacobians
  void initShapeFunctions(const GhostType & ghost_type = _not_ghost) override;
  void initShapeFunctions(const Array<Real> & nodes,
                          const GhostType & ghost_type = _not_ghost);

  /* ------------------------------------------------------------------------ */
  /* Integration method bridges                                               */
  /* ------------------------------------------------------------------------ */
  /// integrate f for all elements of type "type"
  void
  integrate(const Array<Real> & f, Array<Real> & intf,
            UInt nb_degree_of_freedom, const ElementType & type,
            const GhostType & ghost_type = _not_ghost,
            const Array<UInt> & filter_elements = empty_filter) const override;

  /// integrate a scalar value on all elements of type "type"
  Real
  integrate(const Array<Real> & f, const ElementType & type,
            const GhostType & ghost_type = _not_ghost,
            const Array<UInt> & filter_elements = empty_filter) const override;

  /// integrate one element scalar value on all elements of type "type"
  Real integrate(const Vector<Real> & f, const ElementType & type, UInt index,
                 const GhostType & ghost_type = _not_ghost) const override;

  /// integrate partially around an integration point (@f$ intf_q = f_q * J_q *
  /// w_q @f$)
  void integrateOnIntegrationPoints(
      const Array<Real> & f, Array<Real> & intf, UInt nb_degree_of_freedom,
      const ElementType & type, const GhostType & ghost_type = _not_ghost,
      const Array<UInt> & filter_elements = empty_filter) const override;

  /// interpolate on a phyiscal point inside an element
  void interpolate(const Vector<Real> & real_coords,
                   const Matrix<Real> & nodal_values,
                   Vector<Real> & interpolated,
                   const Element & element) const override;

  /// get the number of integration points
  UInt getNbIntegrationPoints(
      const ElementType & type,
      const GhostType & ghost_type = _not_ghost) const override;

  /// get shapes precomputed
  const Array<Real> & getShapes(const ElementType & type,
                                const GhostType & ghost_type = _not_ghost,
                                UInt id = 0) const override;

  /// get the derivatives of shapes
  const Array<Real> &
  getShapesDerivatives(const ElementType & type,
                       const GhostType & ghost_type = _not_ghost,
                       UInt id = 0) const override;

  /// get integration points
  const inline Matrix<Real> & getIntegrationPoints(
      const ElementType & type,
      const GhostType & ghost_type = _not_ghost) const override;

  /* ------------------------------------------------------------------------ */
  /* Shape method bridges                                                     */
  /* ------------------------------------------------------------------------ */

  /// compute the gradient of a nodal field on the integration points
  void gradientOnIntegrationPoints(
      const Array<Real> & u, Array<Real> & nablauq,
      const UInt nb_degree_of_freedom, const ElementType & type,
      const GhostType & ghost_type = _not_ghost,
      const Array<UInt> & filter_elements = empty_filter) const override;

  /// interpolate a nodal field on the integration points
  void interpolateOnIntegrationPoints(
      const Array<Real> & u, Array<Real> & uq, UInt nb_degree_of_freedom,
      const ElementType & type, const GhostType & ghost_type = _not_ghost,
      const Array<UInt> & filter_elements = empty_filter) const override;

  /// interpolate a nodal field on the integration points given a
  /// by_element_type
  void interpolateOnIntegrationPoints(
      const Array<Real> & u, ElementTypeMapArray<Real> & uq,
      const ElementTypeMapArray<UInt> * filter_elements =
          nullptr) const override;

  /// pre multiplies a tensor by the shapes derivaties
  void
  computeBtD(const Array<Real> & Ds, Array<Real> & BtDs,
             const ElementType & type, const GhostType & ghost_type,
             const Array<UInt> & filter_elements = empty_filter) const override;

  /// left and right  multiplies a tensor by the shapes derivaties
  void computeBtDB(
      const Array<Real> & Ds, Array<Real> & BtDBs, UInt order_d,
      const ElementType & type, const GhostType & ghost_type,
      const Array<UInt> & filter_elements = empty_filter) const override;

  /// left multiples a vector by the shape functions
  void computeNtb(const Array<Real> & bs, Array<Real> & Ntbs,
                  const ElementType & type, const GhostType & ghost_type,
                  const Array<UInt> & filter_elements) const override;

  /// compute the position of integration points given by an element_type_map
  /// from nodes position
  inline void computeIntegrationPointsCoordinates(
      ElementTypeMapArray<Real> & quadrature_points_coordinates,
      const ElementTypeMapArray<UInt> * filter_elements =
          nullptr) const override;

  /// compute the position of integration points from nodes position
  inline void computeIntegrationPointsCoordinates(
      Array<Real> & quadrature_points_coordinates, const ElementType & type,
      const GhostType & ghost_type = _not_ghost,
      const Array<UInt> & filter_elements = empty_filter) const override;

  /// interpolate field at given position (interpolation_points) from given
  /// values of this field at integration points (field)
  inline void interpolateElementalFieldFromIntegrationPoints(
      const ElementTypeMapArray<Real> & field,
      const ElementTypeMapArray<Real> & interpolation_points_coordinates,
      ElementTypeMapArray<Real> & result, const GhostType ghost_type,
      const ElementTypeMapArray<UInt> * element_filter) const override;

  /// Interpolate field at given position from given values of this field at
  /// integration points (field)
  /// using matrices precomputed with
  /// initElementalFieldInterplationFromIntegrationPoints
  inline void interpolateElementalFieldFromIntegrationPoints(
      const ElementTypeMapArray<Real> & field,
      const ElementTypeMapArray<Real> &
          interpolation_points_coordinates_matrices,
      const ElementTypeMapArray<Real> & quad_points_coordinates_inv_matrices,
      ElementTypeMapArray<Real> & result, const GhostType ghost_type,
      const ElementTypeMapArray<UInt> * element_filter) const override;

  /// Build pre-computed matrices for interpolation of field form integration
  /// points at other given positions (interpolation_points)
  inline void initElementalFieldInterpolationFromIntegrationPoints(
      const ElementTypeMapArray<Real> & interpolation_points_coordinates,
      ElementTypeMapArray<Real> & interpolation_points_coordinates_matrices,
      ElementTypeMapArray<Real> & quad_points_coordinates_inv_matrices,
      const ElementTypeMapArray<UInt> * element_filter =
          nullptr) const override;

  /// find natural coords from real coords provided an element
  void inverseMap(const Vector<Real> & real_coords, UInt element,
                  const ElementType & type, Vector<Real> & natural_coords,
                  const GhostType & ghost_type = _not_ghost) const;

  /// return true if the coordinates provided are inside the element, false
  /// otherwise
  inline bool contains(const Vector<Real> & real_coords, UInt element,
                       const ElementType & type,
                       const GhostType & ghost_type = _not_ghost) const;

  /// compute the shape on a provided point
  inline void
  computeShapes(const Vector<Real> & real_coords, UInt element,
                const ElementType & type, Vector<Real> & shapes,
                const GhostType & ghost_type = _not_ghost) const override;

  /// compute the shape derivatives on a provided point
  inline void computeShapeDerivatives(
      const Vector<Real> & real__coords, UInt element, const ElementType & type,
      Matrix<Real> & shape_derivatives,
      const GhostType & ghost_type = _not_ghost) const override;

  /* ------------------------------------------------------------------------ */
  /* Other methods                                                            */
  /* ------------------------------------------------------------------------ */
  /// pre-compute normals on integration points
  void computeNormalsOnIntegrationPoints(
      const GhostType & ghost_type = _not_ghost) override;
  void computeNormalsOnIntegrationPoints(
      const Array<Real> & field,
      const GhostType & ghost_type = _not_ghost) override;
  void computeNormalsOnIntegrationPoints(
      const Array<Real> & field, Array<Real> & normal, const ElementType & type,
      const GhostType & ghost_type = _not_ghost) const override;
  template <ElementType type>
  void computeNormalsOnIntegrationPoints(const Array<Real> & field,
                                         Array<Real> & normal,
                                         const GhostType & ghost_type) const;

private:
  // To avoid a weird full specialization of a method in a non specalized class
  void
  computeNormalsOnIntegrationPointsPoint1(const Array<Real> &,
                                          Array<Real> & normal,
                                          const GhostType & ghost_type) const;

public:
  /// function to print the contain of the class
  void printself(std::ostream & stream, int indent = 0) const override;

  void assembleFieldLumped(
      const std::function<void(Matrix<Real> &, const Element &)> & field_funct,
      const ID & matrix_id, const ID & dof_id, DOFManager & dof_manager,
      ElementType type, const GhostType & ghost_type) const override;

  /// assemble a field as a matrix (ex. rho to mass matrix)
  void assembleFieldMatrix(
      const std::function<void(Matrix<Real> &, const Element &)> & field_funct,
      const ID & matrix_id, const ID & dof_id, DOFManager & dof_manager,
      ElementType type, const GhostType & ghost_type) const override;

  /// assemble a field as a lumped matrix (ex. rho in lumped mass)
  // template <class Functor>
  // void assembleFieldLumped(const Functor & field_funct, const ID & matrix_id,
  //                          const ID & dof_id, DOFManager & dof_manager,
  //                          ElementType type,
  //                          const GhostType & ghost_type) const;

  // /// assemble a field as a matrix (ex. rho to mass matrix)
  // template <class Functor>
  // void assembleFieldMatrix(const Functor & field_funct, const ID & matrix_id,
  //                          const ID & dof_id, DOFManager & dof_manager,
  //                          ElementType type,
  //                          const GhostType & ghost_type) const;

  // #ifdef AKANTU_STRUCTURAL_MECHANICS
  //   /// assemble a field as a matrix (ex. rho to mass matrix)
  //   void assembleFieldMatrix(const Array<Real> & field_1,
  //                            UInt nb_degree_of_freedom, SparseMatrix & M,
  //                            Array<Real> * n,
  //                            ElementTypeMapArray<Real> & rotation_mat,
  //                            const ElementType & type,
  //                            const GhostType & ghost_type = _not_ghost)
  //                            const;

  //   /// compute shapes function in a matrix for structural elements
  //   void
  //   computeShapesMatrix(const ElementType & type, UInt nb_degree_of_freedom,
  //                       UInt nb_nodes_per_element, Array<Real> * n, UInt id,
  //                       UInt degree_to_interpolate, UInt degree_interpolated,
  //                       const bool sign,
  //                       const GhostType & ghost_type = _not_ghost) const
  //                       override;
  // #endif

private:
  friend struct fe_engine::details::AssembleLumpedTemplateHelper<kind>;
  friend struct fe_engine::details::AssembleFieldMatrixHelper<kind>;
  friend struct AssembleFieldMatrixStructHelper<kind, void>;

  /// templated function to compute the scaling to assemble a lumped matrix
  template <ElementType type>
  void assembleFieldLumped(
      const std::function<void(Matrix<Real> &, const Element &)> & field_funct,
      const ID & matrix_id, const ID & dof_id, DOFManager & dof_manager,
      const GhostType & ghost_type) const;

  /// @f$ \tilde{M}_{i} = \sum_j M_{ij} = \sum_j \int \rho \varphi_i \varphi_j
  /// dV = \int \rho \varphi_i dV @f$
  template <ElementType type>
  void assembleLumpedRowSum(const Array<Real> & field, const ID & matrix_id,
                            const ID & dof_id, DOFManager & dof_manager,
                            const GhostType & ghost_type) const;

  /// @f$ \tilde{M}_{i} = c * M_{ii} = \int_{V_e} \rho dV @f$
  template <ElementType type>
  void assembleLumpedDiagonalScaling(const Array<Real> & field,
                                     const ID & matrix_id, const ID & dof_id,
                                     DOFManager & dof_manager,
                                     const GhostType & ghost_type) const;

  /// assemble a field as a matrix (ex. rho to mass matrix)
  template <ElementType type>
  void assembleFieldMatrix(
      const std::function<void(Matrix<Real> &, const Element &)> & field_funct,
      const ID & matrix_id, const ID & dof_id, DOFManager & dof_manager,
      const GhostType & ghost_type) const;

#ifdef AKANTU_STRUCTURAL_MECHANICS

  /// assemble a field as a matrix for structural elements (ex. rho to mass
  /// matrix)
  template <ElementType type>
  void assembleFieldMatrix(const Array<Real> & field_1,
                           UInt nb_degree_of_freedom, SparseMatrix & M,
                           Array<Real> * n,
                           ElementTypeMapArray<Real> & rotation_mat,
                           __attribute__((unused))
                           const GhostType & ghost_type) const;

#endif

  /* ------------------------------------------------------------------------ */
  /* Mesh Event Handler interface                                             */
  /* ------------------------------------------------------------------------ */
public:
  void onElementsAdded(const Array<Element> &,
                       const NewElementsEvent &) override;
  void onElementsRemoved(const Array<Element> &,
                         const ElementTypeMapArray<UInt> &,
                         const RemovedElementsEvent &) override;
  void onElementsChanged(const Array<Element> &, const Array<Element> &,
                         const ElementTypeMapArray<UInt> &,
                         const ChangedElementsEvent &) override;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// get the shape class (probably useless: see getShapeFunction)
  const ShapeFunctions & getShapeFunctionsInterface() const override {
    return shape_functions;
  };
  /// get the shape class
  const Shape & getShapeFunctions() const { return shape_functions; };

  /// get the integrator class (probably useless: see getIntegrator)
  const Integrator & getIntegratorInterface() const override {
    return integrator;
  };
  /// get the integrator class
  const Integ & getIntegrator() const { return integrator; };

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  Integ integrator;
  Shape shape_functions;
};

} // namespace akantu

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */
#include "fe_engine_template_tmpl.hh"
#include "fe_engine_template_tmpl_field.hh"
/* -------------------------------------------------------------------------- */
/* Shape Linked specialization                                                */
/* -------------------------------------------------------------------------- */
#if defined(AKANTU_STRUCTURAL_MECHANICS)
#include "fe_engine_template_tmpl_struct.hh"
#endif
/* -------------------------------------------------------------------------- */
/* Shape IGFEM specialization                                                 */
/* -------------------------------------------------------------------------- */
#if defined(AKANTU_IGFEM)
#include "fe_engine_template_tmpl_igfem.hh"
#endif

#endif /* __AKANTU_FE_ENGINE_TEMPLATE_HH__ */

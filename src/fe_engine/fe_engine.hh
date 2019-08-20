/**
 * @file   fe_engine.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Tue Feb 20 2018
 *
 * @brief  FEM class
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
#include "aka_memory.hh"
#include "element_type_map.hh"
#include "mesh_events.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_FE_ENGINE_HH__
#define __AKANTU_FE_ENGINE_HH__

namespace akantu {
class Mesh;
class Integrator;
class ShapeFunctions;
class DOFManager;
class Element;
} // namespace akantu

/* -------------------------------------------------------------------------- */
namespace akantu {
/* -------------------------------------------------------------------------- */

/**
 * The  generic  FEEngine class  derived  in  a  FEEngineTemplate class
 * containing  the
 * shape functions and the integration method
 */
class FEEngine : protected Memory, public MeshEventHandler {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  FEEngine(Mesh & mesh, UInt spatial_dimension = _all_dimensions,
           const ID & id = "fem", MemoryID memory_id = 0);

  ~FEEngine() override;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// pre-compute all the shape functions, their derivatives and the jacobians
  virtual void
  initShapeFunctions(const GhostType & ghost_type = _not_ghost) = 0;

  /// extract the nodal values and store them per element
  template <typename T>
  static void extractNodalToElementField(
      const Mesh & mesh, const Array<T> & nodal_f, Array<T> & elemental_f,
      const ElementType & type, const GhostType & ghost_type = _not_ghost,
      const Array<UInt> & filter_elements = empty_filter);

  /// filter a field
  template <typename T>
  static void
  filterElementalData(const Mesh & mesh, const Array<T> & quad_f,
                      Array<T> & filtered_f, const ElementType & type,
                      const GhostType & ghost_type = _not_ghost,
                      const Array<UInt> & filter_elements = empty_filter);

  /* ------------------------------------------------------------------------ */
  /* Integration method bridges                                               */
  /* ------------------------------------------------------------------------ */
  /// integrate f for all elements of type "type"
  virtual void
  integrate(const Array<Real> & f, Array<Real> & intf,
            UInt nb_degree_of_freedom, const ElementType & type,
            const GhostType & ghost_type = _not_ghost,
            const Array<UInt> & filter_elements = empty_filter) const = 0;

  /// integrate a scalar value f on all elements of type "type"
  virtual Real
  integrate(const Array<Real> & f, const ElementType & type,
            const GhostType & ghost_type = _not_ghost,
            const Array<UInt> & filter_elements = empty_filter) const = 0;

  /// integrate f for all integration points of type "type" but don't sum over
  /// all integration points
  virtual void integrateOnIntegrationPoints(
      const Array<Real> & f, Array<Real> & intf, UInt nb_degree_of_freedom,
      const ElementType & type, const GhostType & ghost_type = _not_ghost,
      const Array<UInt> & filter_elements = empty_filter) const = 0;

  /// integrate one element scalar value on all elements of type "type"
  virtual Real integrate(const Vector<Real> & f, const ElementType & type,
                         UInt index,
                         const GhostType & ghost_type = _not_ghost) const = 0;

  /* ------------------------------------------------------------------------ */
  /* compatibility with old FEEngine fashion */
  /* ------------------------------------------------------------------------ */
  /// get the number of integration points
  virtual UInt
  getNbIntegrationPoints(const ElementType & type,
                         const GhostType & ghost_type = _not_ghost) const = 0;

  /// get the precomputed shapes
  const virtual Array<Real> &
  getShapes(const ElementType & type, const GhostType & ghost_type = _not_ghost,
            UInt id = 0) const = 0;

  /// get the derivatives of shapes
  const virtual Array<Real> &
  getShapesDerivatives(const ElementType & type,
                       const GhostType & ghost_type = _not_ghost,
                       UInt id = 0) const = 0;

  /// get integration points
  const virtual Matrix<Real> &
  getIntegrationPoints(const ElementType & type,
                       const GhostType & ghost_type = _not_ghost) const = 0;

  /* ------------------------------------------------------------------------ */
  /* Shape method bridges                                                     */
  /* ------------------------------------------------------------------------ */
  /// Compute the gradient nablauq on the integration points of an element type
  /// from nodal values u
  virtual void gradientOnIntegrationPoints(
      const Array<Real> & u, Array<Real> & nablauq,
      const UInt nb_degree_of_freedom, const ElementType & type,
      const GhostType & ghost_type = _not_ghost,
      const Array<UInt> & filter_elements = empty_filter) const = 0;

  /// Interpolate a nodal field u at the integration points of an element type
  /// -> uq
  virtual void interpolateOnIntegrationPoints(
      const Array<Real> & u, Array<Real> & uq, UInt nb_degree_of_freedom,
      const ElementType & type, const GhostType & ghost_type = _not_ghost,
      const Array<UInt> & filter_elements = empty_filter) const = 0;

  /// Interpolate a nodal field u at the integration points of many element
  /// types -> uq
  virtual void interpolateOnIntegrationPoints(
      const Array<Real> & u, ElementTypeMapArray<Real> & uq,
      const ElementTypeMapArray<UInt> * filter_elements = nullptr) const = 0;

  /// pre multiplies a tensor by the shapes derivaties
  virtual void
  computeBtD(const Array<Real> & Ds, Array<Real> & BtDs,
             const ElementType & type,
             const GhostType & ghost_type = _not_ghost,
             const Array<UInt> & filter_elements = empty_filter) const = 0;

  /// left and right  multiplies a tensor by the shapes derivaties
  virtual void
  computeBtDB(const Array<Real> & Ds, Array<Real> & BtDBs, UInt order_d,
              const ElementType & type,
              const GhostType & ghost_type = _not_ghost,
              const Array<UInt> & filter_elements = empty_filter) const = 0;

  /// left multiples a vector by the shape functions
  virtual void
  computeNtb(const Array<Real> & bs, Array<Real> & Ntbs,
             const ElementType & type,
             const GhostType & ghost_type = _not_ghost,
             const Array<UInt> & filter_elements = empty_filter) const = 0;

  /// Compute the interpolation point position in the global coordinates for
  /// many element types
  virtual void computeIntegrationPointsCoordinates(
      ElementTypeMapArray<Real> & integration_points_coordinates,
      const ElementTypeMapArray<UInt> * filter_elements = nullptr) const = 0;

  /// Compute the interpolation point position in the global coordinates for an
  /// element type
  virtual void computeIntegrationPointsCoordinates(
      Array<Real> & integration_points_coordinates, const ElementType & type,
      const GhostType & ghost_type = _not_ghost,
      const Array<UInt> & filter_elements = empty_filter) const = 0;

  /// Build pre-computed matrices for interpolation of field form integration
  /// points at other given positions (interpolation_points)
  virtual void initElementalFieldInterpolationFromIntegrationPoints(
      const ElementTypeMapArray<Real> & interpolation_points_coordinates,
      ElementTypeMapArray<Real> & interpolation_points_coordinates_matrices,
      ElementTypeMapArray<Real> & integration_points_coordinates_inv_matrices,
      const ElementTypeMapArray<UInt> * element_filter) const = 0;

  /// interpolate field at given position (interpolation_points) from given
  /// values of this field at integration points (field)
  virtual void interpolateElementalFieldFromIntegrationPoints(
      const ElementTypeMapArray<Real> & field,
      const ElementTypeMapArray<Real> & interpolation_points_coordinates,
      ElementTypeMapArray<Real> & result, const GhostType ghost_type,
      const ElementTypeMapArray<UInt> * element_filter) const = 0;

  /// Interpolate field at given position from given values of this field at
  /// integration points (field)
  /// using matrices precomputed with
  /// initElementalFieldInterplationFromIntegrationPoints
  virtual void interpolateElementalFieldFromIntegrationPoints(
      const ElementTypeMapArray<Real> & field,
      const ElementTypeMapArray<Real> &
          interpolation_points_coordinates_matrices,
      const ElementTypeMapArray<Real> &
          integration_points_coordinates_inv_matrices,
      ElementTypeMapArray<Real> & result, const GhostType ghost_type,
      const ElementTypeMapArray<UInt> * element_filter) const = 0;

  /// interpolate on a phyiscal point inside an element
  virtual void interpolate(const Vector<Real> & real_coords,
                           const Matrix<Real> & nodal_values,
                           Vector<Real> & interpolated,
                           const Element & element) const = 0;

  /// compute the shape on a provided point
  virtual void
  computeShapes(const Vector<Real> & real_coords, UInt elem,
                const ElementType & type, Vector<Real> & shapes,
                const GhostType & ghost_type = _not_ghost) const = 0;

  /// compute the shape derivatives on a provided point
  virtual void
  computeShapeDerivatives(const Vector<Real> & real__coords, UInt element,
                          const ElementType & type,
                          Matrix<Real> & shape_derivatives,
                          const GhostType & ghost_type = _not_ghost) const = 0;

  /* ------------------------------------------------------------------------ */
  /* Other methods                                                            */
  /* ------------------------------------------------------------------------ */

  /// pre-compute normals on integration points
  virtual void computeNormalsOnIntegrationPoints(
      const GhostType & ghost_type = _not_ghost) = 0;

  /// pre-compute normals on integration points
  virtual void computeNormalsOnIntegrationPoints(
      __attribute__((unused)) const Array<Real> & field,
      __attribute__((unused)) const GhostType & ghost_type = _not_ghost) {
    AKANTU_TO_IMPLEMENT();
  }

  /// pre-compute normals on integration points
  virtual void computeNormalsOnIntegrationPoints(
      __attribute__((unused)) const Array<Real> & field,
      __attribute__((unused)) Array<Real> & normal,
      __attribute__((unused)) const ElementType & type,
      __attribute__((unused)) const GhostType & ghost_type = _not_ghost) const {
    AKANTU_TO_IMPLEMENT();
  }

  /// function to print the containt of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;

private:
  /// initialise the class
  void init();

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  using ElementTypesIteratorHelper =
      ElementTypeMapArray<Real, ElementType>::ElementTypesIteratorHelper;

  ElementTypesIteratorHelper elementTypes(UInt dim = _all_dimensions,
                                          GhostType ghost_type = _not_ghost,
                                          ElementKind kind = _ek_regular) const;

  /// get the dimension of the element handeled by this fe_engine object
  AKANTU_GET_MACRO(ElementDimension, element_dimension, UInt);

  /// get the mesh contained in the fem object
  AKANTU_GET_MACRO(Mesh, mesh, const Mesh &);
  /// get the mesh contained in the fem object
  AKANTU_GET_MACRO_NOT_CONST(Mesh, mesh, Mesh &);

  /// get the in-radius of an element
  static inline Real getElementInradius(const Matrix<Real> & coord,
                                        const ElementType & type);

  /// get the normals on integration points
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(NormalsOnIntegrationPoints,
                                         normals_on_integration_points, Real);

  /// get cohesive element type for a given facet type
  static inline ElementType
  getCohesiveElementType(const ElementType & type_facet);

  /// get igfem element type for a given regular type
  static inline Vector<ElementType>
  getIGFEMElementTypes(const ElementType & type);

  /// get the interpolation element associated to an element type
  static inline InterpolationType
  getInterpolationType(const ElementType & el_type);

  /// get the shape function class (probably useless: see getShapeFunction in
  /// fe_engine_template.hh)
  virtual const ShapeFunctions & getShapeFunctionsInterface() const = 0;
  /// get the integrator class (probably useless: see getIntegrator in
  /// fe_engine_template.hh)
  virtual const Integrator & getIntegratorInterface() const = 0;

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// spatial dimension of the problem
  UInt element_dimension;

  /// the mesh on which all computation are made
  Mesh & mesh;

  /// normals at integration points
  ElementTypeMapArray<Real> normals_on_integration_points;
};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

/// standard output stream operator
inline std::ostream & operator<<(std::ostream & stream,
                                 const FEEngine & _this) {
  _this.printself(stream);
  return stream;
}

} // namespace akantu

#include "fe_engine_inline_impl.cc"
#include "fe_engine_template.hh"

#endif /* __AKANTU_FE_ENGINE_HH__ */

/**
 * @file   shape_structural.hh
 *
 * @author Fabian Barras <fabian.barras@epfl.ch>
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Tue Feb 15 2011
 * @date last modification: Tue Feb 20 2018
 *
 * @brief  shape class for element with different set of shapes functions
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
#include "shape_functions.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_SHAPE_STRUCTURAL_HH__
#define __AKANTU_SHAPE_STRUCTURAL_HH__

namespace akantu {

template <ElementKind kind> class ShapeStructural : public ShapeFunctions {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
  // Ctors/Dtors should be explicitely implemented for _ek_structural
public:
  ShapeStructural(Mesh & mesh, UInt spatial_dimension,
                  const ID & id = "shape_structural",
                  const MemoryID & memory_id = 0);
  ~ShapeStructural() override;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// function to print the contain of the class
  void printself(std::ostream & stream, int indent = 0) const override {
    std::string space;
    for (Int i = 0; i < indent; i++, space += AKANTU_INDENT)
      ;
    stream << space << "ShapesStructural [" << std::endl;
    rotation_matrices.printself(stream, indent + 1);
    ShapeFunctions::printself(stream, indent + 1);
    stream << space << "]" << std::endl;
  }

  /// compute shape functions on given integration points
  template <ElementType type>
  void computeShapesOnIntegrationPoints(
      const Array<Real> &, const Matrix<Real> & integration_points,
      Array<Real> & shapes, const GhostType & ghost_type,
      const Array<UInt> & filter_elements = empty_filter) const;

  /// initialization function for structural elements
  inline void initShapeFunctions(const Array<Real> & nodes,
                                 const Matrix<Real> & integration_points,
                                 const ElementType & type,
                                 const GhostType & ghost_type);

  /// precompute the rotation matrices for the elements dofs
  template <ElementType type>
  void precomputeRotationMatrices(const Array<Real> & nodes,
                                  const GhostType & ghost_type);

  /// pre compute all shapes on the element integration points from natural
  /// coordinates
  template <ElementType type>
  void precomputeShapesOnIntegrationPoints(const Array<Real> & nodes,
                                           const GhostType & ghost_type);

  /// pre compute all shapes on the element integration points from natural
  /// coordinates
  template <ElementType type>
  void
  precomputeShapeDerivativesOnIntegrationPoints(const Array<Real> & nodes,
                                                const GhostType & ghost_type);

  /// interpolate nodal values on the integration points
  template <ElementType type>
  void interpolateOnIntegrationPoints(
      const Array<Real> & u, Array<Real> & uq, UInt nb_degree_of_freedom,
      const GhostType & ghost_type = _not_ghost,
      const Array<UInt> & filter_elements = empty_filter) const;

  /// compute the gradient of u on the integration points
  template <ElementType type>
  void gradientOnIntegrationPoints(
      const Array<Real> & u, Array<Real> & nablauq, UInt nb_degree_of_freedom,
      const GhostType & ghost_type = _not_ghost,
      const Array<UInt> & filter_elements = empty_filter) const;

  /// interpolate on physical point
  template <ElementType type>
  void interpolate(const Vector<Real> & /*real_coords*/, UInt /*elem*/,
                   const Matrix<Real> & /*nodal_values*/,
                   Vector<Real> & /*interpolated*/,
                   const GhostType & /*ghost_type*/) const {
    AKANTU_TO_IMPLEMENT();
  }

  /// compute the shapes on a provided point
  template <ElementType type>
  void computeShapes(const Vector<Real> & /*real_coords*/, UInt /*elem*/,
                     Vector<Real> & /*shapes*/,
                     const GhostType & /*ghost_type*/) const {
    AKANTU_TO_IMPLEMENT();
  }

  /// compute the shape derivatives on a provided point
  template <ElementType type>
  void computeShapeDerivatives(const Matrix<Real> & /*real_coords*/,
                               UInt /*elem*/, Tensor3<Real> & /*shapes*/,
                               const GhostType & /*ghost_type*/) const {
    AKANTU_TO_IMPLEMENT();
  }

  /// get the rotations vector
  inline const Array<Real> &
  getRotations(const ElementType & el_type,
               __attribute__((unused))
               const GhostType & ghost_type = _not_ghost) const {
    return rotation_matrices(el_type);
  }

  /* ------------------------------------------------------------------------ */
  template <ElementType type>
  void computeBtD(const Array<Real> & /*Ds*/, Array<Real> & /*BtDs*/,
                  GhostType /*ghost_type*/,
                  const Array<UInt> & /*filter_elements*/) const {
    AKANTU_TO_IMPLEMENT();
  }

  template <ElementType type>
  void computeBtDB(const Array<Real> & /*Ds*/, Array<Real> & /*BtDBs*/,
                   UInt /*order_d*/, GhostType /*ghost_type*/,
                   const Array<UInt> & /*filter_elements*/) const {
    AKANTU_TO_IMPLEMENT();
  }

  /// multiply a field by shape functions
  template <ElementType type>
  void
  computeNtb(const Array<Real> & /*bs*/, Array<Real> & /*Ntbs*/,
             const GhostType & /*ghost_type*/,
             const Array<UInt> & /*filter_elements*/ = empty_filter) const {
    AKANTU_TO_IMPLEMENT();
  }

protected:
  ElementTypeMapArray<Real> rotation_matrices;
};

} // namespace akantu

#include "shape_structural_inline_impl.hh"

#endif /* __AKANTU_SHAPE_STRUCTURAL_HH__ */

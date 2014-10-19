/**
 * @file   shape_linked.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Fabian Barras <fabian.barras@epfl.ch>
 *
 * @date creation: Fri Jul 15 2011
 * @date last modification: Fri Jun 13 2014
 *
 * @brief  shape class for element with different set of shapes functions
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

#include "shape_functions.hh"

#ifndef __AKANTU_SHAPE_LINKED_HH__
#define __AKANTU_SHAPE_LINKED_HH__

__BEGIN_AKANTU__

template <ElementKind kind>
class ShapeLinked : public ShapeFunctions {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  typedef ElementTypeMap<Array<Real> **> ElementTypeMapMultiReal;

  ShapeLinked(Mesh & mesh, const ID & id = "shape_linked", const MemoryID & memory_id = 0);
  virtual ~ShapeLinked();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  inline void initShapeFunctions(const Array<Real> & nodes,
				 const Matrix<Real> & control_points,
				 const ElementType & type,
				 const GhostType & ghost_type);

  /// pre compute all shapes on the element control points from natural coordinates
  template <ElementType type>
  void precomputeShapesOnControlPoints(const Array<Real> & nodes,
				       const GhostType & ghost_type);

  /// pre compute all shapes on the element control points from natural coordinates
  template <ElementType type>
  void precomputeShapeDerivativesOnControlPoints(const Array<Real> & nodes,
						 const GhostType & ghost_type);

  /// interpolate nodal values on the control points
  template <ElementType type>
  void interpolateOnControlPoints(const Array<Real> &u,
				  Array<Real> &uq,
				  UInt nb_degree_of_freedom,
				  const GhostType & ghost_type = _not_ghost,
				  const Array<UInt> & filter_elements = empty_filter,
				  bool accumulate = false,
				  UInt id_shape = 0,
				  UInt num_degre_of_freedom_to_interpolate = 0,
				  UInt num_degre_of_freedom_interpolated = 0) const;


  /// compute the gradient of u on the control points
  template <ElementType type>
  void gradientOnControlPoints(const Array<Real> &u,
			       Array<Real> &nablauq,
			       UInt nb_degree_of_freedom,
			       const GhostType & ghost_type = _not_ghost,
			       const Array<UInt> & filter_elements = empty_filter,
			       bool accumulate = false,
			       UInt id_shape = 0,
			       UInt num_degre_of_freedom_to_interpolate = 0,
			       UInt num_degre_of_freedom_interpolated = 0) const;

  /// multiply a field by shape functions
  template <ElementType type>
  void fieldTimesShapes(__attribute__((unused)) const Array<Real> & field,
			__attribute__((unused)) Array<Real> & fiedl_times_shapes,
			__attribute__((unused)) const GhostType & ghost_type) const {
    AKANTU_DEBUG_TO_IMPLEMENT();
  }

private:
  template <ElementType type>
  void extractNodalToElementField(const Array<Real> & nodal_f,
				  Array<Real> & elemental_f,
				  UInt num_degre_of_freedom_to_extract,
				  const GhostType & ghost_type,
				  const Array<UInt> & filter_elements) const;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  /// get a the shapes vector
  inline const Array<Real> & getShapes(const ElementType & type,
					const GhostType & ghost_type,
					UInt id = 0) const;

  /// get a the shapes derivatives vector
  inline const Array<Real> & getShapesDerivatives(const ElementType & type,
						   const GhostType & ghost_type,
						   UInt id = 0) const;

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// shape functions for all elements
  ElementTypeMapMultiReal shapes;

  /// shape derivatives for all elements
  ElementTypeMapMultiReal shapes_derivatives;
};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#if defined (AKANTU_INCLUDE_INLINE_IMPL)
#  include "shape_linked_inline_impl.cc"
#endif

/// standard output stream operator
// inline std::ostream & operator <<(std::ostream & stream, const ShapeLinked & _this)
// {
//   _this.printself(stream);
//   return stream;
// }


__END_AKANTU__

#endif /* __AKANTU_SHAPE_LINKED_HH__ */

/**
 * @file   shape_igfem.hh
 * @author Aurelia Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @date   Sun May 25 10:34:41 2014
 *
 * @brief  shape functions for interface-enriched generalized FEM
 *
 * @section LICENSE
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

#ifndef __AKANTU_SHAPE_IGFEM_HH__
#define __AKANTU_SHAPE_IGFEM_HH__

#include "shape_functions.hh"
#include "igfem_element.hh"

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
template<>
class ShapeLagrange<_ek_igfem> : public ShapeFunctions{
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  ShapeLagrange(const Mesh & mesh,
		const ID & id = "shape_igfem",
		const MemoryID & memory_id = 0);
  virtual ~ShapeLagrange(){};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  inline void initShapeFunctions(const Array<Real> & nodes,
				 const Matrix<Real> & control_points,
				 const ElementType & type,
				 const GhostType & ghost_type);

  /// set the control points for a given element
  template <ElementType type, bool is_sub>
  void setControlPointsByType(const Matrix<Real> & control_points,
			      const GhostType & ghost_type);

  /// pre compute all shapes on the element control points from natural coordinates
  template<ElementType type, bool is_sub>
  void precomputeShapesOnControlPoints(const Array<Real> & nodes,
				       GhostType ghost_type);

  /// pre compute all shape derivatives on the element control points from natural coordinates
  template <ElementType type, bool is_sub>
  void precomputeShapeDerivativesOnControlPoints(const Array<Real> & nodes,
						 GhostType ghost_type);

  /// interpolate nodal values on the control points
  template <ElementType type>
  void interpolateOnControlPoints(const Array<Real> &u,
				  Array<Real> &uq,
				  UInt nb_degree_of_freedom,
				  GhostType ghost_type = _not_ghost,
				  const Array<UInt> & filter_elements = empty_filter) const;

  template <ElementType type, bool is_sub>
  void interpolateOnControlPoints(const Array<Real> &u,
				  Array<Real> &uq,
				  UInt nb_degree_of_freedom,
				  GhostType ghost_type = _not_ghost,
				  const Array<UInt> & filter_elements = empty_filter) const;

  /// compute the gradient of u on the control points
  template <ElementType type>
  void gradientOnControlPoints(const Array<Real> &u,
			       Array<Real> &nablauq,
			       UInt nb_degree_of_freedom,
			       GhostType ghost_type = _not_ghost,
			       const Array<UInt> & filter_elements = empty_filter) const;

  template <ElementType type, bool is_sub>
  void gradientOnControlPoints(const Array<Real> &u,
			       Array<Real> &nablauq,
			       UInt nb_degree_of_freedom,
			       GhostType ghost_type = _not_ghost,
			       const Array<UInt> & filter_elements = empty_filter) const;

  /// multiply a field by shape functions  @f$ fts_{ij} = f_i * \varphi_j @f$
  template <ElementType type>
  void fieldTimesShapes(const Array<Real> & field,
			Array<Real> & field_times_shapes,
			GhostType ghost_type) const;

  template <ElementType type, bool is_sub>
  void fieldTimesShapes(const Array<Real> & field,
			Array<Real> & field_times_shapes,
			GhostType ghost_type) const;

  /// find natural coords from real coords provided an element
  template <ElementType type>
  void inverseMap(const Vector<Real> & real_coords,
		  UInt element,
		  Vector<Real> & natural_coords,
		  const GhostType & ghost_type = _not_ghost) const;

  template <ElementType type, bool is_sub>
  void inverseMap(const Vector<Real> & real_coords,
		  UInt element,
		  Vector<Real> & natural_coords,
		  const GhostType & ghost_type = _not_ghost) const;

  /// return true if the coordinates provided are inside the element, false otherwise
  template <ElementType type>
  bool contains(const Vector<Real> & real_coords,
		UInt elem,
		const GhostType & ghost_type) const;

  template <ElementType type, bool is_sub>
  bool contains(const Vector<Real> & real_coords,
		UInt elem,
		const GhostType & ghost_type) const;

  /// compute the shape on a provided point
  template <ElementType type>
  void computeShapes(const Vector<Real> & real_coords,
		     UInt elem,
		     Vector<Real> & shapes,
		     const GhostType & ghost_type) const;

  template <ElementType type, bool is_sub>
  void computeShapes(const Vector<Real> & real_coords,
		     UInt elem,
		     Vector<Real> & shapes,
		     const GhostType & ghost_type) const;

  /// function to print the containt of the class
  virtual void printself(std::ostream & stream, int indent = 0) const;

  template<bool is_sub>
  void printself(std::ostream & stream, int indent = 0) const;

protected:
  /// compute the shape derivatives on control points for a given element
  template <ElementType type, bool is_sub>
  inline void computeShapeDerivativesOnCPointsByElement(const Matrix<Real> & node_coords,
							const Matrix<Real> & natural_coords,
							Tensor3<Real> & shapesd);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// get a the shapes vector
  inline const Array<Real> & getShapes(const ElementType & el_type,
					const GhostType & ghost_type = _not_ghost) const;

  /// get a the shapes derivatives vector
  inline const Array<Real> & getShapesDerivatives(const ElementType & el_type,
						   const GhostType & ghost_type = _not_ghost) const;

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// shape functions for all elements
  ElementTypeMapArray<Real, InterpolationType> shapes;

  /// shape functions derivatives for all elements
  ElementTypeMapArray<Real, InterpolationType> shapes_derivatives;

  /// for parent interpolation types store boolean to decide wheter element is enriched or not
  ElementTypeMapArray<bool> is_split;

  /// for subelement store the control points of the corresponding parent elements
  ElementTypeMapArray<Real, InterpolationType> parent_control_points;
  
};


/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */
#include "shape_igfem_inline_impl.cc"

__END_AKANTU__

#endif /* __AKANTU_SHAPE_IGFEM_HH__ */

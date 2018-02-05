/**
 * @file   material_reinforcement.hh
 *
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 *
 * @date creation: Fri Mar 13 2015
 * @date last modification: Tue Nov 24 2015
 *
 * @brief  Reinforcement material
 *
 * @section LICENSE
 *
 * Copyright (©) 2015 EPFL (Ecole Polytechnique Fédérale de Lausanne) Laboratory
 * (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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

#ifndef __AKANTU_MATERIAL_REINFORCEMENT_HH__
#define __AKANTU_MATERIAL_REINFORCEMENT_HH__

#include "aka_common.hh"

#include "material.hh"
#include "embedded_interface_model.hh"
#include "embedded_internal_field.hh"

/* -------------------------------------------------------------------------- */

namespace akantu {

/**
 * @brief Material used to represent embedded reinforcements
 *
 * This class is used for computing the reinforcement stiffness matrix
 * along with the reinforcement residual. Room is made for constitutive law,
 * but actual use of contitutive laws is made in MaterialReinforcementTemplate.
 *
 * Be careful with the dimensions in this class :
 *  -  this->spatial_dimension is always 1
 *  -  the template parameter dim is the dimension of the problem
 */

template <class Mat, UInt dim> class MaterialReinforcement : public Mat {

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  /// Constructor
  MaterialReinforcement(EmbeddedInterfaceModel & model, const ID & id = "");

  /// Destructor
  ~MaterialReinforcement() override;

protected:
  void initialize();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// Init the material
  void initMaterial() override;

  /// Init the filter for background elements
  void initBackgroundFilter();

  /// Init the background shape derivatives
  void initBackgroundShapeDerivatives();

  /// Init the cosine matrices
  void initDirectingCosines();

  /// Assemble stiffness matrix
  void assembleStiffnessMatrix(GhostType ghost_type) override;

  /// Compute all the stresses !
  void computeAllStresses(GhostType ghost_type) override;

  /// Compute energy
  Real getEnergy(const std::string & id) override;

  /// Assemble the residual of one type of element (typically _segment_2)
  void assembleInternalForces(GhostType ghost_type) override;

  // virtual ElementTypeMap<UInt> getInternalDataPerElem(const ID & field_name,
  //                                                     const ElementKind &
  //                                                     kind,
  //                                                     const ID &
  //                                                     fe_engine_id) const;

  // /// Reimplementation of Material's function to accomodate for interface
  // mesh
  // virtual void flattenInternal(const std::string & field_id,
  //                              ElementTypeMapArray<Real> & internal_flat,
  //                              const GhostType ghost_type = _not_ghost,
  //                              ElementKind element_kind = _ek_not_defined)
  //                              const;

  /* ------------------------------------------------------------------------ */
  /* Protected methods                                                        */
  /* ------------------------------------------------------------------------ */
protected:
  /// Allocate the background shape derivatives
  void allocBackgroundShapeDerivatives();

  /// Compute the directing cosines matrix for one element type
  void computeDirectingCosines(const ElementType & type, GhostType ghost_type);

  /// Compute the directing cosines matrix on quadrature points.
  inline void computeDirectingCosinesOnQuad(const Matrix<Real> & nodes,
                                            Matrix<Real> & cosines);

  /// Assemble the stiffness matrix for an element type (typically _segment_2)
  void assembleStiffnessMatrix(const ElementType & type, GhostType ghost_type);

  /// Assemble the stiffness matrix for background & interface types
  void assembleStiffnessMatrixInterface(const ElementType & interface_type,
                                        const ElementType & background_type,
                                        GhostType ghost_type);

  /// Compute the background shape derivatives for a type
  void computeBackgroundShapeDerivatives(const ElementType & type,
                                         GhostType ghost_type);

  /// Compute the background shape derivatives for a type pair
  void computeBackgroundShapeDerivatives(const ElementType & interface_type,
					 const ElementType & bg_type,
					 GhostType ghost_type,
					 const Array<UInt> & filter);

  /// Filter elements crossed by interface of a type
  void filterInterfaceBackgroundElements(Array<UInt> & filter,
                                         const ElementType & type,
                                         const ElementType & interface_type,
                                         GhostType ghost_type);

  /// Assemble the residual of one type of element (typically _segment_2)
  void assembleInternalForces(const ElementType & type, GhostType ghost_type);

  /// Assemble the residual for a pair of elements
  void assembleInternalForcesInterface(const ElementType & interface_type,
                                       const ElementType & background_type,
                                       GhostType ghost_type);

  // TODO figure out why voigt size is 4 in 2D
  inline void stressTensorToVoigtVector(const Matrix<Real> & tensor,
                                        Vector<Real> & vector);
  inline void strainTensorToVoigtVector(const Matrix<Real> & tensor,
                                        Vector<Real> & vector);

  /// Compute gradu on the interface quadrature points
  void computeGradU(const ElementType & type, GhostType ghost_type);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// Embedded model
  EmbeddedInterfaceModel & emodel;

  /// Stress in the reinforcement
  InternalField<Real> stress_embedded;

  /// Gradu of concrete on reinforcement
  InternalField<Real> gradu_embedded;

  /// C matrix on quad
  InternalField<Real> directing_cosines;

  /// Prestress on quad
  InternalField<Real> pre_stress;

  /// Cross-sectional area
  Real area;

  template <typename T>
  using CrossMap = ElementTypeMap<std::unique_ptr<ElementTypeMapArray<T>>>;

  /// Background mesh shape derivatives
  CrossMap<Real> shape_derivatives;

  /// Background element filter (contains segment-bg pairs)
  CrossMap<UInt> background_filter;
};

} // akantu

#include "material_reinforcement_tmpl.hh"

#endif // __AKANTU_MATERIAL_REINFORCEMENT_HH__

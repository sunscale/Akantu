/**
 * @file   material_reinforcement.hh
 *
 * @author Lucas Frérot <lucas.frerot@epfl.ch>
 *
 * @date creation: Thu Mar 12 2015
 * @date last modification: Thu Mar 12 2015
 *
 * @brief  Reinforcement material
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

#ifndef __AKANTU_MATERIAL_REINFORCEMENT_HH__
#define __AKANTU_MATERIAL_REINFORCEMENT_HH__

#include "aka_common.hh"

#include "material.hh"
#include "embedded_interface_model.hh"

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

template<UInt d>
class MaterialReinforcement : public Material {

public:
  /// Constructor
  MaterialReinforcement(SolidMechanicsModel & model, const ID & id = "");

  /// Destructor
  virtual ~MaterialReinforcement();

public:
  /// Init the background shape derivatives
  virtual void initBackgroundShapeDerivatives();

  /// Assemble stiffness matrix
  virtual void assembleStiffnessMatrix(GhostType ghost_type);

  /// Compute the stiffness parameter for elements of a type
  virtual void computeStiffness(const ElementType & type, GhostType ghost_type);

protected:
  /// Compute the directing cosines matrix for one element type
  void computeDirectingCosines(const ElementType & type, GhostType ghost_type);

  /// Compute the directing cosines matrix on quadrature points
  inline void computeDirectingCosinesOnQuad(const Matrix<Real> & nodes,
                                            Matrix<Real> & cosines);

  /// Assemble the stiffness matrix for an element type (typically _segment_2)
  void assembleStiffnessMatrix(const ElementType & type, GhostType ghost_type);

  /// Compute the background shape derivatives for a type (typically _triangle_3 / _tetrahedron_4)
  void computeBackgroundShapeDerivatives(const ElementType & type, GhostType ghost_type);

  /// Filter elements crossed by interface of a type
  void filterInterfaceBackgroundElements(Array<UInt> & filter,
                                         const ElementType & type,
                                         const ElementType & interface_type,
                                         GhostType ghost_type,
                                         GhostType interface_ghost_type);

protected:
  /// Embedded model
  EmbeddedInterfaceModel * model;

  /// C matrix on quad
  InternalField<Real> directing_cosines;

  /// D on quad
  InternalField<Real> reinforcement_stiffness;

  /// Cross-sectional area
  InternalField<Real> area;

  /// Background mesh shape derivatives
  ElementTypeMap< ElementTypeMapArray<Real> > shape_derivatives;

};

#include "material_reinforcement_inline_impl.cc"

__END_AKANTU__

#endif // __AKANTU_MATERIAL_REINFORCEMENT_HH__

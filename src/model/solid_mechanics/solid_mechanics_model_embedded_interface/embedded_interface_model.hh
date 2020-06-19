/**
 * @file   embedded_interface_model.hh
 *
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Wed Jan 31 2018
 *
 * @brief  Model of Solid Mechanics with embedded interfaces
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

#ifndef __AKANTU_EMBEDDED_INTERFACE_MODEL_HH__
#define __AKANTU_EMBEDDED_INTERFACE_MODEL_HH__

#include "aka_common.hh"

#include "mesh.hh"
#include "solid_mechanics_model.hh"

#include "embedded_interface_intersector.hh"

/* -------------------------------------------------------------------------- */

namespace akantu {

/**
 * @brief Solid mechanics model using the embedded model.
 *
 * This SolidMechanicsModel subclass implements the embedded model,
 * a method used to represent 1D elements in a finite elements model
 * (eg. reinforcements in concrete).
 *
 * In addition to the SolidMechanicsModel properties, this model has
 * a mesh of the 1D elements embedded in the model, and an instance of the
 * EmbeddedInterfaceIntersector class for the computation of the intersections
 * of the
 * 1D elements with the background (bulk) mesh.
 *
 * @see MaterialReinforcement
 */
class EmbeddedInterfaceModel : public SolidMechanicsModel {

  using MyFEEngineType = SolidMechanicsModel::MyFEEngineType;

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  /**
   * @brief Constructor
   *
   * @param mesh main mesh (concrete)
   * @param primitive_mesh mesh of the embedded reinforcement
   * @param spatial_dimension the spatial dimension to be considered by this model
   * @param id the id of the model
   * @param memory_id the id of the memory manager to use
   */
  EmbeddedInterfaceModel(Mesh & mesh, Mesh & primitive_mesh,
                         UInt spatial_dimension = _all_dimensions,
                         const ID & id = "embedded_interface_model",
                         const MemoryID & memory_id = 0);

  /// Destructor
  ~EmbeddedInterfaceModel() override;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// Initialise the model
  void initFullImpl(
      const ModelOptions & options = EmbeddedInterfaceModelOptions()) override;

  /// Initialise the materials
  void
  assignMaterialToElements(const ElementTypeMapArray<UInt> * filter) override;

  /// Initialize the embedded shape functions
  void initModel() override;

  /// Allows filtering of dump fields which need to be dumpes on interface mesh
  void addDumpGroupFieldToDumper(const std::string & dumper_name,
                                 const std::string & field_id,
                                 const std::string & group_name,
                                 const ElementKind & element_kind,
                                 bool padding_flag) override;

  // virtual ElementTypeMap<UInt> getInternalDataPerElem(const std::string &
  // field_name,
  //                                                     const ElementKind &
  //                                                     kind);
  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// Get interface mesh
  AKANTU_GET_MACRO(InterfaceMesh, *interface_mesh, Mesh &);

  /// Get associated elements
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE(
      InterfaceAssociatedElements,
      interface_mesh->getData<Element>("associated_element"), Element);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// Intersector object to build the interface mesh
  EmbeddedInterfaceIntersector intersector;

  /// Interface mesh (weak reference)
  Mesh * interface_mesh;

  /// Mesh used to create the CGAL primitives for intersections
  Mesh & primitive_mesh;

  /// Material selector for interface
  MaterialSelector * interface_material_selector;
};

/// Material selector based on mesh data for interface elements
template <typename T>
class InterfaceMeshDataMaterialSelector
    : public ElementDataMaterialSelector<T> {
public:
  InterfaceMeshDataMaterialSelector(const std::string & name,
                                    const EmbeddedInterfaceModel & model,
                                    UInt first_index = 1)
      : ElementDataMaterialSelector<T>(
            model.getInterfaceMesh().getData<T>(name), model, first_index) {}
};

} // namespace akantu

#endif // __AKANTU_EMBEDDED_INTERFACE_MODEL_HH__

/**
 * @file   embedded_interface_model.hh
 *
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 *
 * @date creation: Fri Jun 18 2010
 * @date last modification: Mon Dec 14 2015
 *
 * @brief  Model of Solid Mechanics with embedded interfaces
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

#ifndef __AKANTU_EMBEDDED_INTERFACE_MODEL_HH__
#define __AKANTU_EMBEDDED_INTERFACE_MODEL_HH__

#include "aka_common.hh"

#include "solid_mechanics_model.hh"
#include "mesh.hh"

#include "embedded_interface_intersector.hh"

/* -------------------------------------------------------------------------- */

namespace akantu {

/// Options for the EmbeddedInterfaceModel
struct EmbeddedInterfaceModelOptions : SolidMechanicsModelOptions {
  /**
   * @brief Constructor for EmbeddedInterfaceModelOptions
   * @param analysis_method see SolidMeechanicsModelOptions
   * @param no_init_intersections ignores the embedded elements
   * @param no_init_materials see SolidMeechanicsModelOptions
   */
  EmbeddedInterfaceModelOptions(AnalysisMethod analysis_method = _explicit_lumped_mass,
                                bool no_init_intersections = false,
                                bool no_init_materials = false):
    SolidMechanicsModelOptions(analysis_method, no_init_materials),
    no_init_intersections(no_init_intersections)
  {}

  /// Ignore the reinforcements
  bool no_init_intersections;
};

extern const EmbeddedInterfaceModelOptions default_embedded_interface_model_options;

/**
 * @brief Solid mechanics model using the embedded model.
 *
 * This SolidMechanicsModel subclass implements the embedded model,
 * a method used to represent 1D elements in a finite elements model
 * (eg. reinforcements in concrete).
 *
 * In addition to the SolidMechanicsModel properties, this model has
 * a mesh of the 1D elements embedded in the model, and an instance of the
 * EmbeddedInterfaceIntersector class for the computation of the intersections of the
 * 1D elements with the background (bulk) mesh.
 *
 * @see MaterialReinforcement
 */
class EmbeddedInterfaceModel : public SolidMechanicsModel {

  /// Same type as SolidMechanicsModel
  typedef FEEngineTemplate<IntegratorGauss, ShapeLagrange, _ek_regular> MyFEEngineType;


  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  /**
   * @brief Constructor
   *
   * @param mesh main mesh (concrete)
   * @param primitive_mesh mesh of the embedded reinforcement
   */
  EmbeddedInterfaceModel(Mesh & mesh,
                         Mesh & primitive_mesh,
                         UInt spatial_dimension = _all_dimensions,
                         const ID & id = "embedded_interface_model",
                         const MemoryID & memory_id = 0);

  /// Destructor
  virtual ~EmbeddedInterfaceModel();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// Initialise the model
  virtual void initFull(const ModelOptions & options = default_embedded_interface_model_options);

  /// Initialise the materials
  virtual void initMaterials();

  /// Allows filtering of dump fields which need to be dumpes on interface mesh
  virtual void addDumpGroupFieldToDumper(const std::string & dumper_name,
                                         const std::string & field_id,
                                         const std::string & group_name,
                                         const ElementKind & element_kind,
                                         bool padding_flag);

  // virtual ElementTypeMap<UInt> getInternalDataPerElem(const std::string & field_name,
  //                                                     const ElementKind & kind);
  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// Get interface mesh
  AKANTU_GET_MACRO(InterfaceMesh, *interface_mesh, Mesh &);

  /// Get associated elements
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE(InterfaceAssociatedElements,
                                   interface_mesh->getData<Element>("associated_element"),
                                   Element);

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
template<typename T>
class InterfaceMeshDataMaterialSelector : public ElementDataMaterialSelector<T> {
public:
  InterfaceMeshDataMaterialSelector(const std::string & name, const EmbeddedInterfaceModel & model, UInt first_index = 1) :
    ElementDataMaterialSelector<T>(model.getInterfaceMesh().getData<T>(name), model, first_index)
  {}
};

} // akantu

#endif // __AKANTU_EMBEDDED_INTERFACE_MODEL_HH__

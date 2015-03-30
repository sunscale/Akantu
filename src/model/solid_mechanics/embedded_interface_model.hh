/**
 * @file   embedded_interface_model.hh
 *
 * @author Lucas Frérot <lucas.frerot@epfl.ch>
 *
 * @date creation: Mon Mar 9 2015
 * @date last modification: Mon Mar 9 2015
 *
 * @brief  Model of Solid Mechanics with embedded interfaces
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

#ifndef __AKANTU_EMBEDDED_INTERFACE_MODEL_HH__
#define __AKANTU_EMBEDDED_INTERFACE_MODEL_HH__

#include "aka_common.hh"

#include "solid_mechanics_model.hh"
#include "mesh.hh"
#include "parsable.hh"

#include "embedded_interface.hh"
#include "mesh_geom_container.hh"

#include <CGAL/Cartesian.h>

/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__

typedef CGAL::Cartesian<Real> K;

class EmbeddedInterfaceModel : public SolidMechanicsModel {

  typedef FEEngineTemplate<IntegratorGauss, ShapeLagrange, _ek_regular> MyFEEngineType;


  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  /// Constructor
  EmbeddedInterfaceModel(Mesh & mesh,
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
  virtual void initModel();

  /// Initialise the materials
  virtual void initMaterials();

  /// Instanciate the interfaces
  virtual void instanciateInterfaces();

  /// Initialise the interface mesh
  void initInterface(const std::list<std::pair<K::Segment_3, std::string> > & interface_list);

  /// Assemble the stiffness matrix of the model
  virtual void assembleStiffnessMatrix();

  /// Dump
  /*virtual ElementTypeMap<UInt> getInternalDataPerElem(const std::string & field_name,
                                                      const ElementKind & ek_kind);*/

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// get interface mesh
  AKANTU_GET_MACRO(InterfaceMesh, *interface_mesh, Mesh &);

  /// get associated elements
  AKANTU_GET_MACRO_BY_ELEMENT_TYPE(InterfaceAssociatedElements,
                                   interface_mesh->getData<Element>("associated_element"),
                                   Element);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// Interface mesh (weak reference)
  Mesh * interface_mesh;

  /// Material selector for interface
  MaterialSelector * interface_material_selector;
  
  /// Geom object to build the interface mesh
  MeshGeomContainer interface_container;

};

/// Material selector for EmbeddedInterfaceModel
template <typename T> class EmbeddedInterfaceMaterialSelector : public DefaultMaterialSelector {};

/// Material selector based on associated material name
template<>
class EmbeddedInterfaceMaterialSelector<std::string> : public DefaultMaterialSelector {

public:
  EmbeddedInterfaceMaterialSelector(const std::string & name, const EmbeddedInterfaceModel & model):
    DefaultMaterialSelector(model.getMaterialByElement()),
    names(model.getInterfaceMesh().getData<std::string>(name)),
    model(model)
  {}

  UInt operator() (const Element & element) {
    try {
      DebugLevel dbl = debug::getDebugLevel();
      debug::setDebugLevel(dblError);

      std::string material_name = names(element.type, element.ghost_type)(element.element);

      debug::setDebugLevel(dbl);

      return model.getMaterialIndex(material_name);
    } catch (debug::Exception& e) {
      /// TODO this returns the last material, which should be a bulk material. Try to find cleaner way
      return model.getNbMaterials() - 1;
    }
  }

protected:
  const ElementTypeMapArray<std::string> & names;
  const EmbeddedInterfaceModel & model;
};

__END_AKANTU__

#endif // __AKANTU_EMBEDDED_INTERFACE_MODEL_HH__

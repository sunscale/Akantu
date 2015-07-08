/**
 * @file   element_class_igfem.hh
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 *
 * @brief  Parent material for IGFEM
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "material.hh"
#include "igfem_internal_field.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MATERIAL_IGFEM_HH__
#define __AKANTU_MATERIAL_IGFEM_HH__

/* -------------------------------------------------------------------------- */
namespace akantu {
  class SolidMechanicsModelIGFEM;
}

__BEGIN_AKANTU__


class MaterialIGFEM : public virtual Material {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  typedef FEEngineTemplate<IntegratorGauss,
		      ShapeLagrange, _ek_igfem> MyFEEngineIGFEMType;

  MaterialIGFEM(SolidMechanicsModel & model, const ID & id = "");

  virtual ~MaterialIGFEM();

protected:
 
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  inline void interpolateInternal(const ElementType & type,
				  const Vector<Real> & internal,
				  Vector<Real> & interpolated,
				  const UInt nb_quads,
				  const UInt sub_element);

  /* ------------------------------------------------------------------------ */
  /* MeshEventHandler inherited members                                       */
  /* ------------------------------------------------------------------------ */
public:
  /* ------------------------------------------------------------------------ */
  virtual void computeQuadraturePointsCoordinates(ElementTypeMapArray<Real> & quadrature_points_coordinates,
						  const GhostType & ghost_type) const;
  // virtual void onElementsAdded(const Array<Element> & element_list,
  //                              const NewElementsEvent & event) {};

  // virtual void onElementsRemoved(const Array<Element> & element_list,
  //                                const ElementTypeMapArray<UInt> & new_numbering,
  //                                const RemovedElementsEvent & event) {};
protected:
  void initialize();

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */

protected:
  /// Link to the igfem fem object in the model
  MyFEEngineIGFEMType * fem_igfem;

  const UInt nb_sub_materials = 2;

  /// pointer to the solid mechanics model for igfem elements
  SolidMechanicsModelIGFEM * model;

  ///internal field of bool to know to which sub-material a quad point belongs
  IGFEMInternalField<UInt> sub_material; 

  /// material name of first sub-material
  std::string name_sub_mat_1;

  /// material name of first sub-material
  std::string name_sub_mat_2;
};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "material_igfem_inline_impl.cc"



__END_AKANTU__

#include "igfem_internal_field_tmpl.hh"

#endif /* __AKANTU_MATERIAL_IGFEM_HH__ */

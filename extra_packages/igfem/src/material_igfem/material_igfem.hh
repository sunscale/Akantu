/**
 * @file   element_class_igfem.hh
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 *
 * @brief  Parent material for IGFEM
 *
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "igfem_internal_field.hh"
#include "material.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_MATERIAL_IGFEM_HH_
#define AKANTU_MATERIAL_IGFEM_HH_

/* -------------------------------------------------------------------------- */
namespace akantu {
class SolidMechanicsModelIGFEM;
}

namespace akantu {

class MaterialIGFEM : public virtual Material {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  typedef FEEngineTemplate<IntegratorGauss, ShapeLagrange, _ek_igfem>
      MyFEEngineIGFEMType;

  MaterialIGFEM(SolidMechanicsModel & model, const ID & id = "");

  virtual ~MaterialIGFEM();

protected:
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  virtual void computeAllStresses(GhostType ghost_type = _not_ghost);

  virtual void extrapolateInternal(const ID & id, const Element & element,
                                   const Matrix<Real> & point,
                                   Matrix<Real> & extrapolated);

  /// apply a constant eigengrad_u everywhere in the material
  virtual void applyEigenGradU(const Matrix<Real> & prescribed_eigen_grad_u,
                               const ID & sub_mat_name,
                               const GhostType = _not_ghost);

  /* ------------------------------------------------------------------------ */
  /* MeshEventHandler inherited members                                       */
  /* ------------------------------------------------------------------------ */
public:
  /* ------------------------------------------------------------------------ */
  virtual void computeQuadraturePointsCoordinates(
      ElementTypeMapArray<Real> & quadrature_points_coordinates,
      GhostType ghost_type) const;
  // virtual void onElementsAdded(const Array<Element> & element_list,
  //                              const NewElementsEvent & event) {};

  // virtual void onElementsRemoved(const Array<Element> & element_list,
  //                                const ElementTypeMapArray<UInt> &
  //                                new_numbering,
  //                                const RemovedElementsEvent & event) {};
protected:
  /// constitutive law
  virtual void computeStress(__attribute__((unused)) ElementType el_type,
                             __attribute__((unused))
                             GhostType ghost_type = _not_ghost) {}
  void initialize();

  template <ElementType type>
  void setSubMaterial(const Array<Element> & element_list,
                      GhostType ghost_type);

  /* ------------------------------------------------------------------------ */
  /* DataAccessor inherited members                                           */
  /* ------------------------------------------------------------------------ */
public:
  virtual inline UInt getNbDataForElements(const Array<Element> & elements,
                                           SynchronizationTag tag) const;

  virtual inline void packElementData(CommunicationBuffer & buffer,
                                      const Array<Element> & elements,
                                      SynchronizationTag tag) const;

  virtual inline void unpackElementData(CommunicationBuffer & buffer,
                                        const Array<Element> & elements,
                                        SynchronizationTag tag);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */

protected:
  const UInt nb_sub_materials;

  /// pointer to the solid mechanics model for igfem elements
  SolidMechanicsModelIGFEM * model;

  /// internal field of bool to know to which sub-material a quad point belongs
  IGFEMInternalField<UInt> sub_material;

  /// material name of first sub-material
  std::string name_sub_mat_1;

  /// material name of first sub-material
  std::string name_sub_mat_2;

  /// map the index of the sub-materials to the names
  std::map<UInt, ID> sub_material_names;
};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "material_igfem_inline_impl.hh"

} // namespace akantu

#include "igfem_internal_field_tmpl.hh"

#endif /* AKANTU_MATERIAL_IGFEM_HH_ */

/**
 * @file   material_igfem_saw_tooth_damage.hh
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 *
 *
 * @brief  Linear saw-tooth softening material model for IGFEM elements
 *
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

/* -------------------------------------------------------------------------- */
#include "material_damage.hh"
#include "material_igfem_elastic.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_MATERIAL_IGFEM_SAW_TOOTH_DAMAGE_HH_
#define AKANTU_MATERIAL_IGFEM_SAW_TOOTH_DAMAGE_HH_

namespace akantu {

template <UInt spatial_dimension>
class MaterialIGFEMSawToothDamage
    : public MaterialDamage<spatial_dimension, MaterialIGFEMElastic> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
private:
  typedef MaterialDamage<spatial_dimension, MaterialIGFEMElastic> Parent;

public:
  typedef std::pair<Element, Element> ElementPair;
  MaterialIGFEMSawToothDamage(SolidMechanicsModel & model, const ID & id = "");
  MaterialIGFEMSawToothDamage(SolidMechanicsModel & model, UInt dim,
                              const Mesh & mesh, FEEngine & fe_engine,
                              const ID & id = "");

  virtual ~MaterialIGFEMSawToothDamage() {}

protected:
  void initialize();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  virtual void initMaterial();

  ///  virtual void updateInternalParameters();

  virtual void computeAllStresses(GhostType ghost_type = _not_ghost);

  /// update internal field damage
  virtual UInt updateDamage();

  UInt updateDamage(UInt quad_index, const Real eq_stress,
                    ElementType el_type, GhostType ghost_type);

  /// update energies after damage has been updated
  //  virtual void updateEnergiesAfterDamage(ElementType el_type, GhostType
  //  ghost_typ);

  virtual void onBeginningSolveStep(const AnalysisMethod & method){};

  virtual void onEndSolveStep(const AnalysisMethod & method){};

protected:
  /// constitutive law for all element of a type
  virtual void computeStress(ElementType el_type,
                             GhostType ghost_type = _not_ghost);

  /// compute the equivalent stress on each Gauss point (i.e. the max prinicpal
  /// stress) and normalize it by the tensile strength
  virtual void
  computeNormalizedEquivalentStress(const Array<Real> & grad_u,
                                    ElementType el_type,
                                    GhostType ghost_type = _not_ghost);

  /// find max normalized equivalent stress
  void findMaxNormalizedEquivalentStress(ElementType el_type,
                                         GhostType ghost_type = _not_ghost);

  inline void computeDamageAndStressOnQuad(Matrix<Real> & sigma, Real & dam);

protected:
  /* ------------------------------------------------------------------------ */
  /* MeshEventHandler inherited members                                       */
  /* ------------------------------------------------------------------------ */
public:
  /* ------------------------------------------------------------------------ */
  virtual void onElementsAdded(const Array<Element> & element_list,
                               const NewElementsEvent & event);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /// get max normalized equivalent stress
  AKANTU_GET_MACRO(NormMaxEquivalentStress, norm_max_equivalent_stress, Real);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// resistance to damage
  IGFEMInternalField<Real> Sc;

  /// internal field to store equivalent stress on each Gauss point
  IGFEMInternalField<Real> equivalent_stress;

  /// damage increment
  Real prescribed_dam;

  /// maximum equivalent stress
  Real norm_max_equivalent_stress;

  /// deviation from max stress at which Gauss point will still get damaged
  Real dam_tolerance;

  /// define damage threshold at which damage will be set to 1
  Real dam_threshold;

  /// maximum damage value
  Real max_damage;
};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "material_igfem_saw_tooth_damage_inline_impl.hh"

} // namespace akantu

#endif /* AKANTU_MATERIAL_IGFEM_SAW_TOOTH_DAMAGE_HH_ */

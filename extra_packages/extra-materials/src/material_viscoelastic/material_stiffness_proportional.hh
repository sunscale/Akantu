/**
 * @file   material_stiffness_proportional.hh
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 *
 * @brief  Material isotropic visco-elastic with viscosity proportional to the
 * stiffness
 *
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "material.hh"
#include "material_elastic.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MATERIAL_STIFFNESS_PROPORTIONAL_HH__
#define __AKANTU_MATERIAL_STIFFNESS_PROPORTIONAL_HH__

namespace akantu {

/**
 * Material visco-elastic @f[\sigma = E\epsilon + \alpha E*
 * \frac{d\epsilon}{dt}@f]
 * it can be seen as a Kelvin-Voigt solid with @f[\eta = \alpha E @f]
 *
 * The material satisfies the Caughey condition, the visco-elastic solid has the
 * same eigen-modes as the elastic one. (T.K. Caughey 1960 - Journal of Applied
 * Mechanics 27, 269-271. Classical normal modes in damped linear systems.)
 *
 * parameters in the material files :
 *   - rho : density (default: 0)
 *   - E   : Young's modulus (default: 0)
 *   - nu  : Poisson's ratio (default: 1/2)
 *   - Plane_Stress : if 0: plane strain, else: plane stress (default: 0)
 *   - alpha : viscous ratio
 */
template <UInt spatial_dimension>
class MaterialStiffnessProportional
    : public MaterialElastic<spatial_dimension> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  MaterialStiffnessProportional(SolidMechanicsModel & model,
                                const ID & id = "");

  virtual ~MaterialStiffnessProportional(){};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  void initMaterial();

  /// constitutive law for all element of a type
  void computeStress(ElementType el_type, GhostType ghost_type = _not_ghost);

  /// compute the potential energy for all elements
  virtual void computePotentialEnergy(ElementType el_type,
                                      GhostType ghost_type = _not_ghost);

protected:
  /// constitutive law for a given quadrature point
  // inline void computeStress(Real * F, Real * sigma);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// stress due to viscosity
  InternalField<Real> stress_viscosity;

  /// stress due to elasticity
  InternalField<Real> stress_elastic;

  /// viscous ratio
  Real alpha;
};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

//#include "material_elastic_caughey_inline_impl.hh"

} // namespace akantu

#endif /* __AKANTU_MATERIAL_STIFFNESS_PROPORTIONAL_HH__ */

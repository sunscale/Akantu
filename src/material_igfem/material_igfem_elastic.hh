/**
 * @file   element_class_igfem.hh
 *
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 *
 * @brief  Material isotropic elastic for IGFEM
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "material_igfem.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MATERIAL_IGFEM_ELASTIC_HH__
#define __AKANTU_MATERIAL_IGFEM_ELASTIC_HH__

__BEGIN_AKANTU__

/**
 * Material elastic isotropic
 *
 * parameters in the material files :
 *   - E   : Young's modulus (default: 0)
 *   - nu  : Poisson's ratio (default: 1/2)
 *   - Plane_Stress : if 0: plane strain, else: plane stress (default: 0)
 */
template<UInt spatial_dimension>
class MaterialIGFEMElastic : public MaterialIGFEM {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
private:
  typedef  Material Parent;
public:

  MaterialIGFEMElastic(SolidMechanicsModel & model, const ID & id = "");
  MaterialIGFEMElastic(SolidMechanicsModel & model,
		       UInt dim,
		       const Mesh & mesh,
		       FEEngine & fe_engine,
		       const ID & id = "");

  virtual ~MaterialIGFEMElastic() {}

protected:
  void initialize();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  virtual void initMaterial();
  /// constitutive law for all element of a type
  virtual void computeStress(ElementType el_type, GhostType ghost_type = _not_ghost);

  /// compute the tangent stiffness matrix for an element type
  virtual void computeTangentModuli(const ElementType & el_type,
				    Array<Real> & tangent_matrix,
				    GhostType ghost_type = _not_ghost);

  /// compute the elastic potential energy
  virtual void computePotentialEnergy(ElementType el_type,
				      GhostType ghost_type = _not_ghost);

  virtual void computePotentialEnergyByElement(ElementType type, UInt index,
					       Vector<Real> & epot_on_quad_points);


protected:
  /// update the default internal parameters
  void updateDefaultInternals(const UInt default_idx = 0);

  /// constitutive law for a given quadrature point
  inline void computeStressOnQuad(const Matrix<Real> & grad_u,
				  Matrix<Real> & sigma,
				  __attribute__((unused)) const Real lambda,
				  const Real mu) const;

  /// compute the tangent stiffness matrix for an element
  inline void computeTangentModuliOnQuad(Matrix<Real> & tangent,
					 __attribute__((unused)) const Real lambda,
					 const Real mu) const;

  static inline void computePotentialEnergyOnQuad(const Matrix<Real> & grad_u,
						  const Matrix<Real> & sigma,
						  Real & epot);

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:

  /// the Young's moduli for the different sub-materials
  Vector<Real> E;
    
  /// Poisson ratios for the different sub-materials
  Vector<Real> nu;
  
  /// First Lamé coefficient
  InternalField<Real> lambda;

  /// Second Lamé coefficient (shear modulus)
  InternalField<Real> mu;

  /// Bulk modulus
  InternalField<Real> kpa;

  /// default Lamé coefficients and bulk modulus
  Real lambda_default;
  Real mu_default;
  Real kpa_default;

  /// Plane stress or plane strain
  bool plane_stress;

};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

#include "material_igfem_elastic_inline_impl.cc"

__END_AKANTU__

#endif /* __AKANTU_MATERIAL_IGFEM_ELASTIC_HH__ */

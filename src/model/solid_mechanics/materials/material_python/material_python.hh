/**
 * @file   material_python.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 *
 */

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "material.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MATERIAL_PYTHON_HH__
#define __AKANTU_MATERIAL_PYTHON_HH__

/* -------------------------------------------------------------------------- */
#include "python_functor.hh"
/* -------------------------------------------------------------------------- */



__BEGIN_AKANTU__

class MaterialPython : public Material, PythonFunctor {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:

  MaterialPython(SolidMechanicsModel & model, PyObject * obj, const ID & id = "");

  virtual ~MaterialPython() {};

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:

  void registerInternals();
  
  virtual void initMaterial();


  /// constitutive law for all element of a type
  void computeStress(ElementType el_type, GhostType ghost_type = _not_ghost);

  /// constitutive law for a given quad point
  template <typename it_type>
  void  computeStress(Matrix<Real> & grad_u,
		      Matrix<Real> & sigma,
		      std::vector<it_type> & internal_iterators);

  /// compute the tangent stiffness matrix for an element type
  virtual void computeTangentModuli(const ElementType & el_type,
				    Array<Real> & tangent_matrix,
				    GhostType ghost_type = _not_ghost);

  /// compute the push wave speed of the material
  Real getPushWaveSpeed(const Element & element) const;

protected:
  /// update the dissipated energy, must be called after the stress have been computed
  //virtual void updateEnergies(ElementType el_type, GhostType ghost_type){};

  /// compute the tangent stiffness matrix for a given quadrature point
  //inline void computeTangentModuliOnQuad(Matrix<Real> & tangent, Real & dam){};

  /* ------------------------------------------------------------------------ */
  /* DataAccessor inherited members                                           */
  /* ------------------------------------------------------------------------ */
public:

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  //virtual Real getEnergy(std::string type){};
  //virtual Real getEnergy(std::string energy_id, ElementType type, UInt index){};

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:

  std::vector<Real> local_params;
  std::vector<InternalField<Real> *> internals;  
};

__END_AKANTU__

#endif /* __AKANTU_MATERIAL_PYTHON_HH__ */


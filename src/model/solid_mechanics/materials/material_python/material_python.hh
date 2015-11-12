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
#include <Python.h>
/* -------------------------------------------------------------------------- */



__BEGIN_AKANTU__

class MaterialPython : public Material {
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

  virtual void initMaterial();

  /// compute the tangent stiffness matrix for an element type
  virtual void computeTangentModuli(const ElementType & el_type,
				    Array<Real> & tangent_matrix,
				    GhostType ghost_type = _not_ghost){};

protected:
  /// update the dissipated energy, must be called after the stress have been computed
  virtual void updateEnergies(ElementType el_type, GhostType ghost_type){};

  /// compute the tangent stiffness matrix for a given quadrature point
  inline void computeTangentModuliOnQuad(Matrix<Real> & tangent, Real & dam){};

  /* ------------------------------------------------------------------------ */
  /* DataAccessor inherited members                                           */
  /* ------------------------------------------------------------------------ */
public:

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:

  virtual Real getEnergy(std::string type){};
  virtual Real getEnergy(std::string energy_id, ElementType type, UInt index){};

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:

};

__END_AKANTU__

#endif /* __AKANTU_MATERIAL_PYTHON_HH__ */


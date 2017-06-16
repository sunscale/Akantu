/**
 * @file   material_python.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 *
 * @date creation: Fri Nov 13 2015
 * @date last modification: Fri Nov 13 2015
 *
 * @brief  Python material
 *
 * @section LICENSE
 *
 * Copyright (©) 2015 EPFL (Ecole Polytechnique Fédérale de Lausanne) Laboratory
 * (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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

/**
 * @file   material_python.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 *
 */

/* -------------------------------------------------------------------------- */
#include "python_functor.hh"
/* -------------------------------------------------------------------------- */
#include "material.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_MATERIAL_PYTHON_HH__
#define __AKANTU_MATERIAL_PYTHON_HH__

/* -------------------------------------------------------------------------- */

namespace akantu {

class MaterialPython : public Material, PythonFunctor {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  MaterialPython(SolidMechanicsModel & model, PyObject * obj,
                 const ID & id = "");

  virtual ~MaterialPython(){};

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
  void computeStress(Matrix<Real> & grad_u, Matrix<Real> & sigma,
                     std::vector<it_type> & internal_iterators);

  /// compute the tangent stiffness matrix for an element type
  virtual void computeTangentModuli(const ElementType & el_type,
                                    Array<Real> & tangent_matrix,
                                    GhostType ghost_type = _not_ghost);

  /// compute the push wave speed of the material
  Real getPushWaveSpeed(const Element & element) const;

protected:
  /// update the dissipated energy, must be called after the stress have been
  /// computed
  // virtual void updateEnergies(ElementType el_type, GhostType ghost_type){};

  /// compute the tangent stiffness matrix for a given quadrature point
  // inline void computeTangentModuliOnQuad(Matrix<Real> & tangent, Real &
  // dam){};

  /* ------------------------------------------------------------------------ */
  /* DataAccessor inherited members                                           */
  /* ------------------------------------------------------------------------ */
public:
  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  // virtual Real getEnergy(std::string type){};
  // virtual Real getEnergy(std::string energy_id, ElementType type, UInt
  // index){};

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  std::map<std::string, Real> local_params;
  std::map<std::string, InternalField<Real> *> internals;
};

}  // akantu

#endif /* __AKANTU_MATERIAL_PYTHON_HH__ */

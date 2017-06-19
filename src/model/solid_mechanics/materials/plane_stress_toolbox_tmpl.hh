/**
 * @file   plane_stress_toolbox_tmpl.hh
 *
 * @author Daniel Pino Muñoz <daniel.pinomunoz@epfl.ch>
 *
 * @date creation: Tue Sep 16 2014
 * @date last modification: Tue Aug 18 2015
 *
 * @brief  2D specialization of the akantu::PlaneStressToolbox class
 *
 * @section LICENSE
 *
 * Copyright  (©)  2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de Lausanne)
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

#ifndef __AKANTU_PLANE_STRESS_TOOLBOX_TMPL_HH__
#define __AKANTU_PLANE_STRESS_TOOLBOX_TMPL_HH__

namespace akantu {

/* -------------------------------------------------------------------------- */
template <class ParentMaterial>
class PlaneStressToolbox<2, ParentMaterial> : public ParentMaterial {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  PlaneStressToolbox(SolidMechanicsModel & model, const ID & id = "");
  PlaneStressToolbox(SolidMechanicsModel & model, UInt dim, const Mesh & mesh,
                     FEEngine & fe_engine, const ID & id = "");

  virtual ~PlaneStressToolbox() {}

  AKANTU_GET_MACRO_BY_ELEMENT_TYPE_CONST(ThirdAxisDeformation,
                                         third_axis_deformation, Real);

protected:
  void initialize() {
    this->registerParam("Plane_Stress", plane_stress, false, _pat_parsmod,
                        "Is plane stress");
  }

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /* ------------------------------------------------------------------------ */
  virtual void initMaterial() {

    ParentMaterial::initMaterial();
    if (this->plane_stress && this->initialize_third_axis_deformation) {
      this->third_axis_deformation.initialize(1);
      this->third_axis_deformation.resize();
    }
  }

  /* ------------------------------------------------------------------------ */
  virtual void computeStress(ElementType el_type, GhostType ghost_type) {
    ParentMaterial::computeStress(el_type, ghost_type);
    if (this->plane_stress)
      computeThirdAxisDeformation(el_type, ghost_type);
  }

  /* ------------------------------------------------------------------------ */
  virtual void computeThirdAxisDeformation(__attribute__((unused))
                                           ElementType el_type,
                                           __attribute__((unused))
                                           GhostType ghost_type) {}

  /// Computation of Cauchy stress tensor in the case of finite deformation
  virtual void computeAllCauchyStresses(GhostType ghost_type = _not_ghost) {
    AKANTU_DEBUG_IN();

    if (this->plane_stress) {

      AKANTU_DEBUG_ASSERT(this->finite_deformation,
                          "The Cauchy stress can only be computed if you are "
                          "working in finite deformation.");

      // resizeInternalArray(stress);

      Mesh::type_iterator it =
          this->model->getFEEngine().getMesh().firstType(2, ghost_type);
      Mesh::type_iterator last_type =
          this->model->getFEEngine().getMesh().lastType(2, ghost_type);

      for (; it != last_type; ++it)
        this->computeCauchyStressPlaneStress(*it, ghost_type);
    } else
      ParentMaterial::computeAllCauchyStresses(ghost_type);

    AKANTU_DEBUG_OUT();
  }

  virtual void
  computeCauchyStressPlaneStress(__attribute__((unused)) ElementType el_type,
                                 __attribute__((unused))
                                 GhostType ghost_type = _not_ghost){};

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  /// third axis strain measure value
  InternalField<Real> third_axis_deformation;

  /// Plane stress or plane strain
  bool plane_stress;

  /// For non linear materials, the \epsilon_{zz} might be required
  bool initialize_third_axis_deformation;
};

template <class ParentMaterial>
inline PlaneStressToolbox<2, ParentMaterial>::PlaneStressToolbox(
    SolidMechanicsModel & model, const ID & id)
    : Material(model, id), ParentMaterial(model, id),
      third_axis_deformation("third_axis_deformation", *this),
      plane_stress(false), initialize_third_axis_deformation(false) {

  /// @todo Plane_Stress should not be possible to be modified after
  /// initMaterial (but before)
  this->initialize();
}

template <class ParentMaterial>
inline PlaneStressToolbox<2, ParentMaterial>::PlaneStressToolbox(
    SolidMechanicsModel & model, UInt dim, const Mesh & mesh,
    FEEngine & fe_engine, const ID & id)
    : Material(model, dim, mesh, fe_engine, id),
      ParentMaterial(model, dim, mesh, fe_engine, id),
      third_axis_deformation("third_axis_deformation", *this, dim, fe_engine,
                             this->element_filter),
      plane_stress(false), initialize_third_axis_deformation(false) {
  this->initialize();
}

template <>
inline PlaneStressToolbox<2, Material>::PlaneStressToolbox(
    SolidMechanicsModel & model, const ID & id)
    : Material(model, id),
      third_axis_deformation("third_axis_deformation", *this),
      plane_stress(false), initialize_third_axis_deformation(false) {

  /// @todo Plane_Stress should not be possible to be modified after
  /// initMaterial (but before)
  this->registerParam("Plane_Stress", plane_stress, false, _pat_parsmod,
                      "Is plane stress");
}

} // akantu

#endif /* __AKANTU_PLANE_STRESS_TOOLBOX_TMPL_HH__ */

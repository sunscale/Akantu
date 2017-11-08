/**
 * @file   plane_stress_toolbox.hh
 *
 * @author Daniel Pino Muñoz <daniel.pinomunoz@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Tue Sep 16 2014
 * @date last modification: Tue Aug 18 2015
 *
 * @brief  Tools to implement the plane stress behavior in a material
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

#ifndef __AKANTU_PLANE_STRESS_TOOLBOX_HH__
#define __AKANTU_PLANE_STRESS_TOOLBOX_HH__

namespace akantu {
class SolidMechanicsModel;
class FEEngine;
}  // akantu


namespace akantu {

/**
 * Empty class in dimensions different from 2
 * This class is only specialized for 2D in the tmpl file
 */
template <UInt dim, class ParentMaterial = Material>
class PlaneStressToolbox : public ParentMaterial {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  PlaneStressToolbox(SolidMechanicsModel & model, const ID & id = "")
      : ParentMaterial(model, id) {}
  PlaneStressToolbox(SolidMechanicsModel & model, UInt spatial_dimension,
                     const Mesh & mesh, FEEngine & fe_engine,
                     const ID & id = "")
      : ParentMaterial(model, spatial_dimension, mesh, fe_engine, id) {}

  ~PlaneStressToolbox() override = default;

protected:
  void initialize();

public:
  void computeAllCauchyStresses(GhostType ghost_type = _not_ghost) override {
    AKANTU_DEBUG_IN();

    ParentMaterial::computeAllCauchyStresses(ghost_type);

    AKANTU_DEBUG_OUT();
  }

  virtual void computeCauchyStressPlaneStress(__attribute__((unused))
                                              ElementType el_type,
                                              __attribute__((unused))
                                              GhostType ghost_type) {
    AKANTU_DEBUG_IN();

    AKANTU_DEBUG_ERROR("The function \"computeCauchyStressPlaneStress\" can "
                       "only be used in 2D Plane stress problems, which means "
                       "that you made a mistake somewhere!! ");

    AKANTU_DEBUG_OUT();
  }

protected:
  bool initialize_third_axis_deformation;
};

#define AKANTU_PLANE_STRESS_TOOL_SPEC(dim)                                     \
  template <>                                                                  \
  inline PlaneStressToolbox<dim, Material>::PlaneStressToolbox(                \
      SolidMechanicsModel & model, const ID & id)                              \
      : Material(model, id) {}

AKANTU_PLANE_STRESS_TOOL_SPEC(1)
AKANTU_PLANE_STRESS_TOOL_SPEC(3)

} // namespace akantu

#include "plane_stress_toolbox_tmpl.hh"

#endif /* __AKANTU_PLANE_STRESS_TOOLBOX_HH__ */

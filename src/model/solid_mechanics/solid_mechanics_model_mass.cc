/**
 * @file   solid_mechanics_model_mass.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Tue Oct 05 2010
 * @date last modification: Wed Nov 08 2017
 *
 * @brief  function handling mass computation
 *
 *
 * Copyright (©)  2010-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#include "integrator_gauss.hh"
#include "material.hh"
#include "model_solver.hh"
#include "shape_lagrange.hh"
#include "solid_mechanics_model.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

class ComputeRhoFunctor {
public:
  explicit ComputeRhoFunctor(const SolidMechanicsModel & model)
      : model(model){};

  void operator()(Matrix<Real> & rho, const Element & element) {
    const Array<UInt> & mat_indexes =
        model.getMaterialByElement(element.type, element.ghost_type);
    Real mat_rho =
        model.getMaterial(mat_indexes(element.element)).getParam("rho");
    rho.set(mat_rho);
  }

private:
  const SolidMechanicsModel & model;
};

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::assembleMassLumped() {
  AKANTU_DEBUG_IN();

  if (not need_to_reassemble_lumped_mass) {
    return;
  }

  this->allocNodalField(this->mass, spatial_dimension, "mass");
  mass->zero();

  if (!this->getDOFManager().hasLumpedMatrix("M")) {
    this->getDOFManager().getNewLumpedMatrix("M");
  }

  this->getDOFManager().zeroLumpedMatrix("M");

  assembleMassLumped(_not_ghost);
  assembleMassLumped(_ghost);

  this->getDOFManager().getLumpedMatrixPerDOFs("displacement", "M",
                                               *(this->mass));

/// for not connected nodes put mass to one in order to avoid
#if !defined(AKANTU_NDEBUG)
  bool has_unconnected_nodes = false;
  auto mass_it = mass->begin_reinterpret(mass->size() * mass->getNbComponent());
  auto mass_end = mass->end_reinterpret(mass->size() * mass->getNbComponent());
  for (; mass_it != mass_end; ++mass_it) {
    if (std::abs(*mass_it) < std::numeric_limits<Real>::epsilon() ||
        Math::isnan(*mass_it)) {
      has_unconnected_nodes = true;
      break;
    }
  }

  if (has_unconnected_nodes) {
    AKANTU_DEBUG_WARNING("There are nodes that seem to not be connected to any "
                         "elements, beware that they have lumped mass of 0.");
  }
#endif

  this->synchronize(SynchronizationTag::_smm_mass);

  need_to_reassemble_lumped_mass = false;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::assembleMass() {
  AKANTU_DEBUG_IN();

  if (not need_to_reassemble_mass) {
    return;
  }

  this->getDOFManager().zeroMatrix("M");
  assembleMass(_not_ghost);

  need_to_reassemble_mass = false;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::assembleMassLumped(GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  auto & fem = getFEEngineClass<MyFEEngineType>();
  ComputeRhoFunctor compute_rho(*this);

  for (auto type :
       mesh.elementTypes(Model::spatial_dimension, ghost_type, _ek_regular)) {
    fem.assembleFieldLumped(compute_rho, "M", "displacement",
                            this->getDOFManager(), type, ghost_type);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::assembleMass(GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  auto & fem = getFEEngineClass<MyFEEngineType>();
  ComputeRhoFunctor compute_rho(*this);

  for (auto type :
       mesh.elementTypes(Model::spatial_dimension, ghost_type, _ek_regular)) {
    fem.assembleFieldMatrix(compute_rho, "M", "displacement",
                            this->getDOFManager(), type, ghost_type);
  }

  AKANTU_DEBUG_OUT();
}

} // namespace akantu

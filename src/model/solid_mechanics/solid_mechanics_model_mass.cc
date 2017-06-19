/**
 * @file   solid_mechanics_model_mass.cc
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Tue Oct 05 2010
 * @date last modification: Fri Oct 16 2015
 *
 * @brief  function handling mass computation
 *
 * @section LICENSE
 *
 * Copyright (©)  2010-2012, 2014,  2015 EPFL  (Ecole Polytechnique  Fédérale de
 * Lausanne)  Laboratory (LSMS  -  Laboratoire de  Simulation  en Mécanique  des
 * Solides)
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
#include "material.hh"
#include "model_solver.hh"
#include "solid_mechanics_model.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::assembleMassLumped() {
  AKANTU_DEBUG_IN();

  UInt nb_nodes = mesh.getNbNodes();

  if (this->mass == NULL) {
    std::stringstream sstr_mass;
    sstr_mass << id << ":mass";
    mass = &(alloc<Real>(sstr_mass.str(), nb_nodes, spatial_dimension, 0));
  } else {
    mass->clear();
  }

  if (!this->getDOFManager().hasLumpedMatrix("M")) {
    this->getDOFManager().getNewLumpedMatrix("M");
  }

  this->getDOFManager().clearLumpedMatrix("M");

  assembleMassLumped(_not_ghost);
  assembleMassLumped(_ghost);

  this->getDOFManager().getLumpedMatrixPerDOFs("displacement", "M",
                                               *(this->mass));

/// for not connected nodes put mass to one in order to avoid
#if !defined(AKANTU_NDEBUG)
  bool has_unconnected_nodes = false;
  auto mass_it =
      mass->begin_reinterpret(mass->getSize() * mass->getNbComponent());
  auto mass_end =
      mass->end_reinterpret(mass->getSize() * mass->getNbComponent());
  for (; mass_it != mass_end; ++mass_it) {
    if (std::abs(*mass_it) < std::numeric_limits<Real>::epsilon() ||
        Math::isnan(*mass_it)) {
      has_unconnected_nodes = true;
      break;
    }
  }

  if(has_unconnected_nodes)
    AKANTU_DEBUG_WARNING("There are nodes that seem to not be connected to any "
                         "elements, beware that they have lumped mass of 0.");
#endif

  this->synchronize(_gst_smm_mass);
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::assembleMassLumped(GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  FEEngine & fem = getFEEngine();

  Array<Real> rho(0, spatial_dimension);

  Mesh::type_iterator it = mesh.firstType(spatial_dimension, ghost_type);
  Mesh::type_iterator end = mesh.lastType(spatial_dimension, ghost_type);
  for (; it != end; ++it) {
    ElementType type = *it;

    computeRho(rho, type, ghost_type);

    fem.assembleFieldLumped(rho, "M", "displacement", this->getDOFManager(),
                            type, ghost_type);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::assembleMass() {
  AKANTU_DEBUG_IN();

  if (!this->getDOFManager().hasMatrix("M")) {
    this->getDOFManager().getNewMatrix("M", "J");
  }

  this->getDOFManager().clearMatrix("M");
  assembleMass(_not_ghost);

  AKANTU_DEBUG_OUT();
}

class ComputeRhoFunctor {
public:
  explicit ComputeRhoFunctor(const SolidMechanicsModel & model) : model(model){};

  void operator()(Matrix<Real> & rho, const Element & element,
                  __attribute__((unused)) const Matrix<Real> quad_coords) const {
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
void SolidMechanicsModel::assembleMass(GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  MyFEEngineType & fem = getFEEngineClass<MyFEEngineType>();

  ComputeRhoFunctor compute_rho(*this);

  Mesh::type_iterator it = mesh.firstType(spatial_dimension, ghost_type);
  Mesh::type_iterator end = mesh.lastType(spatial_dimension, ghost_type);
  for (; it != end; ++it) {
    ElementType type = *it;
    fem.assembleFieldMatrix(compute_rho, "M", "displacement",
                            this->getDOFManager(), type, ghost_type);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void SolidMechanicsModel::computeRho(Array<Real> & rho, ElementType type,
                                     GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  FEEngine & fem = this->getFEEngine();
  UInt nb_element = fem.getMesh().getNbElement(type, ghost_type);

  Array<UInt> & mat_indexes = this->material_index(type, ghost_type);

  UInt nb_quadrature_points = fem.getNbIntegrationPoints(type);

  rho.resize(nb_element * nb_quadrature_points);
  Array<Real>::vector_iterator rho_it = rho.begin(spatial_dimension);

  /// compute @f$ rho @f$ for each nodes of each element
  for (UInt el = 0; el < nb_element; ++el) {
    /// here rho is constant in an element
    Real mat_rho = this->materials[mat_indexes(el)]->getParam("rho");

    for (UInt n = 0; n < nb_quadrature_points; ++n, ++rho_it) {
      (*rho_it).set(mat_rho);
    }
  }

  AKANTU_DEBUG_OUT();
}

} // akantu

/**
 * @file   structural_mechanics_model_mass.cc
 *
 * @author Lucas Frerot <lucas.frerot@epfl.ch>
 * @author Sébastien Hartmann <sebastien.hartmann@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Mon Jul 07 2014
 * @date last modification: Fri Dec 15 2017
 *
 * @brief  function handling mass computation
 *
 *
 * Copyright (©) 2014-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "shape_structural.hh"
#include "structural_mechanics_model.hh"

/* -------------------------------------------------------------------------- */
namespace akantu {

class ComputeRhoFunctorStruct {
public:
  explicit ComputeRhoFunctorStruct(const StructuralMechanicsModel & model)
      : model(model){};

  void operator()(Matrix<Real> & rho, const Element & element) const {
    Real mat_rho = model.getMaterial(element).rho;
    rho.set(mat_rho);
  }

private:
  const StructuralMechanicsModel & model;
};

/* -------------------------------------------------------------------------- */
void StructuralMechanicsModel::assembleMass() {
  AKANTU_DEBUG_IN();

  assembleMass(_not_ghost);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void StructuralMechanicsModel::assembleMass(GhostType ghost_type) {
  AKANTU_DEBUG_IN();
  MyFEEngineType & fem = getFEEngineClass<MyFEEngineType>();
  ComputeRhoFunctorStruct compute_rho(*this);

  for (auto type :
       mesh.elementTypes(spatial_dimension, ghost_type, _ek_structural)) {
    fem.assembleFieldMatrix(compute_rho, "M", "displacement",
                            this->getDOFManager(), type, ghost_type);
  }
  AKANTU_DEBUG_OUT();
}

} // namespace akantu

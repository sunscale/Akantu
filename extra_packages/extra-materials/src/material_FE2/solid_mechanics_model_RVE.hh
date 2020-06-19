/**
 * @file   solid_mechanics_model_RVE.hh
 * @author Aurelia Isabel Cuba Ramos <aurelia.cubaramos@epfl.ch>
 * @date   Wed Jan 13 14:54:18 2016
 *
 * @brief  SMM for RVE computations in FE2 simulations
 *
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#ifndef __AKANTU_SOLID_MECHANICS_MODEL_RVE_HH__
#define __AKANTU_SOLID_MECHANICS_MODEL_RVE_HH__

/* -------------------------------------------------------------------------- */
#include "aka_grid_dynamic.hh"
#include "solid_mechanics_model.hh"
#include <unordered_set>
/* -------------------------------------------------------------------------- */

namespace akantu {

class SolidMechanicsModelRVE : public SolidMechanicsModel {

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */

public:
  SolidMechanicsModelRVE(Mesh & mesh, bool use_RVE_mat_selector = true,
                         UInt nb_gel_pockets = 400,
                         UInt spatial_dimension = _all_dimensions,
                         const ID & id = "solid_mechanics_model",
                         const MemoryID & memory_id = 0);

  virtual ~SolidMechanicsModelRVE();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
protected:
  void initFullImpl(const ModelOptions & option) override;

  /// initialize the materials
  void initMaterials() override;

public:
  /// apply boundary contions based on macroscopic deformation gradient
  virtual void
  applyBoundaryConditions(const Matrix<Real> & displacement_gradient);

  /// apply homogeneous temperature field from the macroscale level to the RVEs
  virtual void applyHomogeneousTemperature(const Real & temperature);

  /// advance the reactions -> grow gel and apply homogenized properties
  void advanceASR(const Matrix<Real> & prestrain);

  /// compute average stress or strain in the model
  Real averageTensorField(UInt row_index, UInt col_index,
                          const ID & field_type);

  /// compute effective stiffness of the RVE
  void homogenizeStiffness(Matrix<Real> & C_macro);

  /// compute average eigenstrain
  void homogenizeEigenGradU(Matrix<Real> & eigen_gradu_macro);

  /* ------------------------------------------------------------------------ */
  /* Data Accessor inherited members                                          */
  /* ------------------------------------------------------------------------ */

  inline void unpackData(CommunicationBuffer & buffer,
                         const Array<UInt> & index,
                         const SynchronizationTag & tag) override;

  /* ------------------------------------------------------------------------ */
  /* Accessors */
  /* ------------------------------------------------------------------------ */
public:
  AKANTU_GET_MACRO(CornerNodes, corner_nodes, const Array<UInt> &);
  AKANTU_GET_MACRO(Volume, volume, Real);

private:
  /// find the corner nodes
  void findCornerNodes();

  /// perform virtual testing
  void performVirtualTesting(const Matrix<Real> & H,
                             Matrix<Real> & eff_stresses,
                             Matrix<Real> & eff_strains, const UInt test_no);

  void fillCracks(ElementTypeMapReal & saved_damage);
  void drainCracks(const ElementTypeMapReal & saved_damage);
  /* ------------------------------------------------------------------------ */
  /* Members */
  /* ------------------------------------------------------------------------ */

  /// volume of the RVE
  Real volume;

  /// corner nodes 1, 2, 3, 4 (see Leonardo's thesis, page 98)
  Array<UInt> corner_nodes;

  /// bottom nodes
  std::unordered_set<UInt> bottom_nodes;

  /// left nodes
  std::unordered_set<UInt> left_nodes;

  /// standard mat selector or user one
  bool use_RVE_mat_selector;

  /// the number of gel pockets inside the RVE
  UInt nb_gel_pockets;

  /// dump counter
  UInt nb_dumps;
};

inline void SolidMechanicsModelRVE::unpackData(CommunicationBuffer & buffer,
                                               const Array<UInt> & index,
                                               const SynchronizationTag & tag) {
  SolidMechanicsModel::unpackData(buffer, index, tag);

  //  if (tag == SynchronizationTag::_smm_uv) {
  //    auto disp_it = displacement->begin(spatial_dimension);
  //
  //    for (auto node : index) {
  //      Vector<Real> current_disp(disp_it[node]);
  //
  //      // if node is at the bottom, u_bottom = u_top +u_2 -u_3
  //      if (bottom_nodes.count(node)) {
  //        current_disp += Vector<Real>(disp_it[corner_nodes(1)]);
  //        current_disp -= Vector<Real>(disp_it[corner_nodes(2)]);
  //      }
  //      // if node is at the left, u_left = u_right +u_4 -u_3
  //      else if (left_nodes.count(node)) {
  //        current_disp += Vector<Real>(disp_it[corner_nodes(3)]);
  //        current_disp -= Vector<Real>(disp_it[corner_nodes(2)]);
  //      }
  //    }
  //  }
}

/* -------------------------------------------------------------------------- */
/* ASR material selector                                                      */
/* -------------------------------------------------------------------------- */
class GelMaterialSelector : public MeshDataMaterialSelector<std::string> {
public:
  GelMaterialSelector(SolidMechanicsModel & model, const Real box_size,
                      const std::string & gel_material,
                      const UInt nb_gel_pockets, Real /*tolerance*/ = 0.)
      : MeshDataMaterialSelector<std::string>("physical_names", model),
        model(model), gel_material(gel_material),
        nb_gel_pockets(nb_gel_pockets), nb_placed_gel_pockets(0),
        box_size(box_size) {
    Mesh & mesh = this->model.getMesh();
    UInt spatial_dimension = model.getSpatialDimension();
    Element el{_triangle_3, 0, _not_ghost};
    UInt nb_element = mesh.getNbElement(el.type, el.ghost_type);
    Array<Real> barycenter(nb_element, spatial_dimension);

    for (auto && data : enumerate(make_view(barycenter, spatial_dimension))) {
      el.element = std::get<0>(data);
      auto & bary = std::get<1>(data);
      mesh.getBarycenter(el, bary);
    }

    /// generate the gel pockets
    srand(0.);
    Vector<Real> center(spatial_dimension);
    UInt placed_gel_pockets = 0;
    std::set<int> checked_baries;
    while (placed_gel_pockets != nb_gel_pockets) {
      /// get a random bary center
      UInt bary_id = rand() % nb_element;
      if (checked_baries.find(bary_id) != checked_baries.end())
        continue;
      checked_baries.insert(bary_id);
      el.element = bary_id;
      if (MeshDataMaterialSelector<std::string>::operator()(el) == 1)
        continue; /// element belongs to paste
      gel_pockets.push_back(el);
      placed_gel_pockets += 1;
    }
  }

  UInt operator()(const Element & elem) {
    UInt temp_index = MeshDataMaterialSelector<std::string>::operator()(elem);
    if (temp_index == 1)
      return temp_index;
    std::vector<Element>::const_iterator iit = gel_pockets.begin();
    std::vector<Element>::const_iterator eit = gel_pockets.end();
    if (std::find(iit, eit, elem) != eit) {
      nb_placed_gel_pockets += 1;
      std::cout << nb_placed_gel_pockets << " gelpockets placed" << std::endl;
      return model.getMaterialIndex(gel_material);
      ;
    }
    return 0;
  }

protected:
  SolidMechanicsModel & model;
  std::string gel_material;
  std::vector<Element> gel_pockets;
  UInt nb_gel_pockets;
  UInt nb_placed_gel_pockets;
  Real box_size;
};

} // namespace akantu

///#include "material_selector_tmpl.hh"

#endif /* __AKANTU_SOLID_MECHANICS_MODEL_RVE_HH__ */

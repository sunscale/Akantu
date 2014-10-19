/**
 * @file   dumper_material_padders.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Tue Sep 02 2014
 * @date last modification: Fri Sep 19 2014
 *
 * @brief  Material padders for plane stress/ plane strain
 *
 * @section LICENSE
 *
 * Copyright (©) 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

#ifndef __AKANTU_DUMPER_MATERIAL_PADDERS_HH__
#define __AKANTU_DUMPER_MATERIAL_PADDERS_HH__
/* -------------------------------------------------------------------------- */
#include "dumper_padding_helper.hh"
/* -------------------------------------------------------------------------- */
__BEGIN_AKANTU__
__BEGIN_AKANTU_DUMPER__
/* -------------------------------------------------------------------------- */


template<class T, class R>
class MaterialPadder : public PadderGeneric<Vector<T>, R > {

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */

public:
  MaterialPadder(const SolidMechanicsModel & model) :
    model(model),
    element_index_by_material(model.getElementIndexByMaterial()) { }

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */

  /// return the material from the global element index
  const Material & getMaterialFromGlobalIndex(Element global_index){

    UInt index = global_index.getIndex();
    UInt material_id = element_index_by_material(global_index.getType())(index);
    const Material & material = model.getMaterial(material_id);
    return material;
  }

  /// return the type of the element from global index
  ElementType getElementTypeFromGlobalIndex(Element global_index){
    return global_index.getType();
  }

protected:

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */

  /// all material padders probably need access to solid mechanics model
  const SolidMechanicsModel & model;
  /// they also need an access to the map from global ids to material id and local ids
  const ElementTypeMapArray<UInt>  & element_index_by_material;
  /// the number of data per element
  const ElementTypeMapArray<UInt>  nb_data_per_element;

};

/* -------------------------------------------------------------------------- */

template <UInt spatial_dimension>
class StressPadder :
  public MaterialPadder<Real,Matrix<Real> > {

public:
  StressPadder(const SolidMechanicsModel & model) :
    MaterialPadder<Real, Matrix<Real> >(model){
    this->setPadding(3,3);
  }

  inline Matrix<Real> func(const Vector<Real> & in, Element global_element_id){

    UInt nrows = spatial_dimension;
    UInt ncols = in.size() / nrows;
    UInt nb_data = in.size() / (ncols*ncols);

    Matrix<Real> stress = this->pad(in, nrows,ncols, nb_data);
    const Material & material = this->getMaterialFromGlobalIndex(global_element_id);
    bool plane_strain = true;
    if(spatial_dimension == 2)
      plane_strain = !material.getParam<bool>("Plane_Stress");

    if(plane_strain) {
      Real nu = material.getParam<Real>("nu");
      for (UInt d = 0; d < nb_data; ++d) {
	stress(2, 2 + 3*d) = nu * (stress(0, 0 + 3*d) + stress(1, 1 + 3*d));
      }
    }
    return stress;
  }

  UInt getDim(){return 9;};

  UInt getNbComponent(UInt old_nb_comp){
    return this->getDim();
  };

};

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
class StrainPadder : public MaterialPadder<Real, Matrix<Real> > {
public:
  StrainPadder(const SolidMechanicsModel & model) :
    MaterialPadder<Real, Matrix<Real> >(model) {
    this->setPadding(3,3);
  }

  inline Matrix<Real> func(const Vector<Real> & in, Element global_element_id){

    UInt nrows = spatial_dimension;
    UInt ncols = in.size() / nrows;
    UInt nb_data = in.size() / ncols;

    Matrix<Real> strain = this->pad(in, nrows,ncols, nb_data);
    const Material & material = this->getMaterialFromGlobalIndex(global_element_id);
    bool plane_stress = material.getParam<bool>("Plane_Stress");
    if(plane_stress) {
      Real nu = material.getParam<Real>("nu");
      for (UInt d = 0; d < nb_data; ++d) {
	strain(2, 2 + 3*d) = nu / (nu - 1) * (strain(0, 0 + 3*d) + strain(1, 1 + 3*d));
      }
    }
    return strain;
  }

  UInt getDim(){return 9;};

  UInt getNbComponent(UInt old_nb_comp){
    return this->getDim();
  };

};

/* -------------------------------------------------------------------------- */

__END_AKANTU_DUMPER__
__END_AKANTU__


#endif /* __AKANTU_DUMPER_MATERIAL_PADDERS_HH__ */

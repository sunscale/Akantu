/**
 * @file   dumper_material_padders.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 *
 * @date creation: Tue Sep 02 2014
 * @date last modification: Wed Nov 29 2017
 *
 * @brief  Material padders for plane stress/ plane strain
 *
 * @section LICENSE
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

#ifndef __AKANTU_DUMPER_MATERIAL_PADDERS_HH__
#define __AKANTU_DUMPER_MATERIAL_PADDERS_HH__
/* -------------------------------------------------------------------------- */
#include "dumper_padding_helper.hh"
/* -------------------------------------------------------------------------- */
namespace akantu {
__BEGIN_AKANTU_DUMPER__
/* -------------------------------------------------------------------------- */
class MaterialFunctor {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  MaterialFunctor(const SolidMechanicsModel & model)
      : model(model), material_index(model.getMaterialByElement()),
        nb_data_per_element("nb_data_per_element", model.getID(),
                            model.getMemoryID()),
        spatial_dimension(model.getSpatialDimension()) {}

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
  /// return the material from the global element index
  const Material & getMaterialFromGlobalIndex(Element global_index) {
    UInt index = global_index.element;
    UInt material_id = material_index(global_index.type)(index);
    const Material & material = model.getMaterial(material_id);
    return material;
  }

  /// return the type of the element from global index
  ElementType getElementTypeFromGlobalIndex(Element global_index) {
    return global_index.type;
  }

protected:
  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */

  /// all material padders probably need access to solid mechanics model
  const SolidMechanicsModel & model;

  /// they also need an access to the map from global ids to material id and
  /// local ids
  const ElementTypeMapArray<UInt> & material_index;

  /// the number of data per element
  const ElementTypeMapArray<UInt> nb_data_per_element;

  UInt spatial_dimension;
};

/* -------------------------------------------------------------------------- */
template <class T, class R>
class MaterialPadder : public MaterialFunctor,
                       public PadderGeneric<Vector<T>, R> {
public:
  MaterialPadder(const SolidMechanicsModel & model) : MaterialFunctor(model) {}
};

/* -------------------------------------------------------------------------- */

template <UInt spatial_dimension>
class StressPadder : public MaterialPadder<Real, Matrix<Real>> {

public:
  StressPadder(const SolidMechanicsModel & model)
      : MaterialPadder<Real, Matrix<Real>>(model) {
    this->setPadding(3, 3);
  }

  inline Matrix<Real> func(const Vector<Real> & in,
                           Element global_element_id) override {
    UInt nrows = spatial_dimension;
    UInt ncols = in.size() / nrows;
    UInt nb_data = in.size() / (nrows * nrows);

    Matrix<Real> stress = this->pad(in, nrows, ncols, nb_data);
    const Material & material =
        this->getMaterialFromGlobalIndex(global_element_id);
    bool plane_strain = true;
    if (spatial_dimension == 2)
      plane_strain = !((bool)material.getParam("Plane_Stress"));

    if (plane_strain) {
      Real nu = material.getParam("nu");
      for (UInt d = 0; d < nb_data; ++d) {
        stress(2, 2 + 3 * d) =
            nu * (stress(0, 0 + 3 * d) + stress(1, 1 + 3 * d));
      }
    }
    return stress;
  }

  UInt getDim() override { return 9; };

  UInt getNbComponent(UInt /*old_nb_comp*/) override { return this->getDim(); };
};

/* -------------------------------------------------------------------------- */
template <UInt spatial_dimension>
class StrainPadder : public MaterialFunctor,
                     public PadderGeneric<Matrix<Real>, Matrix<Real>> {
public:
  StrainPadder(const SolidMechanicsModel & model) : MaterialFunctor(model) {
    this->setPadding(3, 3);
  }

  inline Matrix<Real> func(const Matrix<Real> & in,
                           Element global_element_id) override {
    UInt nrows = spatial_dimension;
    UInt nb_data = in.size() / (nrows * nrows);

    Matrix<Real> strain = this->pad(in, nb_data);
    const Material & material =
        this->getMaterialFromGlobalIndex(global_element_id);
    bool plane_stress = material.getParam("Plane_Stress");

    if (plane_stress) {
      Real nu = material.getParam("nu");
      for (UInt d = 0; d < nb_data; ++d) {
        strain(2, 2 + 3 * d) =
            nu / (nu - 1) * (strain(0, 0 + 3 * d) + strain(1, 1 + 3 * d));
      }
    }

    return strain;
  }

  UInt getDim() override { return 9; };

  UInt getNbComponent(UInt /*old_nb_comp*/) override { return this->getDim(); };
};

/* -------------------------------------------------------------------------- */
template <bool green_strain>
class ComputeStrain : public MaterialFunctor,
                      public ComputeFunctor<Vector<Real>, Matrix<Real>> {
public:
  ComputeStrain(const SolidMechanicsModel & model) : MaterialFunctor(model) {}

  inline Matrix<Real> func(const Vector<Real> & in,
                           Element /*global_element_id*/) override {
    UInt nrows = spatial_dimension;
    UInt ncols = in.size() / nrows;
    UInt nb_data = in.size() / (nrows * nrows);

    Matrix<Real> ret_all_strain(nrows, ncols);
    Tensor3<Real> all_grad_u(in.storage(), nrows, nrows, nb_data);
    Tensor3<Real> all_strain(ret_all_strain.storage(), nrows, nrows, nb_data);

    for (UInt d = 0; d < nb_data; ++d) {
      Matrix<Real> grad_u = all_grad_u(d);
      Matrix<Real> strain = all_strain(d);

      if (spatial_dimension == 2) {
        if (green_strain)
          Material::gradUToE<2>(grad_u, strain);
        else
          Material::gradUToEpsilon<2>(grad_u, strain);
      } else if (spatial_dimension == 3) {
        if (green_strain)
          Material::gradUToE<3>(grad_u, strain);
        else
          Material::gradUToEpsilon<3>(grad_u, strain);
      }
    }

    return ret_all_strain;
  }

  UInt getDim() override { return spatial_dimension * spatial_dimension; };

  UInt getNbComponent(UInt /*old_nb_comp*/) override { return this->getDim(); };
};

/* -------------------------------------------------------------------------- */
template <bool green_strain>
class ComputePrincipalStrain
    : public MaterialFunctor,
      public ComputeFunctor<Vector<Real>, Matrix<Real>> {
public:
  ComputePrincipalStrain(const SolidMechanicsModel & model)
      : MaterialFunctor(model) {}

  inline Matrix<Real> func(const Vector<Real> & in,
                           Element /*global_element_id*/) override {
    UInt nrows = spatial_dimension;
    UInt nb_data = in.size() / (nrows * nrows);

    Matrix<Real> ret_all_strain(nrows, nb_data);
    Tensor3<Real> all_grad_u(in.storage(), nrows, nrows, nb_data);
    Matrix<Real> strain(nrows, nrows);

    for (UInt d = 0; d < nb_data; ++d) {
      Matrix<Real> grad_u = all_grad_u(d);

      if (spatial_dimension == 2) {
        if (green_strain)
          Material::gradUToE<2>(grad_u, strain);
        else
          Material::gradUToEpsilon<2>(grad_u, strain);
      } else if (spatial_dimension == 3) {
        if (green_strain)
          Material::gradUToE<3>(grad_u, strain);
        else
          Material::gradUToEpsilon<3>(grad_u, strain);
      }

      Vector<Real> principal_strain(ret_all_strain(d));
      strain.eig(principal_strain);
    }

    return ret_all_strain;
  }

  UInt getDim() override { return spatial_dimension; };

  UInt getNbComponent(UInt /*old_nb_comp*/) override { return this->getDim(); };
};

/* -------------------------------------------------------------------------- */
class ComputeVonMisesStress
    : public MaterialFunctor,
      public ComputeFunctor<Vector<Real>, Vector<Real>> {
public:
  ComputeVonMisesStress(const SolidMechanicsModel & model)
      : MaterialFunctor(model) {}

  inline Vector<Real> func(const Vector<Real> & in,
                           Element /*global_element_id*/) override {
    UInt nrows = spatial_dimension;
    UInt nb_data = in.size() / (nrows * nrows);

    Vector<Real> von_mises_stress(nb_data);
    Matrix<Real> deviatoric_stress(3, 3);

    for (UInt d = 0; d < nb_data; ++d) {
      Matrix<Real> cauchy_stress(in.storage() + d * nrows * nrows, nrows,
                                 nrows);
      von_mises_stress(d) = Material::stressToVonMises(cauchy_stress);
    }

    return von_mises_stress;
  }

  UInt getDim() override { return 1; };

  UInt getNbComponent(UInt /*old_nb_comp*/) override { return this->getDim(); };
};

/* -------------------------------------------------------------------------- */

__END_AKANTU_DUMPER__
} // namespace akantu

#endif /* __AKANTU_DUMPER_MATERIAL_PADDERS_HH__ */

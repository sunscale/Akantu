/**
 * @file   material_inline_impl.cc
 *
 * @author Marco Vocialta <marco.vocialta@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Daniel Pino Muñoz <daniel.pinomunoz@epfl.ch>
 *
 * @date creation: Tue Jul 27 2010
 * @date last modification: Tue Sep 16 2014
 *
 * @brief  Implementation of the inline functions of the class material
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

__END_AKANTU__

#include "solid_mechanics_model.hh"
#include <iostream>

__BEGIN_AKANTU__


/* -------------------------------------------------------------------------- */
inline UInt Material::addElement(const ElementType & type,
				 UInt element,
				 const GhostType & ghost_type) {
  Array<UInt> & el_filter = element_filter(type, ghost_type);
  el_filter.push_back(element);
  return el_filter.getSize()-1;
}

/* -------------------------------------------------------------------------- */
inline UInt Material::getTangentStiffnessVoigtSize(UInt dim) const {
  return (dim * (dim - 1) / 2 + dim);
}
/* -------------------------------------------------------------------------- */
inline UInt Material::getCauchyStressMatrixSize(UInt dim) const {
  return (dim * dim);
}
/* -------------------------------------------------------------------------- */
template<UInt dim>
inline void Material::gradUToF(const Matrix<Real> & grad_u,
			       Matrix<Real> & F) const {
  AKANTU_DEBUG_ASSERT(F.size() >= grad_u.size() && grad_u.size() == dim*dim,
            "The dimension of the tensor F should be greater or equal to the dimension of the tensor grad_u.");

  F.eye();

  for (UInt i = 0; i < dim; ++i)
    for (UInt j = 0; j < dim; ++j)
      F(i, j) += grad_u(i, j);
}

/* -------------------------------------------------------------------------- */
template<UInt dim >
inline void Material::computeCauchyStressOnQuad(const Matrix<Real> & F,
						const Matrix<Real> & piola,
						Matrix<Real> & sigma,
                                                const Real & C33 ) const {

  Real J = F.det() * sqrt(C33);

  Matrix<Real> F_S(dim, dim);
  F_S.mul<false, false>(F, piola);
  Real constant = J ? 1./J : 0;
  sigma.mul<false, true>(F_S, F, constant);
}

/* -------------------------------------------------------------------------- */
inline void Material::rightCauchy(const Matrix<Real> & F,
				  Matrix<Real> & C) const {
  C.mul<true, false>(F, F);
}

/* -------------------------------------------------------------------------- */
inline void Material::leftCauchy(const Matrix<Real> & F,
				 Matrix<Real> & B) const {
  B.mul<false, true>(F, F);
}

/* -------------------------------------------------------------------------- */
template<UInt dim>
inline void Material::gradUToEpsilon(const Matrix<Real> & grad_u,
				     Matrix<Real> & epsilon) const {
  for (UInt i = 0; i < dim; ++i)
    for (UInt j = 0; j < dim; ++j)
      epsilon(i, j) = 0.5*(grad_u(i, j) + grad_u(j, i));
}

/* -------------------------------------------------------------------------- */
template<UInt dim>
inline void Material::gradUToGreenStrain(const Matrix<Real> & grad_u,
					 Matrix<Real> & epsilon) const {
  epsilon.mul<true, false>(grad_u, grad_u, .5);

  for (UInt i = 0; i < dim; ++i)
    for (UInt j = 0; j < dim; ++j)
      epsilon(i, j) += 0.5 * (grad_u(i, j) + grad_u(j, i));
}

/* ---------------------------------------------------------------------------*/
template<UInt dim>
inline void Material::SetCauchyStressArray(const Matrix<Real> & S_t, Matrix<Real> & Stress_vect) {

    AKANTU_DEBUG_IN();

    Stress_vect.clear();

    //UInt cauchy_matrix_size = getCauchyStressArraySize(dim);

    //see Finite ekement formulations for large deformation dynamic analysis, Bathe et al. IJNME vol 9, 1975, page 364 ^t\tau

    /*
     * 1d: [ s11 ]'
     * 2d: [ s11 s22 s12 ]'
     * 3d: [ s11 s22 s33 s23 s13 s12 ]
     */
    for (UInt i = 0; i < dim; ++i)//diagonal terms
        Stress_vect(i, 0) = S_t(i, i);

    for (UInt i = 1; i < dim; ++i)// term s12 in 2D and terms s23 s13 in 3D
        Stress_vect(dim+i-1, 0) = S_t(dim-i-1, dim-1);

    for (UInt i = 2; i < dim; ++i)//term s13 in 3D
        Stress_vect(dim+i, 0) = S_t(0, 1);

    AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt dim>
inline void Material::setCauchyStressMatrix(const Matrix<Real> & S_t, Matrix<Real> & Stress_matrix) {

    AKANTU_DEBUG_IN();

    Stress_matrix.clear();

    /// see Finite ekement formulations for large deformation dynamic analysis,
    /// Bathe et al. IJNME vol 9, 1975, page 364 ^t\tau

    for (UInt i = 0; i < dim; ++i) {
        for (UInt m = 0; m < dim; ++m) {
            for (UInt n = 0; n < dim; ++n) {
                Stress_matrix(i * dim + m, i * dim + n) = S_t(m, n);
            }
        }
    }

    //other terms from the diagonal
    /*for (UInt i = 0; i < 3 - dim; ++i) {
        Stress_matrix(dim * dim + i, dim * dim + i) = S_t(dim + i, dim + i);
    }*/


    AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<ElementType type>
inline void Material::buildElementalFieldInterpolationCoodinates(__attribute__((unused)) const Matrix<Real> & coordinates,
								 __attribute__((unused)) Matrix<Real> & coordMatrix) {
  AKANTU_DEBUG_TO_IMPLEMENT();
}

/* -------------------------------------------------------------------------- */
inline void Material::buildElementalFieldInterpolationCoodinatesLinear(const Matrix<Real> & coordinates,
								       Matrix<Real> & coordMatrix) {

  for (UInt i = 0; i < coordinates.cols(); ++i)
    coordMatrix(i, 0) = 1;
}

/* -------------------------------------------------------------------------- */
inline void Material::buildElementalFieldInterpolationCoodinatesQuadratic(const Matrix<Real> & coordinates,
									  Matrix<Real> & coordMatrix) {

  UInt nb_quadrature_points = coordMatrix.cols();

  for (UInt i = 0; i < coordinates.cols(); ++i) {
    coordMatrix(i, 0) = 1;
    for (UInt j = 1; j < nb_quadrature_points; ++j)
      coordMatrix(i, j) = coordinates(j-1, i);
  }
}

/* -------------------------------------------------------------------------- */
template<>
inline void Material::buildElementalFieldInterpolationCoodinates<_segment_2>(const Matrix<Real> & coordinates,
									     Matrix<Real> & coordMatrix) {
  buildElementalFieldInterpolationCoodinatesLinear(coordinates, coordMatrix);
}

/* -------------------------------------------------------------------------- */
template<>
inline void Material::buildElementalFieldInterpolationCoodinates<_segment_3>(const Matrix<Real> & coordinates,
									     Matrix<Real> & coordMatrix) {

  buildElementalFieldInterpolationCoodinatesQuadratic(coordinates, coordMatrix);
}

/* -------------------------------------------------------------------------- */
template<>
inline void Material::buildElementalFieldInterpolationCoodinates<_triangle_3>(const Matrix<Real> & coordinates,
									      Matrix<Real> & coordMatrix) {
  buildElementalFieldInterpolationCoodinatesLinear(coordinates, coordMatrix);
}

/* -------------------------------------------------------------------------- */
template<>
inline void Material::buildElementalFieldInterpolationCoodinates<_triangle_6>(const Matrix<Real> & coordinates,
									      Matrix<Real> & coordMatrix) {

  buildElementalFieldInterpolationCoodinatesQuadratic(coordinates, coordMatrix);
}


/* -------------------------------------------------------------------------- */
template<>
inline void Material::buildElementalFieldInterpolationCoodinates<_tetrahedron_4>(const Matrix<Real> & coordinates,
										 Matrix<Real> & coordMatrix) {
  buildElementalFieldInterpolationCoodinatesLinear(coordinates, coordMatrix);
}

/* -------------------------------------------------------------------------- */
template<>
inline void Material::buildElementalFieldInterpolationCoodinates<_tetrahedron_10>(const Matrix<Real> & coordinates,
										  Matrix<Real> & coordMatrix) {

  buildElementalFieldInterpolationCoodinatesQuadratic(coordinates, coordMatrix);
}

/**
 * @todo Write a more efficient interpolation for quadrangles by
 * dropping unnecessary quadrature points
 *
 */

/* -------------------------------------------------------------------------- */
template<>
inline void Material::buildElementalFieldInterpolationCoodinates<_quadrangle_4>(const Matrix<Real> & coordinates,
										Matrix<Real> & coordMatrix) {

  for (UInt i = 0; i < coordinates.cols(); ++i) {
    Real x = coordinates(0, i);
    Real y = coordinates(1, i);

    coordMatrix(i, 0) = 1;
    coordMatrix(i, 1) = x;
    coordMatrix(i, 2) = y;
    coordMatrix(i, 3) = x * y;
  }
}

/* -------------------------------------------------------------------------- */
template<>
inline void Material::buildElementalFieldInterpolationCoodinates<_quadrangle_8>(const Matrix<Real> & coordinates,
										Matrix<Real> & coordMatrix) {

  for (UInt i = 0; i < coordinates.cols(); ++i) {

    UInt j = 0;
    Real x = coordinates(0, i);
    Real y = coordinates(1, i);

    for (UInt e = 0; e <= 2; ++e) {
      for (UInt n = 0; n <= 2; ++n) {
	coordMatrix(i, j) = std::pow(x, e) * std::pow(y, n);
	++j;
      }
    }

  }
}

/* -------------------------------------------------------------------------- */
template<ElementType type>
inline UInt Material::getSizeElementalFieldInterpolationCoodinates(GhostType ghost_type) {
  return model->getFEEngine().getNbQuadraturePoints(type, ghost_type);
}

/* -------------------------------------------------------------------------- */
inline UInt Material::getNbDataForElements(const Array<Element> & elements,
					   SynchronizationTag tag) const {
  if(tag == _gst_smm_stress) {
    return (this->isFiniteDeformation() ? 3 : 1) * spatial_dimension * spatial_dimension *
      sizeof(Real) * this->getModel().getNbQuadraturePoints(elements);
  }
  return 0;
}

/* -------------------------------------------------------------------------- */
inline void Material::packElementData(CommunicationBuffer & buffer,
				      const Array<Element> & elements,
				      SynchronizationTag tag) const {
  if(tag == _gst_smm_stress) {
    if(this->isFiniteDeformation()) {
      packElementDataHelper(piola_kirchhoff_2, buffer, elements);
      packElementDataHelper(gradu, buffer, elements);
    }
    packElementDataHelper(stress, buffer, elements);
  }
}

/* -------------------------------------------------------------------------- */
inline void Material::unpackElementData(CommunicationBuffer & buffer,
					const Array<Element> & elements,
					SynchronizationTag tag) {
  if(tag == _gst_smm_stress) {
    if(this->isFiniteDeformation()) {
      unpackElementDataHelper(piola_kirchhoff_2, buffer, elements);
      unpackElementDataHelper(gradu, buffer, elements);
    }
    unpackElementDataHelper(stress, buffer, elements);
  }
}

/* -------------------------------------------------------------------------- */
template <typename T>
inline const T & Material::getParam(const ID & param) const {
  try {
    return get<T>(param);
  } catch (...) {
    AKANTU_EXCEPTION("No parameter " << param << " in the material " << getID());
  }
}

/* -------------------------------------------------------------------------- */
template <typename T>
inline void Material::setParam(const ID & param, T value) {
  try {
    set<T>(param, value);
  } catch(...) {
    AKANTU_EXCEPTION("No parameter " << param << " in the material " << getID());
  }
  updateInternalParameters();
}

/* -------------------------------------------------------------------------- */
template<typename T>
inline void Material::packElementDataHelper(const ElementTypeMapArray<T> & data_to_pack,
					    CommunicationBuffer & buffer,
					    const Array<Element> & elements,
					    const ID & fem_id) const {
  DataAccessor::packElementalDataHelper<T>(data_to_pack, buffer, elements, true, 
					   model->getFEEngine(fem_id));
}

/* -------------------------------------------------------------------------- */
template<typename T>
inline void Material::unpackElementDataHelper(ElementTypeMapArray<T> & data_to_unpack,
					      CommunicationBuffer & buffer,
					      const Array<Element> & elements,
					      const ID & fem_id) {
  DataAccessor::unpackElementalDataHelper<T>(data_to_unpack, buffer, elements, true,
					     model->getFEEngine(fem_id));
}

/* -------------------------------------------------------------------------- */
template<> inline void Material::registerInternal<Real>(InternalField<Real> & vect) {
  internal_vectors_real[vect.getID()] = &vect;
}

template<> inline void Material::registerInternal<UInt>(InternalField<UInt> & vect) {
  internal_vectors_uint[vect.getID()] = &vect;
}

/* -------------------------------------------------------------------------- */
template<> inline void Material::unregisterInternal<Real>(InternalField<Real> & vect) {
  internal_vectors_real.erase(vect.getID());
}

template<> inline void Material::unregisterInternal<UInt>(InternalField<UInt> & vect) {
  internal_vectors_uint.erase(vect.getID());
}

/* -------------------------------------------------------------------------- */
inline bool Material::isInternal(const ID & id, const ElementKind & element_kind) const {

  std::map<ID, InternalField<Real> *>::const_iterator internal_array =
    internal_vectors_real.find(this->getID()+":"+id);

  if (internal_array == internal_vectors_real.end())    return false;
  if (internal_array->second->getElementKind() != element_kind) return false;
  return true;
}

/* -------------------------------------------------------------------------- */

inline ElementTypeMap<UInt> Material::getInternalDataPerElem(const ID & id, const ElementKind & element_kind) const {

  std::map<ID, InternalField<Real> *>::const_iterator internal_array =
    internal_vectors_real.find(this->getID()+":"+id);

  if (internal_array == internal_vectors_real.end())  AKANTU_EXCEPTION("cannot find internal " << id);
  if (internal_array->second->getElementKind() != element_kind) AKANTU_EXCEPTION("cannot find internal " << id);
  
  InternalField<Real> & internal = *internal_array->second;
  InternalField<Real>::type_iterator it = internal.firstType(spatial_dimension, _not_ghost,element_kind);
  InternalField<Real>::type_iterator last_type = internal.lastType(spatial_dimension, _not_ghost,element_kind);


  ElementTypeMap<UInt> res;
  for(; it != last_type; ++it) {
    UInt nb_quadrature_points = 0;
    if (element_kind == _ek_regular)
      nb_quadrature_points = model->getFEEngine().getNbQuadraturePoints(*it);
#if defined(AKANTU_COHESIVE_ELEMENT)
    else if (element_kind == _ek_cohesive)
      nb_quadrature_points = model->getFEEngine("CohesiveFEEngine").getNbQuadraturePoints(*it);
#endif
    res(*it) = internal.getNbComponent() * nb_quadrature_points;
  }
    return res;
}

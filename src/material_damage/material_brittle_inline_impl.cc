/**
 * @file   material_brittle_inline_impl.cc
 *
 * @author Josué Aranda <josue.arandaruiz@epfl.ch>
 * @author Daniel Pino Muñoz <daniel.pinomunoz@epfl.ch>
 *
 * @date   Wed Feb 12 11:09:36 2014
 *
 * @brief  Implementation of the inline functions of the material brittle
 *
 * @section LICENSE
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
template<UInt spatial_dimension>
inline void
 MaterialBrittle<spatial_dimension>::computeStressOnQuad(Matrix<Real> & grad_u,
                                                         Matrix<Real> & grad_v,
                                                         Matrix<Real> & sigma,
                                                         Real & dam,
                                                         Real & sigma_equivalent,
                                                         Real & fracture_stress) {

  MaterialElastic<spatial_dimension>::computeStressOnQuad(grad_u, sigma);

  Real equiv_strain_rate=0.;
  Real volume_change_rate=grad_v.trace();
  if(spatial_dimension==2){
//    if(this->plane_stress){
//      Real e_dot_33 = this->lambda*(volume_change_rate)/(2. * this->mu - this->lambda);
//      volume_change_rate+=e_dot_33;
//      equiv_strain_rate+=2./3.*pow( e_dot_33 - volume_change_rate/3. , 2. );
//    }
//    else
      equiv_strain_rate+=2./3.*pow( volume_change_rate/3. , 2. );
  }

  for(UInt i=0; i<spatial_dimension; ++i)
    for(UInt j=0; j<spatial_dimension; ++j)
      equiv_strain_rate += 2./3. * pow( 0.5*(grad_v(i,j)+grad_v(j,i))-(i==j)*volume_change_rate/3. , 2. );

  equiv_strain_rate=sqrt(equiv_strain_rate);

  fracture_stress=S_0;
  if(equiv_strain_rate>E_0)
    //fracture_stress=A*pow(equiv_strain_rate,3)+B*pow(equiv_strain_rate,2)+C*equiv_strain_rate+D;
    fracture_stress=A;

  Vector<Real> principal_stress(spatial_dimension);
  sigma.eig(principal_stress);
  sigma_equivalent = principal_stress(0);
  for(UInt i=1; i<spatial_dimension; ++i)
    sigma_equivalent = std::max(sigma_equivalent, principal_stress(i));



  /*  Y = 0;
  for (UInt i = 0; i < spatial_dimension; ++i) {
    for (UInt j = 0; j < spatial_dimension; ++j) {
      Y += sigma(i,j) * grad_u(i,j);
    }
  }
  Y *= 0.5;

  if(damage_in_y) Y *= (1 - dam);

  if(yc_limit) Y = std::min(Y, Yc);
  */
  if(!this->is_non_local) {
    computeDamageAndStressOnQuad(sigma, dam, sigma_equivalent, fracture_stress);
  }
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
inline void
MaterialBrittle<spatial_dimension>::computeDamageAndStressOnQuad(Matrix<Real> & sigma,
                                                                Real & dam,
                                                                Real & sigma_c,
                                                                Real & fracture_stress) {
  if ( sigma_c > fracture_stress )
    dam = 1.;

  dam = std::min(dam,1.);

  sigma *= 1-dam;
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
inline UInt MaterialBrittle<spatial_dimension>::getNbDataForElements(const Array<Element> & elements,
                                                                    SynchronizationTag tag) const {
  AKANTU_DEBUG_IN();

  /*UInt size = 0;
  if(tag == _gst_smm_init_mat) {
    size += sizeof(Real) * this->getModel().getNbQuadraturePoints(elements);
  }
  */ 
  UInt size = MaterialDamage<spatial_dimension>::getNbDataForElements(elements, tag);

  AKANTU_DEBUG_OUT();
  return size;
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
inline void MaterialBrittle<spatial_dimension>::packElementData(CommunicationBuffer & buffer,
                                                               const Array<Element> & elements,
                                                               SynchronizationTag tag) const {
  AKANTU_DEBUG_IN();

  /*if(tag == _gst_smm_init_mat) {
    this->packElementDataHelper(Yd, buffer, elements);
    }*/

  MaterialDamage<spatial_dimension>::packElementData(buffer, elements, tag);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template<UInt spatial_dimension>
inline void MaterialBrittle<spatial_dimension>::unpackElementData(CommunicationBuffer & buffer,
                                                                 const Array<Element> & elements,
                                                                 SynchronizationTag tag) {
  AKANTU_DEBUG_IN();

  /*if(tag == _gst_smm_init_mat) {
    this->unpackElementDataHelper(Yd, buffer, elements);
    }*/

  MaterialDamage<spatial_dimension>::unpackElementData(buffer, elements, tag);

  AKANTU_DEBUG_OUT();
}

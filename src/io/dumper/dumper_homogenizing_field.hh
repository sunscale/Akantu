/**
 * @file   dumper_homogenizing_field.hh
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Tue Sep 02 2014
 * @date last modification: Tue Sep 02 2014
 *
 * @brief  description of field homogenizing field
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

#ifndef __AKANTU_DUMPER_HOMOGENIZING_FIELD_HH__
#define __AKANTU_DUMPER_HOMOGENIZING_FIELD_HH__
/* -------------------------------------------------------------------------- */
#include "dumper_compute.hh"
/* -------------------------------------------------------------------------- */

__BEGIN_AKANTU__
__BEGIN_AKANTU_DUMPER__

/* -------------------------------------------------------------------------- */


template <typename type> 
inline type typeConverter(const type & input, 
			  Vector<typename type::value_type> & res,UInt nb_data){

  throw;
  return input;
}

/* -------------------------------------------------------------------------- */


template <typename type> 
inline Matrix<type> typeConverter(const Matrix<type> & input, 
				  Vector<type> & res,UInt nb_data){

  Matrix<type> tmp(res.storage(),input.rows(),nb_data/input.rows());
  Matrix<type> tmp2(tmp,true);
  return tmp2;
}

/* -------------------------------------------------------------------------- */


template <typename type> 
inline Vector<type> typeConverter(const Vector<type> & input, 
				  Vector<type> & res,UInt nb_data){

  return res;
}

/* -------------------------------------------------------------------------- */

template <typename type>
class AvgHomogenizingFunctor : public ComputeFunctor<type,type> {
 
  /* ------------------------------------------------------------------------ */
  /* Typedefs                                                                 */
  /* ------------------------------------------------------------------------ */  
  
  typedef typename type::value_type value_type;

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
 
public:
 
  AvgHomogenizingFunctor(ElementTypeMap<UInt> & nb_datas){

    ElementTypeMap<UInt>::type_iterator tit = nb_datas.firstType();
    ElementTypeMap<UInt>::type_iterator end = nb_datas.lastType();

    nb_data = nb_datas(*tit);

    for (;tit != end; ++tit) if (nb_data != nb_datas(*tit)) throw;
  }
  
  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
 
  virtual type func(const type & d, Element global_index){

    Vector<value_type> res(this->nb_data);
    
    if (d.size() % this->nb_data) throw;
    UInt nb_to_average = d.size()/this->nb_data;

    value_type * ptr =  d.storage();
    for (UInt i = 0; i < nb_to_average; ++i) {
      Vector<value_type> tmp(ptr,this->nb_data);
      res += tmp;
    }
    res /= nb_to_average;
    return typeConverter(d,res,this->nb_data);
  };


  UInt getDim(){return nb_data;};
  UInt getNbComponent(UInt old_nb_comp){
    throw;
  };

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
 
  /// The size of data: i.e. the size of the vector to be returned 
  UInt nb_data;
};
/* -------------------------------------------------------------------------- */



/* -------------------------------------------------------------------------- */


class HomogenizerProxy {

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */

public:

  HomogenizerProxy(){};

public:

  inline static ComputeFunctorInterface * createHomogenizer(Field & field);

  template <typename T>
  inline ComputeFunctorInterface * connectToField(T * field){
    ElementTypeMap<UInt> nb_components = field->getNbComponents();

    typedef typename T::types::return_type ret_type;
    return this->instantiateHomogenizer<ret_type>(nb_components);
  }

  template <typename ret_type>
  inline ComputeFunctorInterface * 
  instantiateHomogenizer(ElementTypeMap<UInt> & nb_components);


};

/* -------------------------------------------------------------------------- */

template <typename ret_type>
inline ComputeFunctorInterface * HomogenizerProxy
::instantiateHomogenizer(ElementTypeMap<UInt> &  nb_components){

  typedef dumper::AvgHomogenizingFunctor<ret_type> Homogenizer;
  Homogenizer * foo = new Homogenizer(nb_components);
  return foo;
}

template <>
inline ComputeFunctorInterface * HomogenizerProxy
::instantiateHomogenizer<Vector<iohelper::ElemType> > (ElementTypeMap<UInt> & nb_components){
  throw;
  return NULL;
}


/* -------------------------------------------------------------------------- */


/// for connection to a FieldCompute
template <typename SubFieldCompute, typename return_type>
inline ComputeFunctorInterface * FieldCompute<SubFieldCompute,return_type>
::connect(HomogenizerProxy & proxy){

  return proxy.connectToField(this);
}
/* -------------------------------------------------------------------------- */

inline ComputeFunctorInterface * HomogenizerProxy::createHomogenizer(Field & field){
    
  HomogenizerProxy homogenizer_proxy;
  return field.connect(homogenizer_proxy);
}

/* -------------------------------------------------------------------------- */


// inline ComputeFunctorInterface & createHomogenizer(Field & field){
  
//   HomogenizerProxy::createHomogenizer(field);
//   throw;
//   ComputeFunctorInterface * ptr = NULL;
//   return *ptr;
// }

// /* -------------------------------------------------------------------------- */



__END_AKANTU_DUMPER__
__END_AKANTU__


#endif /* __AKANTU_DUMPER_HOMOGENIZING_FIELD_HH__ */

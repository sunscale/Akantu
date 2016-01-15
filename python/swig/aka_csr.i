/**
 * @file   aka_csr.i
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Mon Aug 03 2015
 * @date last modification: Mon Nov 16 2015
 *
 * @brief  
 *
 * @section LICENSE
 *
 * Copyright (©) 2015 EPFL (Ecole Polytechnique Fédérale de Lausanne) Laboratory
 * (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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

%{
  #include "aka_csr.hh"
%}

namespace akantu {
  %ignore CSR::begin;
}

%inline %{
namespace akantu {
  template <typename T>
  class CSRIterator{

  public:
  CSRIterator(CSR<T> & csr,UInt row) {
    this->it =  csr.begin(row);
    this->end =  csr.end(row);
  };

  ~CSRIterator(){
  };

  T & __next_cpp(){
    if (this->it == this->end) AKANTU_SILENT_EXCEPTION("StopIteration");
    T & ref = *(this->it);
    ++this->it;
    return ref;
  }

  private:

  typename CSR<T>::iterator it;
  typename CSR<T>::iterator end;
  };
}
%}

%extend akantu::CSRIterator<akantu::Element>
{
  %insert("python") %{
    def __iter__(self):
       return self

    def __next__(self):
       try:
         return self.__next_cpp()
       except Exception as e:
         raise StopIteration


    def next(self):
       return self.__next__()

%}
}

%extend akantu::CSR<akantu::Element>
{
  akantu::CSRIterator<akantu::Element> row(akantu::UInt row){
    return akantu::CSRIterator<akantu::Element>(*$self,row);
  }
}

%include "aka_csr.hh"
namespace akantu {
  %template (CSRUInt) CSR<UInt>;
  %template (CSRElement) CSR<Element>;
  %template (CSRIteratorElement) CSRIterator<Element>;
 }

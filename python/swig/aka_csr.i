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

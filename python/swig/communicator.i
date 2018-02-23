namespace akantu {
%ignore CommunicationBuffer::operator=;
%template(DataAccessorElement) DataAccessor<Element>;
}

%include "data_accessor.hh"

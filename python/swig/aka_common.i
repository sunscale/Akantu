%{
  #include "aka_common.hh"
%}

namespace akantu {
  %ignore getStaticParser;
  %ignore getUserParser;
  %ignore initialize(int & argc, char ** & argv);

  %ignore initialize(const std::string & input_file, int & argc, char ** & argv);
}


%inline %{
  namespace akantu {
    void initialize(const std::string & input_file) {
      int argc = 0;
      char ** argv = NULL;
      initialize(input_file, argc, argv);
    }
    void initialize() {
      int argc = 0;
      char ** argv = NULL;
      initialize(argc, argv);
    }  
  }
  %}


%include "aka_common.hh"

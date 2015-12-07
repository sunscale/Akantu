%{
#include "cohesive_element_inserter.hh"
#include "solid_mechanics_model_cohesive.hh"
%}

namespace akantu {
  %ignore SolidMechanicsModelCohesive::initFacetFilter;
  %ignore SolidMechanicsModelCohesive::initParallel;
  %ignore CohesiveElementInserter::initParallel;
}

%include "cohesive_element_inserter.hh"
%include "solid_mechanics_model_cohesive.hh"


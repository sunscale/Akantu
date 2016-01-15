%{
#include "cohesive_element_inserter.hh"
#include "solid_mechanics_model_cohesive.hh"
#include "material_cohesive.hh"
%}

namespace akantu {
  %ignore SolidMechanicsModelCohesive::initFacetFilter;
  %ignore SolidMechanicsModelCohesive::initParallel;
  %ignore CohesiveElementInserter::initParallel;
}

%extend akantu::SolidMechanicsModelCohesive {

  Array<Real> & getRealInternalCohesiveField(const std::string field_name) {
    akantu::Mesh & mesh = $self->getMesh();
    akantu::ElementType type = *(mesh.firstType(mesh.getSpatialDimension(), akantu::_not_ghost, akantu::_ek_cohesive));
    return ($self->flattenInternal(field_name,akantu::_ek_cohesive, akantu::_not_ghost))(type);
  }
}

%include "cohesive_element_inserter.hh"
%include "solid_mechanics_model_cohesive.hh"


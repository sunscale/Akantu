/**
 * @file   material_python.cc
 *
 * @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
 *
 */

/* -------------------------------------------------------------------------- */
#include "material_python.hh"
#include "solid_mechanics_model.hh"

__BEGIN_AKANTU__

/* -------------------------------------------------------------------------- */
MaterialPython::MaterialPython(SolidMechanicsModel & model,
			       PyObject * obj,
			       const ID & id)  :
  
  Material(model, id),
  PythonFunctor(obj){
  AKANTU_DEBUG_IN();

  this->registerInternals();
  
  std::vector<std::string> param_names = this->callFunctor<std::vector<std::string> >("registerParam");

  this->local_params.resize(param_names.size());

  for (UInt i = 0; i < param_names.size(); ++i) {
    std::stringstream sstr;
    sstr << "PythonParameter" << i; 
    this->registerParam(param_names[i],local_params[i],0., _pat_parsable,sstr.str() );
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */

void MaterialPython::registerInternals() {

  std::vector<std::string> internal_names = this->callFunctor<std::vector<std::string> >("registerInternals");

  this->internals.resize(internal_names.size());

  for (UInt i = 0; i < internal_names.size(); ++i) {
    std::stringstream sstr;
    sstr << "PythonInternal" << i; 
    this->internals[i] = new InternalField<Real>(internal_names[i],*this);
    this->internals[i]->initialize(1);
  }

  
}
/* -------------------------------------------------------------------------- */


void MaterialPython::initMaterial() {
  AKANTU_DEBUG_IN();

  Material::initMaterial();

  // initInternalArray(this->damage, 1);
  // resizeInternalArray(this->damage);


  // lambda = nu * E / ((1 + nu) * (1 - 2*nu));
  // mu     = E / (2 * (1 + nu));

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MaterialPython::computeStress(ElementType el_type, GhostType ghost_type) {
  AKANTU_DEBUG_IN();

  typedef Array<Real>::iterator<Real> it_type;
  std::vector<it_type> its;
  for (auto i: this->internals) {
    its.push_back((*i)(el_type,ghost_type).begin());
  }
  
  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type);

  computeStress(grad_u, sigma, its);

  for (auto b: its) ++(*b);

  MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <typename it_type>
void  MaterialPython::computeStress( Matrix<Real> & grad_u,
				    Matrix<Real> & sigma,
				    std::vector<it_type> & internal_iterators){

  std::vector<Real> inputs;
  for (auto i: internal_iterators){
    inputs.push_back(*i);
  }
  this->callFunctor<void>("computeStress",grad_u,sigma,inputs);
  
  for (UInt i = 0; i < inputs.size(); ++i) {
    *internal_iterators[i] = inputs[i];
   }
}

/* -------------------------------------------------------------------------- */
void MaterialPython::computeTangentModuli(const ElementType & el_type,
					  Array<Real> & tangent_matrix,
					  GhostType ghost_type){

  this->callFunctor<void>("computeTangentModuli",el_type,tangent_matrix,ghost_type);
}

/* -------------------------------------------------------------------------- */

Real MaterialPython::getPushWaveSpeed(const Element & element) const{
  return this->callFunctor<Real>("getPushWaveSpeed");
}

__END_AKANTU__

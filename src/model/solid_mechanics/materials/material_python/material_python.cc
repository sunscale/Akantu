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
  Material(model, id) {
  AKANTU_DEBUG_IN();

  // this->registerParam("E"           , E           , 0.   , _pat_parsable, "Young's modulus"        );
  // this->registerParam("nu"          , nu          , 0.5  , _pat_parsable, "Poisson's ratio"        );
  // this->registerParam("lambda"      , lambda             , _pat_readable, "First Lamé coefficient" );
  // this->registerParam("mu"          , mu                 , _pat_readable, "Second Lamé coefficient");
  // this->registerParam("Yd"          , Yd          ,   50., _pat_parsmod);
  // this->registerParam("Sd"          , Sd          , 5000., _pat_parsmod);

  AKANTU_DEBUG_OUT();
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
// void MaterialPython::computeStress(ElementType el_type, GhostType ghost_type) {
//   AKANTU_DEBUG_IN();

//   // Array<Real>::iterator<Real> dam_it = damage(el_type, ghost_type).begin();

//   // MATERIAL_STRESS_QUADRATURE_POINT_LOOP_BEGIN(el_type, ghost_type);
//   // Real & dam = *dam_it;
//   // computeStress(grad_u, sigma, dam);
//   // ++dam_it;
//   // MATERIAL_STRESS_QUADRATURE_POINT_LOOP_END;

//   AKANTU_DEBUG_OUT();
// }

/* -------------------------------------------------------------------------- */

__END_AKANTU__

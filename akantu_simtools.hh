// ast common
#include "ast_common.hh"
#include "synchronized_array.hh"
#include "parameter_reader.hh"
#include "manual_restart.hh"

// functions
#include "boundary_functions.hh"
#include "node_filter.hh"

// boundary conditions
#include "force_based_dirichlet.hh"
#include "spring_bc.hh"

// ntn/ntrf contact
#include "ntn_base_contact.hh"
#include "ntn_contact.hh"
#include "ntrf_contact.hh"

// ntn/ntrf friction
#include "ntn_base_friction.hh"
#include "ntn_friction.hh"
#include "ntrf_friction.hh"

// friction regularisations
#include "ntn_fricreg_no_regularisation.hh"
#include "ntn_fricreg_rubin_ampuero.hh"
#include "ntn_fricreg_simplified_prakash_clifton.hh"

// friction laws
#include "ntn_friclaw_coulomb.hh"
#include "ntn_friclaw_linear_slip_weakening.hh"

// initiation of friction
#include "ntn_initiation_function.hh"

/*
#include "ntn_friction_coulomb.hh"
#include "ntn_friction_linear_slip_weakening.hh"
#include "ntrf_friction_coulomb.hh"
#include "ntrf_friction_regularized_coulomb.hh"
#include "ntrf_friction_linear_slip_weakening.hh"
#include "ntrf_friction_mathilde.hh"
*/

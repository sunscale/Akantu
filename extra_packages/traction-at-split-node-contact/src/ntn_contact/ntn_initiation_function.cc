/**
 * @file   ntn_initiation_function.cc
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 *
 *
 * @brief  implementation of initializing ntn and ntrf friction
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2012, 2014 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 */

/* -------------------------------------------------------------------------- */
// simtools
#include "ntn_initiation_function.hh"
#include "mIIasym_contact.hh"
#include "ntn_friction.hh"
#include "ntrf_friction.hh"

// friction regularisations
#include "ntn_fricreg_rubin_ampuero.hh"
#include "ntn_fricreg_simplified_prakash_clifton.hh"

// friction laws
#include "ntn_friclaw_linear_cohesive.hh"
#include "ntn_friclaw_linear_slip_weakening.hh"
#include "ntn_friclaw_linear_slip_weakening_no_healing.hh"

namespace akantu {

NTNBaseFriction * initializeNTNFriction(NTNBaseContact * contact,
                                        ParameterReader & data) {
  AKANTU_DEBUG_IN();

  const std::string & friction_law = data.get<std::string>("friction_law");
  const std::string & friction_reg =
      data.get<std::string>("friction_regularisation");

  NTNBaseFriction * friction;

  bool is_ntn_contact = true;
  if (dynamic_cast<NTRFContact *>(contact) != NULL ||
      dynamic_cast<MIIASYMContact *>(contact) != NULL) {
    is_ntn_contact = false;
  }

  if (friction_law == "coulomb") {
    if (friction_reg == "no_regularisation") {
      if (is_ntn_contact)
        friction =
            new NTNFriction<NTNFricLawCoulomb, NTNFricRegNoRegularisation>(
                contact);
      else
        friction =
            new NTRFFriction<NTNFricLawCoulomb, NTNFricRegNoRegularisation>(
                contact);
    } else if (friction_reg == "rubin_ampuero") {
      if (is_ntn_contact)
        friction =
            new NTNFriction<NTNFricLawCoulomb, NTNFricRegRubinAmpuero>(contact);
      else
        friction = new NTRFFriction<NTNFricLawCoulomb, NTNFricRegRubinAmpuero>(
            contact);

      friction->setMixed<SynchronizedArray<Real>>("t_star",
                                                  data.get<Real>("t_star"));
    } else if (friction_reg == "simplified_prakash_clifton") {
      if (is_ntn_contact)
        friction = new NTNFriction<NTNFricLawCoulomb,
                                   NTNFricRegSimplifiedPrakashClifton>(contact);
      else
        friction =
            new NTRFFriction<NTNFricLawCoulomb,
                             NTNFricRegSimplifiedPrakashClifton>(contact);

      friction->setMixed<SynchronizedArray<Real>>("t_star",
                                                  data.get<Real>("t_star"));
    } else {
      AKANTU_ERROR("Do not know the following friction regularisation: "
                   << friction_reg);
    }

    friction->setMixed<SynchronizedArray<Real>>("mu_s", data.get<Real>("mu_s"));
  }

  // Friction Law: Linear Slip Weakening
  else if (friction_law == "linear_slip_weakening") {
    if (friction_reg == "no_regularisation") {
      if (is_ntn_contact)
        friction = new NTNFriction<NTNFricLawLinearSlipWeakening,
                                   NTNFricRegNoRegularisation>(contact);
      else
        friction = new NTRFFriction<NTNFricLawLinearSlipWeakening,
                                    NTNFricRegNoRegularisation>(contact);
    } else if (friction_reg == "rubin_ampuero") {
      if (is_ntn_contact)
        friction = new NTNFriction<NTNFricLawLinearSlipWeakening,
                                   NTNFricRegRubinAmpuero>(contact);
      else
        friction = new NTRFFriction<NTNFricLawLinearSlipWeakening,
                                    NTNFricRegRubinAmpuero>(contact);

      friction->setMixed<SynchronizedArray<Real>>("t_star",
                                                  data.get<Real>("t_star"));
    } else if (friction_reg == "simplified_prakash_clifton") {
      if (is_ntn_contact)
        friction = new NTNFriction<NTNFricLawLinearSlipWeakening,
                                   NTNFricRegSimplifiedPrakashClifton>(contact);
      else
        friction =
            new NTRFFriction<NTNFricLawLinearSlipWeakening,
                             NTNFricRegSimplifiedPrakashClifton>(contact);

      friction->setMixed<SynchronizedArray<Real>>("t_star",
                                                  data.get<Real>("t_star"));
    } else {
      AKANTU_ERROR("Do not know the following friction regularisation: "
                   << friction_reg);
    }

    friction->setMixed<SynchronizedArray<Real>>("mu_s", data.get<Real>("mu_s"));
    friction->setMixed<SynchronizedArray<Real>>("mu_k", data.get<Real>("mu_k"));
    friction->setMixed<SynchronizedArray<Real>>("d_c", data.get<Real>("d_c"));
  }

  // Friction Law: Linear Slip Weakening No Healing
  else if (friction_law == "linear_slip_weakening_no_healing") {
    if (friction_reg == "no_regularisation") {
      if (is_ntn_contact)
        friction = new NTNFriction<NTNFricLawLinearSlipWeakeningNoHealing,
                                   NTNFricRegNoRegularisation>(contact);
      else
        friction = new NTRFFriction<NTNFricLawLinearSlipWeakeningNoHealing,
                                    NTNFricRegNoRegularisation>(contact);
    } else if (friction_reg == "rubin_ampuero") {
      if (is_ntn_contact)
        friction = new NTNFriction<NTNFricLawLinearSlipWeakeningNoHealing,
                                   NTNFricRegRubinAmpuero>(contact);
      else
        friction = new NTRFFriction<NTNFricLawLinearSlipWeakeningNoHealing,
                                    NTNFricRegRubinAmpuero>(contact);

      friction->setMixed<SynchronizedArray<Real>>("t_star",
                                                  data.get<Real>("t_star"));
    } else if (friction_reg == "simplified_prakash_clifton") {
      if (is_ntn_contact)
        friction = new NTNFriction<NTNFricLawLinearSlipWeakeningNoHealing,
                                   NTNFricRegSimplifiedPrakashClifton>(contact);
      else
        friction =
            new NTRFFriction<NTNFricLawLinearSlipWeakeningNoHealing,
                             NTNFricRegSimplifiedPrakashClifton>(contact);

      friction->setMixed<SynchronizedArray<Real>>("t_star",
                                                  data.get<Real>("t_star"));
    } else {
      AKANTU_ERROR("Do not know the following friction regularisation: "
                   << friction_reg);
    }

    friction->setMixed<SynchronizedArray<Real>>("mu_s", data.get<Real>("mu_s"));
    friction->setMixed<SynchronizedArray<Real>>("mu_k", data.get<Real>("mu_k"));
    friction->setMixed<SynchronizedArray<Real>>("d_c", data.get<Real>("d_c"));
  }

  // Friction Law: Linear Cohesive
  else if (friction_law == "linear_cohesive") {
    if (friction_reg == "no_regularisation") {
      if (is_ntn_contact)
        friction = new NTNFriction<NTNFricLawLinearCohesive,
                                   NTNFricRegNoRegularisation>(contact);
      else
        friction = new NTRFFriction<NTNFricLawLinearCohesive,
                                    NTNFricRegNoRegularisation>(contact);
    } else if (friction_reg == "rubin_ampuero") {
      if (is_ntn_contact)
        friction =
            new NTNFriction<NTNFricLawLinearCohesive, NTNFricRegRubinAmpuero>(
                contact);
      else
        friction =
            new NTRFFriction<NTNFricLawLinearCohesive, NTNFricRegRubinAmpuero>(
                contact);

      friction->setMixed<SynchronizedArray<Real>>("t_star",
                                                  data.get<Real>("t_star"));
    } else if (friction_reg == "simplified_prakash_clifton") {
      if (is_ntn_contact)
        friction = new NTNFriction<NTNFricLawLinearCohesive,
                                   NTNFricRegSimplifiedPrakashClifton>(contact);
      else
        friction =
            new NTRFFriction<NTNFricLawLinearCohesive,
                             NTNFricRegSimplifiedPrakashClifton>(contact);

      friction->setMixed<SynchronizedArray<Real>>("t_star",
                                                  data.get<Real>("t_star"));
    } else {
      AKANTU_ERROR("Do not know the following friction regularisation: "
                   << friction_reg);
    }

    friction->setMixed<SynchronizedArray<Real>>("G_c", data.get<Real>("G_c"));
    friction->setMixed<SynchronizedArray<Real>>("tau_c",
                                                data.get<Real>("tau_c"));
    friction->setMixed<SynchronizedArray<Real>>("tau_r",
                                                data.get<Real>("tau_r"));
  }

  else {
    AKANTU_ERROR("Do not know the following friction law: " << friction_law);
  }

  AKANTU_DEBUG_OUT();
  return friction;
}

/* -------------------------------------------------------------------------- */
NTNBaseFriction * initializeNTNFriction(NTNBaseContact * contact) {
  AKANTU_DEBUG_IN();
  std::pair<Parser::const_section_iterator, Parser::const_section_iterator>
      sub_sect = getStaticParser().getSubSections(_st_friction);

  Parser::const_section_iterator it = sub_sect.first;
  const ParserSection & section = *it;

  std::string friction_law = section.getName();
  std::string friction_reg = section.getOption();
  if (friction_reg == "") {
    std::string friction_reg = "no_regularisation";
  }

  NTNBaseFriction * friction =
      initializeNTNFriction(contact, friction_law, friction_reg);

  friction->parseSection(section);

  if (++it != sub_sect.second) {
    AKANTU_DEBUG_WARNING("There were several friction sections in input file. "
                         << "Only first one was used and all others ignored.");
  }

  AKANTU_DEBUG_OUT();
  return friction;
}

/* -------------------------------------------------------------------------- */
NTNBaseFriction * initializeNTNFriction(NTNBaseContact * contact,
                                        const std::string & friction_law,
                                        const std::string & friction_reg) {
  AKANTU_DEBUG_IN();

  NTNBaseFriction * friction;

  // check whether is is node-to-rigid-flat contact or mIIasym (which is also
  // ntrf)
  bool is_ntn_contact = true;
  if (dynamic_cast<NTRFContact *>(contact) != NULL ||
      dynamic_cast<MIIASYMContact *>(contact) != NULL) {
    is_ntn_contact = false;
  }

  if (friction_law == "coulomb") {
    if (friction_reg == "no_regularisation") {
      if (is_ntn_contact)
        friction =
            new NTNFriction<NTNFricLawCoulomb, NTNFricRegNoRegularisation>(
                contact);
      else
        friction =
            new NTRFFriction<NTNFricLawCoulomb, NTNFricRegNoRegularisation>(
                contact);
    } else if (friction_reg == "rubin_ampuero") {
      if (is_ntn_contact)
        friction =
            new NTNFriction<NTNFricLawCoulomb, NTNFricRegRubinAmpuero>(contact);
      else
        friction = new NTRFFriction<NTNFricLawCoulomb, NTNFricRegRubinAmpuero>(
            contact);
    } else if (friction_reg == "simplified_prakash_clifton") {
      if (is_ntn_contact)
        friction = new NTNFriction<NTNFricLawCoulomb,
                                   NTNFricRegSimplifiedPrakashClifton>(contact);
      else
        friction =
            new NTRFFriction<NTNFricLawCoulomb,
                             NTNFricRegSimplifiedPrakashClifton>(contact);
    } else {
      AKANTU_ERROR("Do not know the following friction regularisation: "
                   << friction_reg);
    }
  }

  // Friction Law: Linear Slip Weakening
  else if (friction_law == "linear_slip_weakening") {
    if (friction_reg == "no_regularisation") {
      if (is_ntn_contact)
        friction = new NTNFriction<NTNFricLawLinearSlipWeakening,
                                   NTNFricRegNoRegularisation>(contact);
      else
        friction = new NTRFFriction<NTNFricLawLinearSlipWeakening,
                                    NTNFricRegNoRegularisation>(contact);
    } else if (friction_reg == "rubin_ampuero") {
      if (is_ntn_contact)
        friction = new NTNFriction<NTNFricLawLinearSlipWeakening,
                                   NTNFricRegRubinAmpuero>(contact);
      else
        friction = new NTRFFriction<NTNFricLawLinearSlipWeakening,
                                    NTNFricRegRubinAmpuero>(contact);
    } else if (friction_reg == "simplified_prakash_clifton") {
      if (is_ntn_contact)
        friction = new NTNFriction<NTNFricLawLinearSlipWeakening,
                                   NTNFricRegSimplifiedPrakashClifton>(contact);
      else
        friction =
            new NTRFFriction<NTNFricLawLinearSlipWeakening,
                             NTNFricRegSimplifiedPrakashClifton>(contact);
    } else {
      AKANTU_ERROR("Do not know the following friction regularisation: "
                   << friction_reg);
    }
  }

  // Friction Law: Linear Slip Weakening No Healing
  else if (friction_law == "linear_slip_weakening_no_healing") {
    if (friction_reg == "no_regularisation") {
      if (is_ntn_contact)
        friction = new NTNFriction<NTNFricLawLinearSlipWeakeningNoHealing,
                                   NTNFricRegNoRegularisation>(contact);
      else
        friction = new NTRFFriction<NTNFricLawLinearSlipWeakeningNoHealing,
                                    NTNFricRegNoRegularisation>(contact);
    } else if (friction_reg == "rubin_ampuero") {
      if (is_ntn_contact)
        friction = new NTNFriction<NTNFricLawLinearSlipWeakeningNoHealing,
                                   NTNFricRegRubinAmpuero>(contact);
      else
        friction = new NTRFFriction<NTNFricLawLinearSlipWeakeningNoHealing,
                                    NTNFricRegRubinAmpuero>(contact);
    } else if (friction_reg == "simplified_prakash_clifton") {
      if (is_ntn_contact)
        friction = new NTNFriction<NTNFricLawLinearSlipWeakeningNoHealing,
                                   NTNFricRegSimplifiedPrakashClifton>(contact);
      else
        friction =
            new NTRFFriction<NTNFricLawLinearSlipWeakeningNoHealing,
                             NTNFricRegSimplifiedPrakashClifton>(contact);
    } else {
      AKANTU_ERROR("Do not know the following friction regularisation: "
                   << friction_reg);
    }
  }

  // Friction Law: Linear Cohesive
  else if (friction_law == "linear_cohesive") {
    if (friction_reg == "no_regularisation") {
      if (is_ntn_contact)
        friction = new NTNFriction<NTNFricLawLinearCohesive,
                                   NTNFricRegNoRegularisation>(contact);
      else
        friction = new NTRFFriction<NTNFricLawLinearCohesive,
                                    NTNFricRegNoRegularisation>(contact);
    } else if (friction_reg == "rubin_ampuero") {
      if (is_ntn_contact)
        friction =
            new NTNFriction<NTNFricLawLinearCohesive, NTNFricRegRubinAmpuero>(
                contact);
      else
        friction =
            new NTRFFriction<NTNFricLawLinearCohesive, NTNFricRegRubinAmpuero>(
                contact);
    } else if (friction_reg == "simplified_prakash_clifton") {
      if (is_ntn_contact)
        friction = new NTNFriction<NTNFricLawLinearCohesive,
                                   NTNFricRegSimplifiedPrakashClifton>(contact);
      else
        friction =
            new NTRFFriction<NTNFricLawLinearCohesive,
                             NTNFricRegSimplifiedPrakashClifton>(contact);
    } else {
      AKANTU_ERROR("Do not know the following friction regularisation: "
                   << friction_reg);
    }
  } else {
    AKANTU_ERROR("Do not know the following friction law: " << friction_law);
  }

  AKANTU_DEBUG_OUT();
  return friction;
}

} // namespace akantu

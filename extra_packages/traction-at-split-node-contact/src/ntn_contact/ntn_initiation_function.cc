/**
 * @file   ntn_initiation_function.cc
 *
 * @author David Simon Kammer <david.kammer@epfl.ch>
 *
 * @date creation: Tue Dec 02 2014
 * @date last modification: Fri Feb 23 2018
 *
 * @brief  implementation of initializing ntn and ntrf friction
 *
 *
 * Copyright (©) 2015-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
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

#include "aka_factory.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
std::unique_ptr<NTNBaseFriction>
initializeNTNFriction(NTNBaseContact & contact) {
  AKANTU_DEBUG_IN();
  auto sub_sect = getStaticParser().getSubSections(ParserType::_friction);
  auto it = sub_sect.first;
  const ParserSection & section = *it;

  std::string friction_law = section.getName();
  std::string friction_reg = section.getOption("no_regularisation");

  std::unique_ptr<NTNBaseFriction> friction =
      initializeNTNFriction(contact, friction_law, friction_reg);

  friction->parseSection(section);

  if (++it != sub_sect.second) {
    AKANTU_DEBUG_WARNING("There were several friction sections in input file. "
                         << "Only first one was used and all others ignored.");
  }

  AKANTU_DEBUG_OUT();
  return friction;
}

namespace {
  using NTNFactory =
      Factory<NTNBaseFriction, std::tuple<bool, ID, ID>, NTNBaseContact &>;

  // std::ostream & operator<<(std::ostream & stream,
  //                           const std::tuple<bool, ID, ID> & tuple) {
  //   stream << "[" << std::get<0>(tuple) << ", " << std::get<1>(tuple) << ", "
  //          << std::get<2>(tuple) << ", "
  //          << "]" << std::endl;
  //   return stream;
  // }

  template <bool is_ntn, template <class> class FrictionLaw, class FrictionReg>
  bool registerFriction(const ID & friction_law, const ID & friction_reg) {
    NTNFactory::getInstance().registerAllocator(
        std::make_tuple(is_ntn, friction_law, friction_reg),
        [](NTNBaseContact & contact) -> std::unique_ptr<NTNBaseFriction> {
          return std::make_unique<
              std::conditional_t<is_ntn, NTNFriction<FrictionLaw, FrictionReg>,
                                 NTRFFriction<FrictionLaw, FrictionReg>>>(
              contact);
        });
    return true;
  }

  template <template <class> class FrictionLaw, class FrictionReg>
  bool registerFrictionNTNandNTRF(const ID & friction_law,
                                  const ID & friction_reg) {
    registerFriction<true, FrictionLaw, FrictionReg>(friction_law,
                                                     friction_reg);
    registerFriction<false, FrictionLaw, FrictionReg>(friction_law,
                                                      friction_reg);
    return true;
  }

  template <template <class> class FrictionLaw>
  bool registerFrictionRegs(const ID & friction_law) {
    registerFrictionNTNandNTRF<FrictionLaw, NTNFricRegRubinAmpuero>(
        friction_law, "no_regularisation");
    registerFrictionNTNandNTRF<FrictionLaw, NTNFricRegRubinAmpuero>(
        friction_law, "rubin_ampuero");
    registerFrictionNTNandNTRF<FrictionLaw, NTNFricRegSimplifiedPrakashClifton>(
        friction_law, "simplified_prakash_clifton");
    return true;
  }

  bool registerFrictionLaws() {
    registerFrictionRegs<NTNFricLawCoulomb>("coulomb");
    registerFrictionRegs<NTNFricLawLinearSlipWeakening>(
        "linear_slip_weakening");
    registerFrictionRegs<NTNFricLawLinearSlipWeakeningNoHealing>(
        "linear_slip_weakening_no_healing");
    registerFrictionRegs<NTNFricLawLinearCohesive>("linear_cohesive");
    return true;
  }

  static bool _ = registerFrictionLaws();
} // namespace

/* -------------------------------------------------------------------------- */
std::unique_ptr<NTNBaseFriction>
initializeNTNFriction(NTNBaseContact & contact,
                      const std::string & friction_law,
                      const std::string & friction_reg) {
  bool is_ntn_contact = contact.isNTNContact();
  return NTNFactory::getInstance().allocate(
      std::make_tuple(is_ntn_contact, friction_law, friction_reg), contact);
}

} // namespace akantu

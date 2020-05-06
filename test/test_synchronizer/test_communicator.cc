/**
 * @file   test_communicator.cc
 *
 * @author Nicolas Richart
 *
 * @date creation  Thu Feb 21 2019
 *
 * @brief A Documented file.
 *
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as  published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */
/* -------------------------------------------------------------------------- */
#include <aka_iterators.hh>
#include <communication_tag.hh>
#include <communicator.hh>
/* -------------------------------------------------------------------------- */
#include <gtest/gtest.h>
#include <random>
/* -------------------------------------------------------------------------- */

using namespace akantu;

TEST(Communicator, Bcast) {
  auto r = 0xdeadbeef;

  auto & c = Communicator::getStaticCommunicator();
  c.broadcast(r);

  EXPECT_EQ(r, 0xdeadbeef);
}

TEST(Communicator, ReceiveAny) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<> dis(1, 10);
  std::vector<CommunicationRequest> reqs;

  auto & c = Communicator::getStaticCommunicator();
  auto && rank = c.whoAmI();
  auto && size = c.getNbProc();

  for (auto n : arange(100)) {
    AKANTU_DEBUG_INFO("ROUND " << n);
    auto tag = Tag::genTag(0, 1, 0);

    if (rank == 0) {
      std::vector<int> sends(size - 1);
      for (auto & s : sends) {
        s = dis(gen);
      }

      c.broadcast(sends);
      AKANTU_DEBUG_INFO("Messages " << [&]() {
        std::string msgs;
        for (auto s : enumerate(sends)) {
          if (std::get<0>(s) != 0)
            msgs += ", ";
          msgs += std::to_string(std::get<0>(s) + 1) + ": " +
                  std::to_string(std::get<1>(s));
        }
        return msgs;
      }());

      int nb_recvs = 0;
      for (auto && data : enumerate(sends)) {
        auto & send = std::get<1>(data);
        int p = std::get<0>(data) + 1;

        if (send > 5) {
          reqs.push_back(
              c.asyncSend(send, p, tag, CommunicationMode::_synchronous));
        }

        if (p <= send) {
          ++nb_recvs;
        }
      }

      c.receiveAnyNumber<int>(reqs,
                              [&](auto && proc, auto && msg) {
                                EXPECT_EQ(msg[0], sends[proc - 1] + 100 * proc);
                                EXPECT_LE(proc, sends[proc - 1]);
                                --nb_recvs;
                              },
                              tag);
      EXPECT_EQ(nb_recvs, 0);
    } else {
      std::vector<int> recv(size - 1);
      c.broadcast(recv);

      AKANTU_DEBUG_INFO("Messages " << [&]() {
        std::string msgs;
        for (auto s : enumerate(recv)) {
          if (std::get<0>(s) != 0)
            msgs += ", ";
          msgs += std::to_string(std::get<0>(s) + 1) + ": " +
                  std::to_string(std::get<1>(s));
        }
        return msgs;
      }());
      auto send = recv[rank - 1] + 100 * rank;
      if (rank <= recv[rank - 1]) {
        reqs.push_back(
            c.asyncSend(send, 0, tag, CommunicationMode::_synchronous));
      }

      bool has_recv = false;
      c.receiveAnyNumber<int>(reqs,
                              [&](auto && proc, auto && msg) {
                                EXPECT_EQ(msg[0], recv[rank - 1]);
                                EXPECT_EQ(proc, 0);
                                has_recv = true;
                              },
                              tag);

      bool should_recv = (recv[rank - 1] > 5);
      EXPECT_EQ(has_recv, should_recv);
    }
    reqs.clear();
  }
}

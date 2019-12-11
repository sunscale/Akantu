/* -------------------------------------------------------------------------- */
#include "aka_str_hash.hh"
#include "aka_tuple_tools.hh"
/* -------------------------------------------------------------------------- */
#include <gtest/gtest.h>
/* -------------------------------------------------------------------------- */

using namespace aka;
using namespace aka::hash;

TEST(STRHash, CRC32) {
  auto crc = "stack-overflow"_crc32;

  EXPECT_EQ(crc, 0x335CC04A);
}

TEST(STRHash, fnv1a32) {
  auto hash = "cyclist"_fnv1a32;

  EXPECT_EQ(hash, 0xc441a122);
}

TEST(STRHash, fnv1a64) {
  auto hash = "cyclist"_fnv1a64;

  EXPECT_EQ(hash, 0x4e376afdf5d98282);
}

TEST(STRHash, Litteral) {
  auto h = "cyclist"_h;
  EXPECT_EQ(h, 0x4e376afdf5d98282);
}

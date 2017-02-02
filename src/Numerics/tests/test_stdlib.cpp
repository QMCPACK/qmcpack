
#include "catch.hpp"

// Put tests for src/config/stdlib functions here

// Ensure that the backup "round" function is always defined
#include <config/stdlib/math.h>


namespace qmcplusplus
{

TEST_CASE("stdlib round", "[numerics]")
{
  REQUIRE(round(0.4f) == 0.0f);
  REQUIRE(round(0.6f) == 1.0f);
  REQUIRE(round(1.2f) == 1.0f);
  REQUIRE(round(1.8f) == 2.0f);
  REQUIRE(round(-1.4f) == -1.0f);
}

}

//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//
//////////////////////////////////////////////////////////////////////////////////////

#include "catch.hpp"

#include "Numerics/RotationMatrix3D.h"
#include <string>
#include <vector>

namespace qmcplusplus
{

TEMPLATE_TEST_CASE("RandomRotationMatrix", "[numerics]", float, double)
{
  // TestType is defined by Catch
  using TensorType = Tensor<TestType, 3>;

  TensorType rmat = generateRotationMatrix<TestType>(0.0, 0.0, 0.0);

  // Default rotation matrix should be the identity
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      if (i == j)
        CHECK(rmat(i, j) == Approx(1.0));
      else
        CHECK(rmat(i, j) == Approx(0.0));


  TensorType rmat2 = generateRotationMatrix<TestType>(0.1, 0.2, 0.3);

  CHECK(rmat2(0, 0) == Approx(-0.459016994374947));
  CHECK(rmat2(0, 1) == Approx(0.842075137094350));
  CHECK(rmat2(0, 2) == Approx(-0.283218753550188));
  CHECK(rmat2(1, 0) == Approx(-0.489403985718866));
  CHECK(rmat2(1, 1) == Approx(0.0263932022500210));
  CHECK(rmat2(1, 2) == Approx(0.871657695220709));
  CHECK(rmat2(2, 0) == Approx(0.741476323045772));
  CHECK(rmat2(2, 1) == Approx(0.538714082201795));
  CHECK(rmat2(2, 2) == Approx(0.400000000000000));


  TensorType rmat3 = generateRotationMatrix<TestType>(0.9, 0.5, 0.8);

  CHECK(rmat3(0, 0) == Approx(0.485410196624969));
  CHECK(rmat3(0, 1) == Approx(-0.352671151375484));
  CHECK(rmat3(0, 2) == Approx(0.800000000000000));
  CHECK(rmat3(1, 0) == Approx(-0.587785252292473));
  CHECK(rmat3(1, 1) == Approx(-0.809016994374947));
  CHECK(rmat3(1, 2) == Approx(9.79717439317882e-17));
  CHECK(rmat3(2, 0) == Approx(0.647213595499958));
  CHECK(rmat3(2, 1) == Approx(-0.470228201833979));
  CHECK(rmat3(2, 2) == Approx(-0.600000000000000));
}

} // namespace qmcplusplus

//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2023 QMCPACK developers
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
#include <catch2/catch_test_macros.hpp>
#include "type_traits/complex_help.hpp"

namespace qmcplusplus
{
template<typename P>
class TestComplexHelper
{
  using Cmplx        = std::complex<P>;
  using Real         = RealAlias<Cmplx>;
  using CmplxRebuild = ValueAlias<P, Cmplx>;
  using RealRebuild  = ValueAlias<P, Real>;

public:
  void run()
  {
    Cmplx aa;
    CmplxRebuild bb;
    aa = bb;
    static_cast<void>(aa);
    static_cast<void>(bb);

    Real cc;
    RealRebuild dd(0);
    cc = dd;
    static_cast<void>(cc);
    static_cast<void>(dd);
  }
};

TEST_CASE("complex_helper", "[type_traits]")
{
  TestComplexHelper<float> float_test;
  float_test.run();
  TestComplexHelper<double> double_test;
  double_test.run();
}

} // namespace qmcplusplus

//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "OhmmsPETE/OhmmsMatrix.h"
#include "QMCWaveFunctions/Fermion/MultiDiracDeterminant.h"

//#include <stdio.h>
#include <string>

using std::string;


namespace qmcplusplus
{
template<class T>
class TestMultiDiracDeterminantCalculator
{
public:
  template<typename DT>
  using OffloadMatrix = Matrix<DT, OffloadPinnedAllocator<DT>>;
  T default_evaluate(int power_of_two)
  {
    MultiDiracDeterminantCalculator<double> MDDC;
    int power2 = std::pow(2, power_of_two);
    MDDC.resize(power2);
    OffloadMatrix<double> dots(2 * power2); //This is an 2n by 2n matrix if you don't reduce pairs
    double n = 0.0;
    int i    = 0;
    //Just making some non trivial data
    for (auto& m : dots)
    {
      if (++i % 2 != 0)
        m = -n / (T)power2;
      else
        m = n / (T)power2;
      if (n > (T)power2 - 0.5)
        n = 0.0;
      else
        n += 1.0;
    }
    std::vector<int> it_things(power2 * power2);
    i = 0;
    for (auto& itt : it_things)
      itt = i++;
    std::vector<int>::const_iterator it = it_things.begin();
    return MDDC.evaluate(dots, it, power2);
  }
};

/** Simple synthetic test case will trip on changes in this method.
 */
TEST_CASE("MultiDiracDeterminantCalculator::evaluate-Small", "[wavefunction][fermion][multidet]")
{
  TestMultiDiracDeterminantCalculator<double> double_test;
  double det_value_expect = -1.1086723208;
  REQUIRE(double_test.default_evaluate(3) == Approx(det_value_expect));
  det_value_expect = -1.3432116824;
  REQUIRE(double_test.default_evaluate(7) == Approx(det_value_expect));
  det_value_expect = -1.3586431786;
  REQUIRE(double_test.default_evaluate(12) == Approx(det_value_expect));
}

} // namespace qmcplusplus

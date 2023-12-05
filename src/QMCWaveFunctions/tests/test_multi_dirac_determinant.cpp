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
class TestSmallMatrixDetCalculator
{
  template<typename DT>
  using OffloadMatrix = Matrix<DT, OffloadPinnedAllocator<DT>>;

  SmallMatrixDetCalculator<double> MDDC;
  OffloadMatrix<double> dots;
  std::vector<int> it_things;

public:
  void build_interal_data(int dim_size)
  {
    MDDC.resize(dim_size);
    dots.resize(dim_size, dim_size);
    it_things.resize(2 * dim_size);

    double n = 0.0;
    int i    = 0;
    //Just making some non trivial data
    for (auto& m : dots)
    {
      if (++i % 2 != 0)
        m = -n / (T)dim_size;
      else
        m = n / (T)dim_size;
      if (n > (T)dim_size - 0.5)
        n = 0.0;
      else
        n += 1.0;
    }

    for (size_t i = 0; i < dim_size; i++)
      it_things[i] = it_things[i + dim_size] = i;
  }

  T generic_evaluate(int dim_size)
  {
    build_interal_data(dim_size);
    return MDDC.evaluate(dots, it_things.data(), dim_size);
  }

  template<unsigned EXT_LEVEL>
  T customized_evaluate()
  {
    build_interal_data(EXT_LEVEL);
    return calcSmallDeterminant(EXT_LEVEL, dots.data(), it_things.data(),dots.cols());
  }
};

/** Simple synthetic test case will trip on changes in this method.
 */
TEST_CASE("SmallMatrixDetCalculator::evaluate-Small", "[wavefunction][fermion][multidet]")
{
  TestSmallMatrixDetCalculator<double> double_test;
  CHECK(double_test.generic_evaluate(1) == Approx(0.0));
  CHECK(double_test.generic_evaluate(2) == Approx(0.5));
  CHECK(double_test.generic_evaluate(3) == Approx(-0.7407407407));
  CHECK(double_test.generic_evaluate(4) == Approx(-0.87890625));
  CHECK(double_test.generic_evaluate(5) == Approx(-0.96768));

  CHECK(double_test.customized_evaluate<1>() == Approx(0.0));
  CHECK(double_test.customized_evaluate<2>() == Approx(0.5));
  CHECK(double_test.customized_evaluate<3>() == Approx(-0.7407407407));
  CHECK(double_test.customized_evaluate<4>() == Approx(-0.87890625));
  CHECK(double_test.customized_evaluate<5>() == Approx(-0.96768));

  CHECK(double_test.generic_evaluate(1<<3) == Approx(-1.1086723208));
  CHECK(double_test.generic_evaluate(1<<7) == Approx(-1.3432116824));
  CHECK(double_test.generic_evaluate(1<<12) == Approx(-1.3586431786));
}

} // namespace qmcplusplus

//////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License.  See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:
//
// File created by:
// Miguel A. Morales, moralessilva2@llnl.gov
//    Lawrence Livermore National Laboratory
////////////////////////////////////////////////////////////////////////////////

#include "catch.hpp"
#include "Configuration.h"

#undef APP_ABORT
#define APP_ABORT(x) \
  {                  \
    std::cout << x;  \
    throw;           \
  }

#include <vector>
#include <iostream>


#include "AFQMC/Matrix/tests/matrix_helpers.h"
#include "AFQMC/Numerics/ma_blas.hpp"

#include "multi/array.hpp"
#include "multi/array_ref.hpp"

using boost::multi::array;
using boost::multi::array_ref;
using std::vector;
template<std::ptrdiff_t D>
using iextensions = typename boost::multi::iextensions<D>;

namespace qmcplusplus
{
void ma_tensor_tests()
{
  vector<double> v = {1., 2., 3.};
  {
    array_ref<double, 1> V(v.data(), iextensions<1u>{v.size()});
    ma::scal(2., V);
    {
      vector<double> v2 = {2., 4., 6.};
      array_ref<double, 1> V2(v2.data(), iextensions<1u>{v2.size()});
      verify_approx(V, V2);
    }
  }
}

TEST_CASE("test_tensor", "[tensor_operations]") { ma_tensor_tests(); }

} // namespace qmcplusplus

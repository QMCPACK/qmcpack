//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2018 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//
// File created by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "OhmmsSoA/VectorSoaContainer.h"
#include "spline2/MultiBsplineMPIShared.hpp"
#include "spline2/SingleBsplineAllocator.hpp"
#include "spline2/MultiBsplineEval.hpp"
#include "QMCWaveFunctions/BsplineFactory/contraction_helper.hpp"
#include "config/stdlib/Constants.h"

namespace qmcplusplus
{

/** Supports testing many sizes of splines for benchmarking
 *  modified from einspline/tests/test_3d.cpp
 */

template<typename T, int GRID_SIZE>
class test_splines_base
{
protected:
  BCtype_d bc[3];
  Ugrid grid[3];

  const int N = GRID_SIZE;
  double delta;
  std::vector<double> data;

public:
  test_splines_base()
  {
    data.resize(N * N * N);

    grid[0].start = 0.0;
    grid[0].end   = 1.0;
    grid[0].num   = N;

    grid[1].start = 0.0;
    grid[1].end   = 1.0;
    grid[1].num   = N;

    grid[2].start = 0.0;
    grid[2].end   = 1.0;
    grid[2].num   = N;

    delta = (grid[0].end - grid[0].start) / grid[0].num;

    bc[0].lCode = PERIODIC;
    bc[0].rCode = PERIODIC;
    bc[0].lVal  = 0.0;
    bc[0].rVal  = 0.0;
    bc[1].lCode = PERIODIC;
    bc[1].rCode = PERIODIC;
    bc[1].lVal  = 0.0;
    bc[1].rVal  = 0.0;
    bc[2].lCode = PERIODIC;
    bc[2].rCode = PERIODIC;
    bc[2].lVal  = 0.0;
    bc[2].rVal  = 0.0;

    double tpi = 2 * M_PI;
    // Generate the data in double precision regardless of the target spline precision
    for (int i = 0; i < N; i++)
      for (int j = 0; j < N; j++)
        for (int k = 0; k < N; k++)
        {
          double x                    = delta * i;
          double y                    = delta * j;
          double z                    = delta * k;
          data[i * N * N + j * N + k] = std::sin(tpi * x) + std::sin(3 * tpi * y) + std::sin(4 * tpi * z);
        }
  }
};

/** Unspecialized test_splines doesn't test eval values
 */
template<typename T, int GRID_SIZE = 5>
struct test_splines : public test_splines_base<T, GRID_SIZE>
{
  using base = test_splines_base<T, GRID_SIZE>;
  using base::bc;
  using base::data;
  using base::delta;
  using base::grid;
  using base::N;

  void test(size_t num_splines)
  {
    auto comm_distributed = std::make_unique<Communicate>(*OHMMS::Controller, OHMMS::Controller->size());
    auto& comm(*comm_distributed);
    MultiBsplineMPIShared<T> bs(grid, bc, num_splines, std::move(comm_distributed));

    const size_t npad = getAlignedSize<T>(num_splines);
    REQUIRE(bs.num_splines_padded() == getAlignedSize<T>(num_splines));

    SingleBsplineAllocator<double> mAllocator;
    UBspline_3d_d* aspline = mAllocator.allocateUBspline(grid[0], grid[1], grid[2], bc[0], bc[1], bc[2], data.data());

    auto offsets = FairDivideAligned<std::vector<size_t>>(num_splines, getAlignment<T>(), comm.size());
    for (int i = offsets[comm.rank()]; i < offsets[comm.rank() + 1]; i++)
      bs.set_spline(*aspline, i);
    comm.barrier();
    mAllocator.destroy(aspline);

    //  The values for N=5 are not good enough for finer grids so by default we don't do those checks

    TinyVector<T, 3> pos = {0, 0, 0};

    aligned_vector<T> v(npad);
    bs.evaluate_v(pos, v);

    VectorSoaContainer<T, 3> dv(npad);
    VectorSoaContainer<T, 6> hess(npad);
    bs.evaluate_vgh(pos, v, dv, hess);

    pos = {0.1, 0.2, 0.3};
    bs.evaluate_v(pos, v);

    bs.evaluate_vgh(pos, v, dv, hess);

    VectorSoaContainer<T, 3> lap(npad);
    bs.evaluate_vgl(pos, v, dv, lap);

    VectorSoaContainer<T, 10> ghess(npad);
    bs.evaluate_vghgh(pos, v, dv, hess, ghess);
  }
};

/** Partially specialized test_splines does test eval values
 * See gen_bspline_values.py
 */
template<typename T>
struct test_splines<T, 5> : public test_splines_base<T, 5>
{
  using base = test_splines_base<T, 5>;
  using base::bc;
  using base::data;
  using base::delta;
  using base::grid;
  using base::N;

  void test(size_t num_splines)
  {
    auto comm_distributed = std::make_unique<Communicate>(*OHMMS::Controller, OHMMS::Controller->size());
    auto& comm(*comm_distributed);
    MultiBsplineMPIShared<T> bs(grid, bc, num_splines, std::move(comm_distributed));

    const size_t npad = getAlignedSize<T>(num_splines);
    REQUIRE(bs.num_splines_padded() == getAlignedSize<T>(num_splines));

    SingleBsplineAllocator<double> mAllocator;
    UBspline_3d_d* aspline = mAllocator.allocateUBspline(grid[0], grid[1], grid[2], bc[0], bc[1], bc[2], data.data());

    auto offsets = FairDivideAligned<std::vector<size_t>>(num_splines, getAlignment<T>(), comm.size());
    for (int i = offsets[comm.rank()]; i < offsets[comm.rank() + 1]; i++)
      bs.set_spline(*aspline, i);
    comm.barrier();
    mAllocator.destroy(aspline);

    //  Code from here to the end of the function is generated by gen_bspline_values.py

    TinyVector<T, 3> pos = {0, 0, 0};

    // symbolic value at pos =  (cx[0]/6 + 2*cx[1]/3 + cx[2]/6)*(cy[0]/6 + 2*cy[1]/3 + cy[2]/6)*(cz[0]/6 + 2*cz[1]/3 + cz[2]/6)
    aligned_vector<T> v(npad);
    bs.evaluate_v(pos, v);
    CHECK(v[0] == Approx(-3.529930688e-12));
    return;

    VectorSoaContainer<T, 3> dv(npad);
    VectorSoaContainer<T, 6> hess(npad);
    bs.evaluate_vgh(pos, v, dv, hess);
    // Gradient
    CHECK(dv[0][0] == Approx(6.178320809));
    CHECK(dv[0][1] == Approx(-7.402942564));
    CHECK(dv[0][2] == Approx(-6.178320809));

    // Hessian
    for (int i = 0; i < 6; i++)
    {
      CHECK(hess[0][i] == Approx(0.0));
    }

    pos = {0.1, 0.2, 0.3};
    bs.evaluate_v(pos, v);

    // Value
    CHECK(v[0] == Approx(-0.9476393279));

    bs.evaluate_vgh(pos, v, dv, hess);
    // Value
    CHECK(v[0] == Approx(-0.9476393279));
    // Gradient
    CHECK(dv[0][0] == Approx(5.111042137));
    CHECK(dv[0][1] == Approx(5.989106342));
    CHECK(dv[0][2] == Approx(1.952244379));
    // Hessian
    CHECK(hess[0][0] == Approx(-21.34557341));
    CHECK(hess[0][1] == Approx(1.174505743e-09));
    CHECK(hess[0][2] == Approx(-1.1483271e-09));
    CHECK(hess[0][3] == Approx(133.9204891));
    CHECK(hess[0][4] == Approx(-2.15319293e-09));
    CHECK(hess[0][5] == Approx(34.53786329));


    VectorSoaContainer<T, 3> lap(npad);
    bs.evaluate_vgl(pos, v, dv, lap);
    // Value
    CHECK(v[0] == Approx(-0.9476393279));
    // Gradient
    CHECK(dv[0][0] == Approx(5.111042137));
    CHECK(dv[0][1] == Approx(5.989106342));
    CHECK(dv[0][2] == Approx(1.952244379));
    // Laplacian
    CHECK(lap[0][0] == Approx(147.1127789));

    VectorSoaContainer<T, 10> ghess(npad);
    bs.evaluate_vghgh(pos, v, dv, hess, ghess);
    // Value
    CHECK(v[0] == Approx(-0.9476393279));
    // Gradient
    CHECK(dv[0][0] == Approx(5.111042137));
    CHECK(dv[0][1] == Approx(5.989106342));
    CHECK(dv[0][2] == Approx(1.952244379));
    // Hessian
    CHECK(hess[0][0] == Approx(-21.34557341));
    CHECK(hess[0][1] == Approx(1.174505743e-09));
    CHECK(hess[0][2] == Approx(-1.1483271e-09));
    CHECK(hess[0][3] == Approx(133.9204891));

    CHECK(hess[0][4] == Approx(-2.15319293e-09));
    CHECK(hess[0][5] == Approx(34.53786329));


    // Catch default is 100*(float epsilson)
    double eps = 2000 * std::numeric_limits<float>::epsilon();

    // Gradient of Hessian
    CHECK(ghess[0][0] == Approx(-213.455734));
    CHECK(ghess[0][1] == Approx(2.311193459e-09).epsilon(eps));
    CHECK(ghess[0][2] == Approx(3.468205279e-09).epsilon(eps));
    CHECK(ghess[0][3] == Approx(1.58092329e-07).epsilon(eps));
    CHECK(ghess[0][4] == Approx(1.255694171e-08).epsilon(eps));
    CHECK(ghess[0][5] == Approx(4.78981157e-08).epsilon(eps));
    CHECK(ghess[0][6] == Approx(-1753.041961));
    CHECK(ghess[0][7] == Approx(-2.575826885e-09).epsilon(eps));
    CHECK(ghess[0][8] == Approx(-4.683496702e-09).epsilon(eps));
    CHECK(ghess[0][9] == Approx(-81.53283531));
  }
};

TEST_CASE("MultiBsplineMPIShared periodic double", "[spline2]") { test_splines<double>().test(13); }

//TEST_CASE("MultiBsplineMPIShared periodic float", "[spline2]") { test_splines<float>().test(11); }

} // namespace qmcplusplus

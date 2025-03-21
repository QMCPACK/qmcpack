/////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2025 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//////////////////////////////////////////////////////////////////////////////////////

#include "SpaceGridTest.hpp"
#include "GenerateRandomParticleSets.h"
#include "catch.hpp"

#include "NESpaceGrid.h"
#include "SpaceGridInput.h"
#include "NEReferencePoints.h"
#include "ReferencePointsInput.h"
#include "ValidReferencePointsInput.h"
#include "ValidSpaceGridInput.h"
#include "OhmmsData/Libxml2Doc.h"
#include "EstimatorTesting.h"
#include "Particle/tests/MinimalParticlePool.h"
#include "NativeInitializerPrint.hpp"
#include <checkVector.hpp>
#include <memory>
#include <numeric>
/** \file
 *  This is a postfacto unit testing written for NESpaceGrid during porting of EnergyDensity
 *  to the batched version of the estimators.
 */
namespace qmcplusplus
{

constexpr bool generate_test_data = false;
using Real                        = double;

const ParticleSet::ParticlePos default_start_pos{{1.451870349, 1.381521229, 1.165202269},
                                                 {1.244515371, 1.382273176, 1.21105285},
                                                 {0.000459944, 1.329603408, 1.265030556},
                                                 {0.748660329, 1.63420622, 1.393637791},
                                                 {0.033228526, 1.391869137, 0.654413566},
                                                 {1.114198787, 1.654334594, 0.231075822},
                                                 {1.657151589, 0.883870516, 1.201243939},
                                                 {0.97317591, 1.245644974, 0.284564732}};

template<typename REAL>
REAL tensorAccessor(const std::vector<REAL>& grid_data, int i, int j, int k, int iv)
{
  return grid_data[1200 * i + 60 * j + 3 * k + iv];
}

template<typename REAL>
REAL tensorAccessor(const Matrix<REAL>& grid_data, int i, int j, int k, int iv)
{
  return grid_data.data()[1200 * i + 60 * j + 3 * k + iv];
}

template<typename REAL>
REAL tensorAccessor(const NESpaceGrid<REAL>& grid, int i, int j, int k, int iv)
{
  auto dm   = grid.getDM();
  auto nvpd = grid.getNValuesPerDomain();
  return grid.getDataVector()[dm[0] * nvpd * i + dm[1] * nvpd * j + nvpd * k + iv];
}

template<typename REAL>
REAL tensorAccessor(const NESpaceGrid<REAL>& grid, std::array<int, 3> gind, int iv)
{
  auto dm   = grid.getDM();
  auto nvpd = grid.getNValuesPerDomain();
  return grid.getDataVector()[dm[0] * nvpd * gind[0] + dm[1] * nvpd * gind[1] + nvpd * gind[2] + iv];
}

TEST_CASE("SpaceGrid::Construction", "[estimators]")
{
  using Input = testing::ValidSpaceGridInput;
  Communicate* comm;
  comm = OHMMS::Controller;

  testing::SpaceGridEnv<Input::valid::ORIGIN> sge(comm);

  // EnergyDensityEstimator gets this from an enum giving indexes into each offset of SOA buffer.
  // It is a smell.
  NESpaceGrid<Real> space_grid(*(sge.sgi_), sge.ref_points_->get_points(), 1, false);

  using NES         = testing::NESpaceGridTests<double>;
  auto buffer_start = NES::getBufferStart(space_grid);
  auto buffer_end   = NES::getBufferEnd(space_grid);
  space_grid.write_description(std::cout, std::string(""));
  auto& sgi = *(sge.sgi_);
  auto& agr = sgi.get_axis_grids();
  for (int id = 0; id < OHMMS_DIM; ++id)
    CHECK(NES::getOdu(space_grid)[id] == agr[id].odu);
}

TEST_CASE("SpaceGrid::CYLINDRICAL", "[estimators]")
{
  using Input = testing::ValidSpaceGridInput;
  Communicate* comm;
  comm = OHMMS::Controller;

  testing::SpaceGridEnv<Input::valid::CYLINDRICAL> sge(comm);

  // EnergyDensityEstimator gets this from an enum giving indexes into each offset of SOA buffer.
  // It is a smell.
  NESpaceGrid<Real> space_grid(*(sge.sgi_), sge.ref_points_->get_points(), 1, false);

  using NES         = testing::NESpaceGridTests<double>;
  auto buffer_start = NES::getBufferStart(space_grid);
  auto buffer_end   = NES::getBufferEnd(space_grid);
  space_grid.write_description(std::cout, std::string(""));
  auto& sgi = *(sge.sgi_);
  auto& agr = sgi.get_axis_grids();
  for (int id = 0; id < OHMMS_DIM; ++id)
    CHECK(NES::getOdu(space_grid)[id] == agr[id].odu);
}

TEST_CASE("SpaceGrid::SPHERICAL", "[estimators]")
{
  using Input = testing::ValidSpaceGridInput;
  Communicate* comm;
  comm = OHMMS::Controller;

  testing::SpaceGridEnv<Input::valid::SPHERICAL> sge(comm);

  // EnergyDensityEstimator gets this from an enum giving indexes into each offset of SOA buffer.
  // It is a smell.
  NESpaceGrid<Real> space_grid(*(sge.sgi_), sge.ref_points_->get_points(), 1, false);

  using NES         = testing::NESpaceGridTests<double>;
  auto buffer_start = NES::getBufferStart(space_grid);
  auto buffer_end   = NES::getBufferEnd(space_grid);
  space_grid.write_description(std::cout, std::string(""));
  auto& sgi = *(sge.sgi_);
  auto& agr = sgi.get_axis_grids();
  for (int id = 0; id < OHMMS_DIM; ++id)
    CHECK(NES::getOdu(space_grid)[id] == agr[id].odu);
}

TEST_CASE("SpaceGrid::Basic", "[estimators]")
{
  using Input = testing::ValidSpaceGridInput;
  Communicate* comm;
  comm = OHMMS::Controller;
  testing::SpaceGridEnv<Input::valid::ORIGIN> sge(comm);
  sge.pset_elec_.R = default_start_pos;
  int num_values   = 3;
  NESpaceGrid<Real> space_grid(*(sge.sgi_), sge.ref_points_->get_points(), num_values, true);
  using NES         = testing::NESpaceGridTests<double>;
  auto buffer_start = NES::getBufferStart(space_grid);
  auto buffer_end   = NES::getBufferEnd(space_grid);
  space_grid.write_description(std::cout, std::string(""));
  auto& sgi = *(sge.sgi_);
  auto& agr = sgi.get_axis_grids();
  for (int id = 0; id < OHMMS_DIM; ++id)
    CHECK(NES::getOdu(space_grid)[id] == agr[id].odu);
  //CHECK(buffer_start == 0);
  //CHECK(buffer_end == 23999);

  CHECK(space_grid.nDomains() == 8000);

  CHECK(space_grid.getDataVector().size() == 24000);

  Matrix<Real> values;
  values.resize(sge.pset_elec_.getTotalNum(), num_values);

  for (int ip = 0; ip < sge.pset_elec_.getTotalNum(); ++ip)
    for (int iv = 0; iv < num_values; ++iv)
      values(ip, iv) = ip + 0.1 * iv;

  const int ei_tid = sge.pset_elec_.addTable(sge.pset_ions_);
  sge.pset_elec_.update();
  sge.pset_ions_.update();

  std::vector<bool> p_outside(8, false);
  space_grid.accumulate(sge.pset_elec_.R, values, p_outside, sge.pset_elec_.getDistTableAB(ei_tid));

  // check that what's in data_pool is what is expected.
  const auto& grid_data = NES::getData(space_grid);

  CHECK(tensorAccessor(grid_data, 10, 17, 9, 0) == Approx(2.0));
  CHECK(tensorAccessor(grid_data, 10, 17, 9, 1) == Approx(2.1));
  CHECK(tensorAccessor(grid_data, 10, 17, 9, 2) == Approx(2.2));
  CHECK(tensorAccessor(grid_data, 12, 15, 7, 0) == Approx(4.0));
  CHECK(tensorAccessor(grid_data, 12, 15, 7, 1) == Approx(4.1));
  CHECK(tensorAccessor(grid_data, 12, 15, 7, 2) == Approx(4.2));
  CHECK(tensorAccessor(grid_data, 12, 16, 11, 0) == Approx(3.0));
  CHECK(tensorAccessor(grid_data, 12, 16, 11, 1) == Approx(3.1));
  CHECK(tensorAccessor(grid_data, 12, 16, 11, 2) == Approx(3.2));
  CHECK(tensorAccessor(grid_data, 13, 11, 15, 0) == Approx(6.0));
  CHECK(tensorAccessor(grid_data, 13, 11, 15, 1) == Approx(6.1));
  CHECK(tensorAccessor(grid_data, 13, 11, 15, 2) == Approx(6.2));
  CHECK(tensorAccessor(grid_data, 14, 13, 13, 0) == Approx(1.0));
  CHECK(tensorAccessor(grid_data, 14, 13, 13, 1) == Approx(1.2));
  CHECK(tensorAccessor(grid_data, 14, 13, 13, 2) == Approx(1.4));
  CHECK(tensorAccessor(grid_data, 15, 11, 10, 0) == Approx(7.0));
  CHECK(tensorAccessor(grid_data, 15, 11, 10, 1) == Approx(7.1));
  CHECK(tensorAccessor(grid_data, 15, 11, 10, 2) == Approx(7.2));
  CHECK(tensorAccessor(grid_data, 17, 12, 9, 0) == Approx(5.0));
  CHECK(tensorAccessor(grid_data, 17, 12, 9, 1) == Approx(5.1));
  CHECK(tensorAccessor(grid_data, 17, 12, 9, 2) == Approx(5.2));
  // new pset R's
  // check again
  auto min_R =
      ParticleSet::ParticlePos{{-0.6759092808, 0.835668385, 1.985307097},   {0.09710352868, -0.76751858, -1.89306891},
                               {-0.5605484247, -0.9578875303, 1.476860642}, {2.585144997, 1.862680197, 3.282609463},
                               {-0.1961335093, 1.111888766, -0.578481257},  {1.794641614, 1.6000278, -0.9474347234},
                               {2.157717228, 0.9254754186, 2.263158321},    {1.883366346, 2.136350632, 3.188981533}};
  sge.pset_elec_.applyMinimumImage(min_R);
  sge.pset_elec_.R = min_R;

  sge.pset_elec_.update();

  std::vector<bool> p_outside_2(8, false);
  space_grid.accumulate(sge.pset_elec_.R, values, p_outside_2, sge.pset_elec_.getDistTableAB(ei_tid));

  CHECK(tensorAccessor(grid_data, 1, 13, 15, 0) == Approx(2.0));
  CHECK(tensorAccessor(grid_data, 1, 13, 15, 1) == Approx(2.1));
  CHECK(tensorAccessor(grid_data, 1, 13, 15, 2) == Approx(2.2));
  CHECK(tensorAccessor(grid_data, 2, 6, 7, 0) == Approx(5.0));
  CHECK(tensorAccessor(grid_data, 2, 6, 7, 1) == Approx(5.1));
  CHECK(tensorAccessor(grid_data, 2, 6, 7, 2) == Approx(5.2));
  CHECK(tensorAccessor(grid_data, 4, 0, 11, 1) == Approx(0.1));
  CHECK(tensorAccessor(grid_data, 4, 0, 11, 2) == Approx(0.2));
  CHECK(tensorAccessor(grid_data, 10, 17, 9, 0) == Approx(2.0));
  CHECK(tensorAccessor(grid_data, 10, 17, 9, 1) == Approx(2.1));
  CHECK(tensorAccessor(grid_data, 10, 17, 9, 2) == Approx(2.2));
  CHECK(tensorAccessor(grid_data, 12, 0, 18, 0) == Approx(7.0));
  CHECK(tensorAccessor(grid_data, 12, 0, 18, 1) == Approx(7.1));
  CHECK(tensorAccessor(grid_data, 12, 0, 18, 2) == Approx(7.2));
  CHECK(tensorAccessor(grid_data, 12, 13, 0, 0) == Approx(6.0));
  CHECK(tensorAccessor(grid_data, 12, 13, 0, 1) == Approx(6.1));
  CHECK(tensorAccessor(grid_data, 12, 13, 0, 2) == Approx(6.2));
  CHECK(tensorAccessor(grid_data, 12, 15, 7, 0) == Approx(4.0));
  CHECK(tensorAccessor(grid_data, 12, 15, 7, 1) == Approx(4.1));
  CHECK(tensorAccessor(grid_data, 12, 15, 7, 2) == Approx(4.2));
  CHECK(tensorAccessor(grid_data, 12, 16, 11, 0) == Approx(3.0));
  CHECK(tensorAccessor(grid_data, 12, 16, 11, 1) == Approx(3.1));
  CHECK(tensorAccessor(grid_data, 12, 16, 11, 2) == Approx(3.2));
  CHECK(tensorAccessor(grid_data, 13, 1, 6, 0) == Approx(1.0));
  CHECK(tensorAccessor(grid_data, 13, 1, 6, 1) == Approx(1.1));
  CHECK(tensorAccessor(grid_data, 13, 1, 6, 2) == Approx(1.2));
  CHECK(tensorAccessor(grid_data, 13, 11, 15, 0) == Approx(6.0));
  CHECK(tensorAccessor(grid_data, 13, 11, 15, 1) == Approx(6.1));
  CHECK(tensorAccessor(grid_data, 13, 11, 15, 2) == Approx(6.2));
  CHECK(tensorAccessor(grid_data, 13, 17, 1, 0) == Approx(3.0));
  CHECK(tensorAccessor(grid_data, 13, 17, 1, 1) == Approx(3.1));
  CHECK(tensorAccessor(grid_data, 13, 17, 1, 2) == Approx(3.2));
  CHECK(tensorAccessor(grid_data, 14, 12, 4, 0) == Approx(4.0));
  CHECK(tensorAccessor(grid_data, 14, 12, 4, 1) == Approx(4.1));
  CHECK(tensorAccessor(grid_data, 14, 12, 4, 2) == Approx(4.2));
  CHECK(tensorAccessor(grid_data, 14, 13, 13, 0) == Approx(1.0));
  CHECK(tensorAccessor(grid_data, 14, 13, 13, 1) == Approx(1.2));
  CHECK(tensorAccessor(grid_data, 14, 13, 13, 2) == Approx(1.4));
  CHECK(tensorAccessor(grid_data, 15, 11, 10, 0) == Approx(7.0));
  CHECK(tensorAccessor(grid_data, 15, 11, 10, 1) == Approx(7.1));
  CHECK(tensorAccessor(grid_data, 15, 11, 10, 2) == Approx(7.2));
  CHECK(tensorAccessor(grid_data, 17, 12, 9, 0) == Approx(5.0));
  CHECK(tensorAccessor(grid_data, 17, 12, 9, 1) == Approx(5.1));
  CHECK(tensorAccessor(grid_data, 17, 12, 9, 2) == Approx(5.2));
}

TEST_CASE("SpaceGrid::Accumulate::outside", "[estimators]")
{
  using Input = testing::ValidSpaceGridInput;
  Communicate* comm;
  comm = OHMMS::Controller;
  testing::SpaceGridEnv<Input::valid::CYLINDRICAL> sge(comm);
  int num_values = 3;
  NESpaceGrid<Real> space_grid(*(sge.sgi_), sge.ref_points_->get_points(), num_values, false);
  using NES         = testing::NESpaceGridTests<double>;
  auto buffer_start = NES::getBufferStart(space_grid);
  auto buffer_end   = NES::getBufferEnd(space_grid);
  space_grid.write_description(std::cout, std::string(""));
  auto& sgi = *(sge.sgi_);
  auto& agr = sgi.get_axis_grids();
  for (int id = 0; id < OHMMS_DIM; ++id)
    CHECK(NES::getOdu(space_grid)[id] == agr[id].odu);

  CHECK(space_grid.nDomains() == 2000);
  CHECK(space_grid.getDataVector().size() == 6000);

  Matrix<Real> values;
  values.resize(sge.pset_elec_.getTotalNum(), num_values);

  for (int ip = 0; ip < sge.pset_elec_.getTotalNum(); ++ip)
    for (int iv = 0; iv < num_values; ++iv)
      values(ip, iv) = ip + 0.1 * iv;

  const int ei_tid = sge.pset_elec_.addTable(sge.pset_ions_);
  sge.pset_elec_.update();
  sge.pset_ions_.update();

  std::vector<bool> p_outside(8, false);
  space_grid.accumulate(sge.pset_elec_.R, values, p_outside, sge.pset_elec_.getDistTableAB(ei_tid));

  // new pset R's
  // check again
  auto min_R =
      ParticleSet::ParticlePos{{1.883366346, 2.136350632, 3.188981533}, {0.09710352868, -0.76751858, -1.89306891}};
  sge.pset_elec_.applyMinimumImage(min_R);
  sge.pset_elec_.R = min_R;

  sge.pset_elec_.update();
  std::cout << NativePrint(p_outside) << '\n';

  std::vector<bool> p_outside_2(8, false);
  space_grid.accumulate(sge.pset_elec_.R, values, p_outside_2, sge.pset_elec_.getDistTableAB(ei_tid));

  std::cout << NativePrint(p_outside_2) << '\n';
}

TEST_CASE("SpaceGrid::BadPeriodic", "[estimators]")
{
  using Input = testing::ValidSpaceGridInput;
  Communicate* comm;
  comm = OHMMS::Controller;
  testing::SpaceGridEnv<Input::valid::ORIGIN> sge(comm);
  int num_values   = 3;
  sge.pset_elec_.R = default_start_pos;
  NESpaceGrid<Real> space_grid(*(sge.sgi_), sge.ref_points_->get_points(), num_values, false);

  using NES         = testing::NESpaceGridTests<double>;
  auto buffer_start = NES::getBufferStart(space_grid);
  auto buffer_end   = NES::getBufferEnd(space_grid);
  space_grid.write_description(std::cout, std::string(""));
  auto& sgi = *(sge.sgi_);
  auto& agr = sgi.get_axis_grids();
  for (int id = 0; id < OHMMS_DIM; ++id)
    CHECK(NES::getOdu(space_grid)[id] == agr[id].odu);
  // CHECK(buffer_start == Approx(0));
  // CHECK(buffer_end == Approx(23999));

  Matrix<Real> values;
  values.resize(sge.pset_elec_.getTotalNum(), num_values);

  for (int ip = 0; ip < sge.pset_elec_.getTotalNum(); ++ip)
    for (int iv = 0; iv < num_values; ++iv)
      values(ip, iv) = ip + 0.1 * iv;

  const int ei_tid = sge.pset_elec_.addTable(sge.pset_ions_);
  sge.pset_elec_.update();
  sge.pset_ions_.update();

  std::vector<bool> p_outside(8, false);
  // set a position a small amount outside relative to the grid density.
  // because numerical inaccuracy in the dot product to go from cartesian to lattice coordinates we need to check that
  // particles that end up within on index of the of the grid bounds aren't an issue.
  // the accumulate in the periodic case has a strong assumption that minimum image has been applied to the incoming coords.

  sge.pset_elec_.R[2] = {2.585144997, 1.862680197, 2.9};
  space_grid.accumulate(sge.pset_elec_.R, values, p_outside, sge.pset_elec_.getDistTableAB(ei_tid));

  // set a position outside of the cell more than numerical error could ever get
  // u = {-0.256273, 1.11711, 0.94714}
  sge.pset_elec_.R[2] = {1.451870349, 3.481521229, 1.165202269};

  CHECK_THROWS_AS(space_grid.accumulate(sge.pset_elec_.R, values, p_outside, sge.pset_elec_.getDistTableAB(ei_tid)),
                  std::runtime_error);

  try
  {
    space_grid.accumulate(sge.pset_elec_.R, values, p_outside, sge.pset_elec_.getDistTableAB(ei_tid));
  }
  catch (const std::exception& exc)
  {
    std::cout << exc.what() << '\n';
  }


  // This is just barely out of the unit cell in a negative direction.
  // u = {-1.01, 0.5, 0.5}
  sge.pset_elec_.R[2] = {-0.860156, 1.68658, -0.860156};
  // So it should not throw
  space_grid.accumulate(sge.pset_elec_.R, values, p_outside, sge.pset_elec_.getDistTableAB(ei_tid));
  // But it does set a value
  const auto& grid_data = NES::getData(space_grid);
  CHECK(tensorAccessor(grid_data, 14, 14, 0, 1) == Approx(2.1));

  // This is just barely out of the unit cell in a negative direction.
  // u = {1.01, 0.5, 0.5}
  sge.pset_elec_.R[2] = {2.54674, 1.68658, 2.54674};
  // So it should not throw
  space_grid.accumulate(sge.pset_elec_.R, values, p_outside, sge.pset_elec_.getDistTableAB(ei_tid));
  // But it does set a value
  CHECK(tensorAccessor(grid_data, 14, 14, 19, 1) == Approx(2.1));
}

// This should preserve some of the rather strange (to me) behavior of cartesian grids.
TEST_CASE("SpaceGrid::WeirdCartesian", "[estimators]")
{
  using Input = testing::ValidSpaceGridInput;
  Communicate* comm;
  comm = OHMMS::Controller;
  testing::SpaceGridEnv<Input::valid::WEIRD_CARTESIAN> sge(comm);
  sge.pset_elec_.R = default_start_pos;

  int num_values = 3;
  NESpaceGrid<Real> space_grid(*(sge.sgi_), sge.ref_points_->get_points(), num_values, true);
  using NES         = testing::NESpaceGridTests<double>;
  auto buffer_start = NES::getBufferStart(space_grid);
  auto buffer_end   = NES::getBufferEnd(space_grid);
  space_grid.write_description(std::cout, std::string(""));
  auto& sgi = *(sge.sgi_);
  auto& agr = sgi.get_axis_grids();
  std::cout << "AxisGridDomains: ";
  for (int id = 0; id < OHMMS_DIM; ++id)
  {
    CHECK(NES::getOdu(space_grid)[id] == agr[id].odu);
    std::cout << NativePrint(agr[id].ndom_int) << ',';
  }
  std::cout << '\n';
  //CHECK(buffer_start == 0);
  //CHECK(buffer_end == 23999);

  CHECK(space_grid.nDomains() == 2000);

  CHECK(space_grid.getDataVector().size() == 6000);


  Matrix<Real> values;
  values.resize(sge.pset_elec_.getTotalNum(), num_values);

  auto dm              = space_grid.getDM();
  auto nvpd            = space_grid.getNValuesPerDomain();
  auto total_size_grid = space_grid.nDomains() * space_grid.getNValuesPerDomain();
  std::vector<Real> expected_grid_vector(total_size_grid);
  auto getLinearIndex = [dm, nvpd](auto gind, int iv) {
    return dm[0] * nvpd * gind[0] + dm[1] * nvpd * gind[1] + nvpd * gind[2] + iv;
  };

  for (int ip = 0; ip < sge.pset_elec_.getTotalNum(); ++ip)
    for (int iv = 0; iv < num_values; ++iv)
    {
      values(ip, iv) = ip + 0.1 * iv;
      auto gind      = space_grid.findGMapIndexes(sge.pset_elec_.R[ip]);
      expected_grid_vector[getLinearIndex(gind, iv)] += values(ip, iv);
    }

  const int ei_tid = sge.pset_elec_.addTable(sge.pset_ions_);
  sge.pset_elec_.update();
  sge.pset_ions_.update();

  auto min_R = sge.pset_elec_.R;
  sge.pset_elec_.applyMinimumImage(min_R);

  std::vector<bool> p_outside(8, false);
  space_grid.accumulate(min_R, values, p_outside, sge.pset_elec_.getDistTableAB(ei_tid));

  // check that what's in data_pool is what is expected.
  const auto& grid_data = NES::getData(space_grid);

  auto check = checkVector(expected_grid_vector, grid_data);
  CHECKED_ELSE(check.result) { FAIL(check.result_message); }

  // new pset R's
  // check again
  min_R =
      ParticleSet::ParticlePos{{-0.6759092808, 0.835668385, 1.985307097},   {0.09710352868, -0.76751858, -1.89306891},
                               {-0.5605484247, -0.9578875303, 1.476860642}, {2.585144997, 1.862680197, 3.282609463},
                               {-0.1961335093, 1.111888766, -0.578481257},  {1.794641614, 1.6000278, -0.9474347234},
                               {2.157717228, 0.9254754186, 2.263158321},    {1.883366346, 2.136350632, 3.188981533}};
  sge.pset_elec_.applyMinimumImage(min_R);
  sge.pset_elec_.R = min_R;

  sge.pset_elec_.update();

  std::vector<bool> p_outside_2(8, false);
  space_grid.accumulate(sge.pset_elec_.R, values, p_outside_2, sge.pset_elec_.getDistTableAB(ei_tid));

  for (int ip = 0; ip < sge.pset_elec_.getTotalNum(); ++ip)
    for (int iv = 0; iv < num_values; ++iv)
    {
      auto gind = space_grid.findGMapIndexes(sge.pset_elec_.R[ip]);
      expected_grid_vector[getLinearIndex(gind, iv)] += values(ip, iv);
    }

  check = checkVector(expected_grid_vector, grid_data);
  CHECKED_ELSE(check.result) { FAIL(check.result_message); }
}

TEST_CASE("SpaceGrid::hdf5", "[estimators]")
{
  using Input = testing::ValidSpaceGridInput;
  Communicate* comm;
  comm = OHMMS::Controller;
  testing::SpaceGridEnv<Input::valid::ORIGIN> sge(comm);
  sge.pset_elec_.R = default_start_pos;
  int num_values   = 3;
  NESpaceGrid<Real> space_grid(*(sge.sgi_), sge.ref_points_->get_points(), num_values, false);
  using NES         = testing::NESpaceGridTests<double>;
  auto buffer_start = NES::getBufferStart(space_grid);
  auto buffer_end   = NES::getBufferEnd(space_grid);
  space_grid.write_description(std::cout, std::string(""));
  auto& sgi = *(sge.sgi_);
  auto& agr = sgi.get_axis_grids();
  for (int id = 0; id < OHMMS_DIM; ++id)
    CHECK(NES::getOdu(space_grid)[id] == agr[id].odu);
  //CHECK(buffer_start == 0);
  //CHECK(buffer_end == 23999);

  Matrix<Real> values;
  values.resize(sge.pset_elec_.getTotalNum(), num_values);

  for (int ip = 0; ip < sge.pset_elec_.getTotalNum(); ++ip)
    for (int iv = 0; iv < num_values; ++iv)
      values(ip, iv) = ip + 0.1 * iv;

  const int ei_tid = sge.pset_elec_.addTable(sge.pset_ions_);
  sge.pset_elec_.update();
  sge.pset_ions_.update();

  std::vector<bool> p_outside(8, false);
  space_grid.accumulate(sge.pset_elec_.R, values, p_outside, sge.pset_elec_.getDistTableAB(ei_tid));

  hdf_archive hd;
  std::string test_file{"spacegrid_test.hdf"};
  bool okay = hd.create(test_file);
  REQUIRE(okay);

  std::vector<ObservableHelper> h5desc;
  space_grid.registerGrid(hd, 0);

  space_grid.write(hd);

  hd.close();

  hdf_archive hd_read;
  bool okay_read = hd.open(test_file);
  hd.push("spacegrid1");
  //hdf5 values always end up as doubles
  Matrix<double> read_values(1, 24000);
  hd.readEntry(read_values, "value");

  auto value = tensorAccessor(read_values, 10, 17, 9, 0);

  CHECK(tensorAccessor(read_values, 10, 17, 9, 0) == Approx(2.0));
  CHECK(tensorAccessor(read_values, 10, 17, 9, 1) == Approx(2.1));
  CHECK(tensorAccessor(read_values, 10, 17, 9, 2) == Approx(2.2));
  CHECK(tensorAccessor(read_values, 12, 15, 7, 0) == Approx(4.0));
  CHECK(tensorAccessor(read_values, 12, 15, 7, 1) == Approx(4.1));
  CHECK(tensorAccessor(read_values, 12, 15, 7, 2) == Approx(4.2));
  CHECK(tensorAccessor(read_values, 12, 16, 11, 0) == Approx(3.0));
  CHECK(tensorAccessor(read_values, 12, 16, 11, 1) == Approx(3.1));
  CHECK(tensorAccessor(read_values, 12, 16, 11, 2) == Approx(3.2));
  CHECK(tensorAccessor(read_values, 13, 11, 15, 0) == Approx(6.0));
  CHECK(tensorAccessor(read_values, 13, 11, 15, 1) == Approx(6.1));
  CHECK(tensorAccessor(read_values, 13, 11, 15, 2) == Approx(6.2));
  CHECK(tensorAccessor(read_values, 14, 13, 13, 0) == Approx(1.0));
  CHECK(tensorAccessor(read_values, 14, 13, 13, 1) == Approx(1.2));
  CHECK(tensorAccessor(read_values, 14, 13, 13, 2) == Approx(1.4));
  CHECK(tensorAccessor(read_values, 15, 11, 10, 0) == Approx(7.0));
  CHECK(tensorAccessor(read_values, 15, 11, 10, 1) == Approx(7.1));
  CHECK(tensorAccessor(read_values, 15, 11, 10, 2) == Approx(7.2));
  CHECK(tensorAccessor(read_values, 17, 12, 9, 0) == Approx(5.0));
  CHECK(tensorAccessor(read_values, 17, 12, 9, 1) == Approx(5.1));
  CHECK(tensorAccessor(read_values, 17, 12, 9, 2) == Approx(5.2));

  /// \todo add additional hdf5 output checks
}

TEST_CASE("SpaceGrid::collect", "[estimators]")
{
  using Input = testing::ValidSpaceGridInput;
  Communicate* comm;
  comm = OHMMS::Controller;
  testing::SpaceGridEnv<Input::valid::ORIGIN> sge(comm);
  int num_values = 3;
  NESpaceGrid<Real> collection_space_grid(*(sge.sgi_), sge.ref_points_->get_points(), num_values, false);
  UPtrVector<NESpaceGrid<Real>> accumulating_space_grids;

  // this could all be factored out into a test helper class or function if we add more tests over multiple grids
  // note its similar to EnergyDensityTest test helper class code.
  testing::SpaceGridTest<Real, testing::ValidSpaceGridInput::valid::ORIGIN> sgt(sge, 4, generate_test_data);

  int num_accumulators = 4;
  for (int ia = 0; ia < num_accumulators; ++ia)
  {
    accumulating_space_grids.emplace_back(std::make_unique<NESpaceGrid<Real>>(collection_space_grid));
  }

  auto psets_elecs = sgt.getPSetList();

  int ei_tid;
  using NES = testing::NESpaceGridTests<double>;

  for (int ia = 0; ia < num_accumulators; ++ia)
  {
    ParticleSet& pset_elec = psets_elecs[ia];
    ei_tid                 = pset_elec.addTable(sge.pset_ions_);
  }

  auto dm              = collection_space_grid.getDM();
  auto nvpd            = collection_space_grid.getNValuesPerDomain();
  auto total_size_grid = collection_space_grid.nDomains() * collection_space_grid.getNValuesPerDomain();
  std::vector<Real> expected_grid_vector(total_size_grid);
  auto getLinearIndex = [dm, nvpd](auto gind, int iv) {
    return dm[0] * nvpd * gind[0] + dm[1] * nvpd * gind[1] + nvpd * gind[2] + iv;
  };

  auto doAccumulate = [&](double nudge) {
    for (int ia = 0; ia < num_accumulators; ++ia)
    {
      //Make a step of fake data
      Matrix<Real> values;
      values.resize(sge.pset_elec_.getTotalNum(), num_values);

      ParticleSet& pset_elec = psets_elecs[ia];
      auto min_R             = pset_elec.R;
      pset_elec.applyMinimumImage(min_R);
      pset_elec.R = min_R;
      pset_elec.update();
      sge.pset_ions_.update();

      auto val_for_part = [nudge](int ip, int iv, int ia) { return ip + 0.01 * iv - ia * 0.05 + nudge; };
      for (int ip = 0; ip < sge.pset_elec_.getTotalNum(); ++ip)
        for (int iv = 0; iv < num_values; ++iv)
        {
          values(ip, iv) = val_for_part(ip, iv, ia);
          auto gind      = accumulating_space_grids[ia]->findGMapIndexes(pset_elec.R[ip]);
          expected_grid_vector[getLinearIndex(gind, iv)] += values(ip, iv);
        }

      std::vector<bool> p_outside(8, false);
      accumulating_space_grids[ia]->accumulate(pset_elec.R, values, p_outside, pset_elec.getDistTableAB(ei_tid));

      //spot check
      auto grid_indexes = accumulating_space_grids[ia]->findGMapIndexes(pset_elec.R[0]);
      std::cout << "grid_indexes: " << NativePrint(grid_indexes) << " value: " << val_for_part(0, 1, ia) << '\n';
      CHECK(tensorAccessor(*(accumulating_space_grids[ia]), grid_indexes, 2) == Approx(val_for_part(0, 2, ia)));
      grid_indexes = accumulating_space_grids[ia]->findGMapIndexes(pset_elec.R[3]);
      std::cout << "grid_indexes: " << NativePrint(grid_indexes) << " value: " << val_for_part(3, 1, ia) << '\n';
      CHECK(tensorAccessor(*(accumulating_space_grids[ia]), grid_indexes, 2) == Approx(val_for_part(3, 2, ia)));
    }
  };

  doAccumulate(0.0);
  // emphasize this is a static call.
  decltype(collection_space_grid)::collect(collection_space_grid,
                                           convertUPtrToRefVector(accumulating_space_grids));
  for (auto& accumulating_space_grid : accumulating_space_grids)
    accumulating_space_grid->zero();
  // here we check our expected grid that is trivially acumulated
  // not as far as the discretized mapping everything depends on NESpaceGride::finGMapIndexes
  auto check = checkVector(expected_grid_vector, collection_space_grid.getDataVector());
  CHECKED_ELSE(check.result) { FAIL(check.result_message); }

  doAccumulate(0.2);
  // emphasize this is a static call.
  decltype(collection_space_grid)::collect(collection_space_grid,
                                           convertUPtrToRefVector(accumulating_space_grids));
  for (auto& accumulating_space_grid : accumulating_space_grids)
    accumulating_space_grid->zero();

  // here we check our expected grid that is trivially acumulated
  // not as far as the discretized mapping everything depends on NESpaceGride::finGMapIndexes
  check = checkVector(expected_grid_vector, collection_space_grid.getDataVector());
  CHECKED_ELSE(check.result) { FAIL(check.result_message); }
}

} // namespace qmcplusplus

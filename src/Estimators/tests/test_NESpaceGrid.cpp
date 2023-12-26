//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2023 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//////////////////////////////////////////////////////////////////////////////////////

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
/** \file
 *  This is a postfacto unit testing written for NESpaceGrid during porting of EnergyDensity
 *  to the batched version of the estimators.
 */

namespace qmcplusplus
{

namespace testing
{
template<typename REAL>
class NESpaceGridTests
{
public:
  static const auto& getData(const NESpaceGrid<REAL>& nesg) { return nesg.data_; }
  static int getBufferStart(const NESpaceGrid<REAL>& nesg) { return nesg.buffer_start_; }
  static int getBufferEnd(const NESpaceGrid<REAL>& nesg) { return nesg.buffer_end_; }
  static auto* getOdu(const NESpaceGrid<REAL>& nesg) { return nesg.odu_; }
};

template<ValidSpaceGridInput::valid VALID>
class SpaceGridEnv
{
public:
  using Input = ValidSpaceGridInput;
  SpaceGridEnv(Communicate* comm)
      : particle_pool_(MinimalParticlePool::make_diamondC_1x1x1(comm)),
        pset_elec_(*(particle_pool_.getParticleSet("e"))),
        pset_ions_(*(particle_pool_.getParticleSet("ion")))
  {
    // Setup particleset
    // particle positions must be inside the unit cell
    pset_elec_.R =
        ParticleSet::ParticlePos{{1.451870349, 1.381521229, 1.165202269}, {1.244515371, 1.382273176, 1.21105285},
                                 {0.000459944, 1.329603408, 1.265030556}, {0.748660329, 1.63420622, 1.393637791},
                                 {0.033228526, 1.391869137, 0.654413566}, {1.114198787, 1.654334594, 0.231075822},
                                 {1.657151589, 0.883870516, 1.201243939}, {0.97317591, 1.245644974, 0.284564732}};

    Libxml2Document doc;
    bool okay       = doc.parseFromString(Input::xml[VALID]);
    xmlNodePtr node = doc.getRoot();
    sgi_            = std::make_unique<SpaceGridInput>(node);

    using RPInput = ValidReferencePointsInputs;
    Libxml2Document doc2;
    bool okay2       = doc.parseFromString(RPInput::xml[RPInput::CELL]);
    xmlNodePtr node2 = doc.getRoot();
    rpi_             = std::make_unique<ReferencePointsInput>(node2);
    ref_psets_.push_back(pset_ions_);
    ref_points_ = std::make_unique<NEReferencePoints>(*rpi_, pset_elec_, ref_psets_);
  }
  UPtr<ReferencePointsInput> rpi_;
  UPtr<SpaceGridInput> sgi_;
  UPtr<NEReferencePoints> ref_points_;
  RefVector<ParticleSet> ref_psets_;
  ParticleSetPool particle_pool_;
  ParticleSet pset_elec_;
  ParticleSet pset_ions_;
};

template<>
class SpaceGridEnv<ValidSpaceGridInput::valid::CYLINDRICAL>
{
public:
  using Input = ValidSpaceGridInput;
  SpaceGridEnv(Communicate* comm)
      : particle_pool_(MinimalParticlePool::make_H2(comm)),
        pset_elec_(*(particle_pool_.getParticleSet("e"))),
        pset_ions_(*(particle_pool_.getParticleSet("ion")))
  {
    Libxml2Document doc;
    bool okay       = doc.parseFromString(Input::xml[ValidSpaceGridInput::valid::CYLINDRICAL]);
    xmlNodePtr node = doc.getRoot();
    sgi_            = std::make_unique<SpaceGridInput>(node);

    using RPInput = ValidReferencePointsInputs;
    Libxml2Document doc2;
    bool okay2       = doc.parseFromString(RPInput::xml[RPInput::CELL]);
    xmlNodePtr node2 = doc.getRoot();
    rpi_             = std::make_unique<ReferencePointsInput>(node2);
    ref_psets_.push_back(pset_ions_);
    ref_points_ = std::make_unique<NEReferencePoints>(*rpi_, pset_elec_, ref_psets_);
  }
  UPtr<ReferencePointsInput> rpi_;
  UPtr<SpaceGridInput> sgi_;
  UPtr<NEReferencePoints> ref_points_;
  RefVector<ParticleSet> ref_psets_;
  ParticleSetPool particle_pool_;
  ParticleSet pset_elec_;
  ParticleSet pset_ions_;
};

} // namespace testing

constexpr bool generate_test_data = false;
using Real                        = double;

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

  // CHECK(buffer_start == 0);
  // CHECK(buffer_end == 7999);
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
  int num_values = 3;
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

  auto tensorAccessor = [](const auto& grid_data, int i, int j, int k, int iv) {
    return grid_data[1200 * i + 60 * j + 3 * k + iv];
  };

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
    ParticleSet::ParticlePos{{1.883366346, 2.136350632, 3.188981533},   {0.09710352868, -0.76751858, -1.89306891}};
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

  // set a position outside of the cell
  sge.pset_elec_.R[2] = {1.451870349, 4.381521229, 1.165202269};

  std::vector<bool> p_outside(8, false);

  CHECK_THROWS_AS(space_grid.accumulate(sge.pset_elec_.R, values, p_outside, sge.pset_elec_.getDistTableAB(ei_tid)),
                  std::runtime_error);
}

TEST_CASE("SpaceGrid::hdf5", "[estimators]")
{
  using Input = testing::ValidSpaceGridInput;
  Communicate* comm;
  comm = OHMMS::Controller;
  testing::SpaceGridEnv<Input::valid::ORIGIN> sge(comm);
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

  auto tensorAccessor = [](const auto& grid_data, int i, int j, int k, int iv) -> double {
    return grid_data.data()[1200 * i + 60 * j + 3 * k + iv];
  };

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
} // namespace qmcplusplus

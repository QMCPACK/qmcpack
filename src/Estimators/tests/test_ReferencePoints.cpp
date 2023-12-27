//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2023 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//////////////////////////////////////////////////////////////////////////////////////

#include "catch.hpp"

#include "NEReferencePoints.h"
#include "ReferencePointsInput.h"
#include "ValidReferencePointsInput.h"
#include "OhmmsData/Libxml2Doc.h"
#include "EstimatorTesting.h"
#include "Particle/tests/MinimalParticlePool.h"

/** \file
 *  This is a postfacto unit testing written for reference points during porting of EnergyDensity
 *  to the batched version of the estimators.
 */

namespace qmcplusplus
{

namespace testing
{
class TestableNEReferencePoints : public NEReferencePoints
{
public:
  TestableNEReferencePoints(const NEReferencePoints& nerp) : NEReferencePoints(nerp) {}
  void write_testable_description(std::ostream& os) const;
};

void TestableNEReferencePoints::write_testable_description(std::ostream& os) const
{
  os << "{" << '\n';
  std::map<std::string, Point>::const_iterator it, end = points_.end();
  for (it = points_.begin(); it != end; ++it)
  {
    os << " {\"" << it->first << "\", {" << std::setw(16) << std::setprecision(16) << it->second[0] << ","
       << it->second[1] << "," << it->second[2] << "}}," << '\n';
  }
  os << "};" << '\n';
  return;
}
} // namespace testing

std::ostream& operator<<(std::ostream& out, const testing::TestableNEReferencePoints& rhs)
{
  rhs.write_testable_description(out);
  return out;
}


std::ostream& operator<<(std::ostream& out, const testing::TestableNEReferencePoints& rhs);

constexpr bool generate_test_data = false;

template<typename T1, typename T2, unsigned D>
bool approxEquality(const TinyVector<T1, D>& val_a, const TinyVector<T2, D>& val_b)
{
  for (int i = 0; i < D; ++i)
    if (val_a[i] != Approx(val_b[i]))
      return false;
  return true;
}

ReferencePointsInput makeTestRPI()
{
  using Input = testing::ValidReferencePointsInputs;
  Libxml2Document doc;
  bool okay       = doc.parseFromString(Input::xml[Input::valid::CELL]);
  xmlNodePtr node = doc.getRoot();
  return {node};
}

struct PSetsAndRefList
{
  ParticleSetPool ppool;
  ParticleSet pset;
  ParticleSet pset_ions;
  RefVector<ParticleSet> ref_psets;
};

PSetsAndRefList makePsets()
{
  auto lattice = testing::makeTestLattice();
  Communicate* comm;
  comm               = OHMMS::Controller;
  auto particle_pool = MinimalParticlePool::make_diamondC_1x1x1(comm);
  auto& pset         = *(particle_pool.getParticleSet("e"));
  auto& pset_ions    = *(particle_pool.getParticleSet("ion"));

  // Setup particleset
  pset.R = ParticleSet::ParticlePos{{1.751870349, 4.381521229, 2.865202269}, {3.244515371, 4.382273176, 4.21105285},
                                    {3.000459944, 3.329603408, 4.265030556}, {3.748660329, 3.63420622, 5.393637791},
                                    {3.033228526, 3.391869137, 4.654413566}, {3.114198787, 2.654334594, 5.231075822},
                                    {3.657151589, 4.883870516, 4.201243939}, {2.97317591, 4.245644974, 4.284564732}};

  RefVector<ParticleSet> ref_psets;
  ref_psets.push_back(pset_ions);
  return {std::move(particle_pool), pset, pset_ions, ref_psets};
}

auto expectedReferencePoints()
{
  typename NEReferencePoints::Points expected_reference_points;
  if constexpr (std::is_same_v<NEReferencePoints::Real, double>)
    expected_reference_points = {
        {"a1", {3.37316107749939, 3.37316107749939, 0}},
        {"a2", {0, 3.37316107749939, 3.37316107749939}},
        {"a3", {3.37316107749939, 0, 3.37316107749939}},
        {"cmmm", {-3.37316107749939, -3.37316107749939, -3.37316107749939}},
        {"cmmp", {0, -3.37316107749939, 0}},
        {"cmpm", {-3.37316107749939, 0, 0}},
        {"cmpp", {0, 0, 3.37316107749939}},
        {"cpmm", {0, 0, -3.37316107749939}},
        {"cpmp", {3.37316107749939, 0, 0}},
        {"cppm", {0, 3.37316107749939, 0}},
        {"cppp", {3.37316107749939, 3.37316107749939, 3.37316107749939}},
        {"f1m", {-1.686580538749695, -1.686580538749695, 0}},
        {"f1p", {1.686580538749695, 1.686580538749695, 0}},
        {"f2m", {0, -1.686580538749695, -1.686580538749695}},
        {"f2p", {0, 1.686580538749695, 1.686580538749695}},
        {"f3m", {-1.686580538749695, 0, -1.686580538749695}},
        {"f3p", {1.686580538749695, 0, 1.686580538749695}},
        {"ion1", {0, 0, 0}},
        {"ion2", {1.686580538749695, 1.686580538749695, 1.686580538749695}},
        {"r1", {3.37316107749939, 3.37316107749939, 0}},
        {"r2", {0, 3.37316107749939, 3.37316107749939}},
        {"r3", {3.37316107749939, 0, 3.37316107749939}},
        {"zero", {0, 0, 0}},
    };
  else
    expected_reference_points = {
        {"a1", {3.37316115, 3.37316115, 0}},
        {"a2", {0, 3.37316115, 3.37316115}},
        {"a3", {3.37316115, 0, 3.37316115}},
        {"cmmm", {-3.37316115, -3.37316115, -3.37316115}},
        {"cmmp", {0, -3.37316115, 0}},
        {"cmpm", {-3.37316115, 0, 0}},
        {"cmpp", {0, 0, 3.37316115}},
        {"cpmm", {0, 0, -3.37316115}},
        {"cpmp", {3.37316115, 0, 0}},
        {"cppm", {0, 3.37316115, 0}},
        {"cppp", {3.37316115, 3.37316115, 3.37316115}},
        {"f1m", {-1.686580575, -1.686580575, 0}},
        {"f1p", {1.686580575, 1.686580575, 0}},
        {"f2m", {0, -1.686580575, -1.686580575}},
        {"f2p", {0, 1.686580575, 1.686580575}},
        {"f3m", {-1.686580575, 0, -1.686580575}},
        {"f3p", {1.686580575, 0, 1.686580575}},
        {"ion1", {0, 0, 0}},
        {"ion2", {1.68658058, 1.68658058, 1.68658058}},
        {"r1", {3.37316115, 3.37316115, 0}},
        {"r2", {0, 3.37316115, 3.37316115}},
        {"r3", {3.37316115, 0, 3.37316115}},
        {"zero", {0, 0, 0}},
    };
  return expected_reference_points;
}

TEST_CASE("ReferencePoints::DefaultConstruction", "[estimators]")
{
  auto rpi = makeTestRPI();
  auto ppr = makePsets();
  NEReferencePoints ref_points(rpi, ppr.pset, ppr.ref_psets);

  if (generate_test_data)
  {
    testing::TestableNEReferencePoints tref_points(ref_points);
    std::cout << "expected_reference_points" << tref_points;
  }

  typename NEReferencePoints::Points expected_reference_points;
  expected_reference_points = expectedReferencePoints();

  for (auto& [key, value] : ref_points.get_points())
  {
    bool coords_match = approxEquality(expected_reference_points[key], value);
    CHECK(coords_match);
  }
}

TEST_CASE("ReferencePoints::Construction", "[estimators]")
{
  auto rpi = makeTestRPI();
  auto ppr = makePsets();
  NEReferencePoints ref_points(rpi, ppr.pset, ppr.ref_psets);

  if (generate_test_data)
  {
    testing::TestableNEReferencePoints tref_points(ref_points);
    std::cout << "expected_reference_points" << tref_points;
  }

  typename NEReferencePoints::Points expected_reference_points;
  expected_reference_points = expectedReferencePoints();

  for (auto& [key, value] : ref_points.get_points())
  {
    bool coords_match = approxEquality(expected_reference_points[key], value);
    CHECK(coords_match);
  }
}

TEST_CASE("ReferencePoints::Description", "[estimators]")
{
  auto rpi = makeTestRPI();
  auto ppr = makePsets();
  NEReferencePoints ref_points(rpi, ppr.pset, ppr.ref_psets);

  std::ostringstream ostr_stream;
  ref_points.write_description(ostr_stream, "  ");

  std::string expected_description;
  if constexpr (std::is_same_v<NEReferencePoints::Real, double>)
    expected_description = R"(  reference_points
    a1:         3.37316115        3.37316115                 0
    a2:                  0        3.37316115        3.37316115
    a3:         3.37316115                 0        3.37316115
    cmmm:        -3.37316115       -3.37316115       -3.37316115
    cmmp:                  0       -3.37316115                 0
    cmpm:        -3.37316115                 0                 0
    cmpp:                  0                 0        3.37316115
    cpmm:                  0                 0       -3.37316115
    cpmp:         3.37316115                 0                 0
    cppm:                  0        3.37316115                 0
    cppp:         3.37316115        3.37316115        3.37316115
    f1m:       -1.686580575      -1.686580575                 0
    f1p:        1.686580575       1.686580575                 0
    f2m:                  0      -1.686580575      -1.686580575
    f2p:                  0       1.686580575       1.686580575
    f3m:       -1.686580575                 0      -1.686580575
    f3p:        1.686580575                 0       1.686580575
    ion1:                  0                 0                 0
    ion2:         1.68658058        1.68658058        1.68658058
    r1:         3.37316115        3.37316115                 0
    r2:                  0        3.37316115        3.37316115
    r3:         3.37316115                 0        3.37316115
    zero:                  0                 0                 0
  end reference_points
)";
  else
    expected_description = R"(  reference_points
    a1:        3.373161077       3.373161077                 0
    a2:                  0       3.373161077       3.373161077
    a3:        3.373161077                 0       3.373161077
    cmmm:       -3.373161077      -3.373161077      -3.373161077
    cmmp:                  0      -3.373161077                 0
    cmpm:       -3.373161077                 0                 0
    cmpp:                  0                 0       3.373161077
    cpmm:                  0                 0      -3.373161077
    cpmp:        3.373161077                 0                 0
    cppm:                  0       3.373161077                 0
    cppp:        3.373161077       3.373161077       3.373161077
    f1m:       -1.686580539      -1.686580539                 0
    f1p:        1.686580539       1.686580539                 0
    f2m:                  0      -1.686580539      -1.686580539
    f2p:                  0       1.686580539       1.686580539
    f3m:       -1.686580539                 0      -1.686580539
    f3p:        1.686580539                 0       1.686580539
    ion1:                  0                 0                 0
    ion2:        1.686580539       1.686580539       1.686580539
    r1:        3.373161077       3.373161077                 0
    r2:                  0       3.373161077       3.373161077
    r3:        3.373161077                 0       3.373161077
    zero:                  0                 0                 0
  end reference_points
)";
  CHECK(ostr_stream.str() == expected_description);
  // This subclass and function are used to generate the test data and may be useful for further test cases in future.
  testing::TestableNEReferencePoints test_ref_points(ref_points);
  std::ostringstream ostr_testing_stream;
  ostr_testing_stream << test_ref_points;
  std::string expected_testable_description;
  if constexpr (std::is_same_v<NEReferencePoints::Real, double>)
  {
    std::cout << "NEReferencePoints::Real == double\n";
    expected_testable_description = R"({
 {"a1", {      3.37316115,3.37316115,0}},
 {"a2", {               0,3.37316115,3.37316115}},
 {"a3", {      3.37316115,0,3.37316115}},
 {"cmmm", {     -3.37316115,-3.37316115,-3.37316115}},
 {"cmmp", {               0,-3.37316115,0}},
 {"cmpm", {     -3.37316115,0,0}},
 {"cmpp", {               0,0,3.37316115}},
 {"cpmm", {               0,0,-3.37316115}},
 {"cpmp", {      3.37316115,0,0}},
 {"cppm", {               0,3.37316115,0}},
 {"cppp", {      3.37316115,3.37316115,3.37316115}},
 {"f1m", {    -1.686580575,-1.686580575,0}},
 {"f1p", {     1.686580575,1.686580575,0}},
 {"f2m", {               0,-1.686580575,-1.686580575}},
 {"f2p", {               0,1.686580575,1.686580575}},
 {"f3m", {    -1.686580575,0,-1.686580575}},
 {"f3p", {     1.686580575,0,1.686580575}},
 {"ion1", {               0,0,0}},
 {"ion2", {      1.68658058,1.68658058,1.68658058}},
 {"r1", {      3.37316115,3.37316115,0}},
 {"r2", {               0,3.37316115,3.37316115}},
 {"r3", {      3.37316115,0,3.37316115}},
 {"zero", {               0,0,0}},
};
)";
  }
  else
  {
    std::cout << "NEReferencePoints::Real == float\n";
    expected_testable_description = R"({
 {"a1", {3.37316107749939,3.37316107749939,0}},
 {"a2", {               0,3.37316107749939,3.37316107749939}},
 {"a3", {3.37316107749939,0,3.37316107749939}},
 {"cmmm", {-3.37316107749939,-3.37316107749939,-3.37316107749939}},
 {"cmmp", {               0,-3.37316107749939,0}},
 {"cmpm", {-3.37316107749939,0,0}},
 {"cmpp", {               0,0,3.37316107749939}},
 {"cpmm", {               0,0,-3.37316107749939}},
 {"cpmp", {3.37316107749939,0,0}},
 {"cppm", {               0,3.37316107749939,0}},
 {"cppp", {3.37316107749939,3.37316107749939,3.37316107749939}},
 {"f1m", {-1.686580538749695,-1.686580538749695,0}},
 {"f1p", {1.686580538749695,1.686580538749695,0}},
 {"f2m", {               0,-1.686580538749695,-1.686580538749695}},
 {"f2p", {               0,1.686580538749695,1.686580538749695}},
 {"f3m", {-1.686580538749695,0,-1.686580538749695}},
 {"f3p", {1.686580538749695,0,1.686580538749695}},
 {"ion1", {               0,0,0}},
 {"ion2", {1.686580538749695,1.686580538749695,1.686580538749695}},
 {"r1", {3.37316107749939,3.37316107749939,0}},
 {"r2", {               0,3.37316107749939,3.37316107749939}},
 {"r3", {3.37316107749939,0,3.37316107749939}},
 {"zero", {               0,0,0}},
};
)";
  }
  CHECK(ostr_testing_stream.str() == expected_testable_description);
}

TEST_CASE("ReferencePoints::HDF5", "[estimators]")
{
  auto rpi = makeTestRPI();
  auto ppr = makePsets();
  NEReferencePoints ref_points(rpi, ppr.pset, ppr.ref_psets);

  hdf_archive hd;
  std::string test_file{"reference_points_test.hdf"};
  bool okay = hd.create(test_file);
  REQUIRE(okay);

  ref_points.write(hd);

  hd.close();

  hdf_archive hd_read;
  bool okay_read = hd.open(test_file);

  hd.push("reference_points");

  typename NEReferencePoints::Points expected_reference_points;
  expected_reference_points = expectedReferencePoints();

  for (auto& map_entry : expected_reference_points)
  {
    std::string key{map_entry.first};
    NEReferencePoints::Point point;
    hd.readEntry(point, key);
    CHECK(approxEquality(point, map_entry.second));
  }
  hd.close();
}

} // namespace qmcplusplus

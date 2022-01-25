//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Lab
//
// File created by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Lab
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "MomentumDistribution.h"
#include "ValidOneBodyDensityMatricesInput.h"
#include "InvalidOneBodyDensityMatricesInput.h"
#include "EstimatorTesting.h"
#include "ParticleSet.h"
#include "TrialWaveFunction.h"
#include "OhmmsData/Libxml2Doc.h"
#include "Message/UniformCommunicateError.h"
#include "Particle/tests/MinimalParticlePool.h"
#include "QMCWaveFunctions/tests/MinimalWaveFunctionPool.h"
#include "Utilities/StdRandom.h"
#include "Utilities/StlPrettyPrint.hpp"

#include <stdio.h>
#include <sstream>


namespace qmcplusplus
{
using RealType = QMCTraits::RealType;
using PosType  = QMCTraits::PosType;

namespace testing
{
/** class to preserve access control in MomentumDistribution
 */
class MomentumDistributionTests
{
public:
  void testCopyConstructor(const MomentumDistribution& md)
  {
    MomentumDistribution md2(md);

    CHECK(md2.M == md.M);
    CHECK(md2.twist[0] == Approx(md.twist[0]));
    CHECK(md2.twist[1] == Approx(md.twist[1]));
    CHECK(md2.twist[2] == Approx(md.twist[2]));
    CHECK(md2.kPoints.size() == md.kPoints.size());
    CHECK(md.data_ != md2.data_);
  }
};
} // namespace testing

TEST_CASE("MomentumDistribution::MomentumDistribution", "[estimators]")
{
  // clang-format: off
  const char* xml = R"(
<estimator type="MomentumDistribution" name="nofk" samples="5" kmax="3"/>
)";
  // clang-format: on

  // Read xml into input object
  Libxml2Document doc;
  bool okay = doc.parseFromString(xml);
  if (!okay)
    throw std::runtime_error("cannot parse OneBodyDensitMatricesInput section");
  xmlNodePtr node = doc.getRoot();
  MomentumDistributionInput mdi;
  mdi.readXML(node);

  // Instantiate other dependencies (internal QMCPACK objects)
  auto lattice = testing::makeTestLattice();
  Communicate* comm;
  comm = OHMMS::Controller;
  outputManager.pause();
  MinimalParticlePool mpp;
  ParticleSetPool particle_pool = mpp(comm);
  MinimalWaveFunctionPool wfp;
  WaveFunctionPool wavefunction_pool = wfp(comm, particle_pool);
  auto& pset                         = *(particle_pool.getParticleSet("e"));
  auto& wf_factory                   = *(wavefunction_pool.getWaveFunctionFactory("wavefunction"));
  DataLocality dl                    = DataLocality::crowd;

  // Build from input
  MomentumDistribution md(std::move(mdi), pset.getTotalNum(), pset.getTwist(), pset.getLattice(), dl);

  CHECK(md.M == 5);
  CHECK(md.twist[0] == Approx(0.0));
  CHECK(md.twist[1] == Approx(0.0));
  CHECK(md.twist[2] == Approx(0.0));
  CHECK(md.kPoints.size() == 27);

  // make sure there is something in mds data
  using namespace testing;
  OEBAccessor oeba(md);
  oeba[0] = 1.0;

  MomentumDistributionTests mdt;
  mdt.testCopyConstructor(md);

  outputManager.resume();
}


TEST_CASE("MomentumDistribution::accumulate", "[estimators]")
{
  using MCPWalker = OperatorEstBase::MCPWalker;

  // clang-format: off
  const char* xml = R"(
<estimator type="MomentumDistribution" name="nofk" samples="5" kmax="3"/>
)";
  // clang-format: on

  // Read xml into input object
  Libxml2Document doc;
  bool okay = doc.parseFromString(xml);
  if (!okay)
    throw std::runtime_error("cannot parse OneBodyDensitMatricesInput section");
  xmlNodePtr node = doc.getRoot();
  MomentumDistributionInput mdi;
  mdi.readXML(node);

  // Instantiate other dependencies (internal QMCPACK objects)
  auto lattice = testing::makeTestLattice();
  Communicate* comm;
  comm = OHMMS::Controller;
  outputManager.pause();
  MinimalParticlePool mpp;
  ParticleSetPool particle_pool = mpp(comm);
  MinimalWaveFunctionPool wfp;
  WaveFunctionPool wavefunction_pool = wfp(comm, particle_pool);
  auto& pset                         = *(particle_pool.getParticleSet("e"));
  auto& wf_factory                   = *(wavefunction_pool.getWaveFunctionFactory("wavefunction"));
  DataLocality dl                    = DataLocality::crowd;

  // Setup particleset
  pset.R = ParticleSet::ParticlePos{{1.751870349, 4.381521229, 2.865202269}, {3.244515371, 4.382273176, 4.21105285},
                                    {3.000459944, 3.329603408, 4.265030556}, {3.748660329, 3.63420622, 5.393637791},
                                    {3.033228526, 3.391869137, 4.654413566}, {3.114198787, 2.654334594, 5.231075822},
                                    {3.657151589, 4.883870516, 4.201243939}, {2.97317591, 4.245644974, 4.284564732}};

  // Build from input
  MomentumDistribution md(std::move(mdi), pset.getTotalNum(), pset.getTwist(), pset.getLattice(), dl);

  // Test accumulate

  //   Setup walker, particleset, wavefunction ref vectors
  //     Make clones
  std::vector<MCPWalker> walkers;
  int nwalkers = 4;
  for (int iw = 0; iw < nwalkers; ++iw)
    walkers.emplace_back(8);

  std::vector<ParticleSet> psets;
  for (int iw = 0; iw < nwalkers; ++iw)
    psets.emplace_back(pset);

  auto& trial_wavefunction = *(wavefunction_pool.getPrimary());
  std::vector<UPtr<TrialWaveFunction>> wfns(nwalkers);
  for (int iw = 0; iw < nwalkers; ++iw)
    wfns[iw] = trial_wavefunction.makeClone(psets[iw]);

  //     Initialize walker, pset, wfn
  for (int iw = 0; iw < nwalkers; ++iw)
  {
    auto& walker = walkers[iw];
    auto& pset   = psets[iw];
    pset.update(true);
    pset.donePbyP();
    wfns[iw]->evaluateLog(pset);
    pset.saveWalker(walker);
  }

  //     Create ref vectors
  auto ref_walkers = makeRefVector<MCPWalker>(walkers);
  auto ref_psets   = makeRefVector<ParticleSet>(psets);
  auto ref_wfns    = convertUPtrToRefVector(wfns);

  //   Setup RNG
  RandomGenerator rng;

  //   Perform accumulate
  md.accumulate(ref_walkers, ref_psets, ref_wfns, rng);

  //   Check data
  std::vector<RealType>& data = md.get_data();

  using Data = MomentumDistribution::Data;
  Data ref_data;

  ref_data = {3.92261216,    -5.752141485, 4.78276286,    8.307662762,   -5.130834919, 0.08942598353, 0.9716326509,
              21.82310933,   -9.177741101, -0.2024849597, -2.520417488,  -9.470020717, -9.4969045,    3.866360129,
              -9.4969045,    -9.470020717, -2.520417488,  -0.2024849597, -9.177741101, 21.82310933,   0.9716326509,
              0.08942598353, -5.130834919, 8.307662762,   4.78276286,    -5.752141485, 3.92261216};

  //std::cout<<"\n\n\nn(k) data:\n{";
  //for(int i=0;i<data.size();++i)
  //  std::cout<<data[i]<<", ";
  //std::cout<<"}\n\n\n";

  for (size_t id = 0; id < ref_data.size(); ++id)
  {
#ifdef MIXED_PRECISION
    CHECK(data[id] == Approx(ref_data[id]).epsilon(2.e-05));
#else
    // default Catch2 epsilon std::numeric_limits<float>::epsilon()*100
    // set value for x86_64
    CHECK(data[id] == Approx(ref_data[id]).epsilon(1.192092896e-05));
#endif
  }

  outputManager.resume();
}


} // namespace qmcplusplus

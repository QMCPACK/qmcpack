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
  auto lattice     = testing::makeTestLattice();
  Communicate* comm;
  comm = OHMMS::Controller;
  outputManager.pause();
  MinimalParticlePool mpp;
  ParticleSetPool particle_pool = mpp(comm);
  MinimalWaveFunctionPool wfp;
  WaveFunctionPool wavefunction_pool = wfp(comm, particle_pool);
  auto& pset                         = *(particle_pool.getParticleSet("e"));
  auto& wf_factory                   = *(wavefunction_pool.getWaveFunctionFactory("wavefunction"));
  DataLocality dl = DataLocality::crowd;
  
  // Build from input
  MomentumDistribution md(std::move(mdi), pset.getTotalNum(), pset.getTwist(), 
                          pset.Lattice, dl);
  
  CHECK(md.M==5);
  CHECK(md.twist[0]==Approx(0.0));
  CHECK(md.twist[1]==Approx(0.0));
  CHECK(md.twist[2]==Approx(0.0));
  CHECK(md.kPoints.size()==27);
  
  // Copy constructor
  MomentumDistribution md2(md);
  
  CHECK(md2.M==5);
  CHECK(md2.twist[0]==Approx(0.0));
  CHECK(md2.twist[1]==Approx(0.0));
  CHECK(md2.twist[2]==Approx(0.0));
  CHECK(md2.kPoints.size()==27);
  
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
  auto lattice     = testing::makeTestLattice();
  Communicate* comm;
  comm = OHMMS::Controller;
  outputManager.pause();
  MinimalParticlePool mpp;
  ParticleSetPool particle_pool = mpp(comm);
  MinimalWaveFunctionPool wfp;
  WaveFunctionPool wavefunction_pool = wfp(comm, particle_pool);
  auto& pset                         = *(particle_pool.getParticleSet("e"));
  auto& wf_factory                   = *(wavefunction_pool.getWaveFunctionFactory("wavefunction"));
  DataLocality dl = DataLocality::crowd;
  
  // Build from input
  MomentumDistribution md(std::move(mdi), pset.getTotalNum(), pset.getTwist(), 
                          pset.Lattice, dl);
  
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
  RandomGenerator_t rng;

  //   Perform accumulate
  md.accumulate(ref_walkers, ref_psets, ref_wfns, rng);

  //   Check data
  std::vector<RealType>& data = md.get_data_ref();

  using Data = MomentumDistribution::Data::element_type;
  Data ref_data;

  ref_data = {46.29362653, 3.159647865, -7.619226601, 27.60997229, -0.2988776457, -37.49398908, -8.803560748, -37.01514494, 16.00502823, 18.17654053, 1.090546098, 12.9452013, -28.36356552, 1.734146548, -28.36356552, 12.9452013, 1.090546098, 18.17654053, 16.00502823, -37.01514494, -8.803560748, -37.49398908, -0.2988776457, 27.60997229, -7.619226601, 3.159647865, 46.29362653 };

  for (size_t id = 0; id < ref_data.size(); ++id)
    CHECK(data[id] == Approx(ref_data[id]));

  //std::cout<<"\n\n\nn(k) data:\n{";
  //for(int i=0;i<data.size();++i)
  //  std::cout<<data[i]<<", ";
  //std::cout<<"}\n\n\n";

  outputManager.resume();

}



} // namespace qmcplusplus

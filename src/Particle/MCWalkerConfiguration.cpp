//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jordan E. Vincent, University of Illinois at Urbana-Champaign
//                    Bryan Clark, bclark@Princeton.edu, Princeton University
//                    Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Cynthia Gu, zg1@ornl.gov, Oak Ridge National Laboratory
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "MCWalkerConfiguration.h"
#include "ParticleBase/RandomSeqGenerator.h"
#include "Message/Communicate.h"
#include "Message/CommOperators.h"
#include "Utilities/IteratorUtility.h"
#include "LongRange/StructFact.h"
#include "Particle/HDFWalkerOutput.h"
#include "Particle/MCSample.h"
#include "Particle/Reptile.h"
#include "hdf/hdf_hyperslab.h"
#include "hdf/HDFVersion.h"
#include <map>

namespace qmcplusplus
{
MCWalkerConfiguration::MCWalkerConfiguration(const SimulationCell& simulation_cell, const DynamicCoordinateKind kind)
    : ParticleSet(simulation_cell, kind),
      ReadyForPbyP(false),
      UpdateMode(Update_Walker),
      reptile(0),
      Polymer(0)
{
}

MCWalkerConfiguration::MCWalkerConfiguration(const MCWalkerConfiguration& mcw)
    : ParticleSet(mcw),
      ReadyForPbyP(false),
      UpdateMode(Update_Walker),
      Polymer(0)
{
  samples.clearEnsemble();
  samples.setMaxSamples(mcw.getMaxSamples());
  GlobalNumWalkers = mcw.GlobalNumWalkers;
  WalkerOffsets    = mcw.WalkerOffsets;
  Properties       = mcw.Properties;
}

MCWalkerConfiguration::~MCWalkerConfiguration() = default;

void MCWalkerConfiguration::createWalkers(int n)
{
  const int old_nw = getActiveWalkers();
  WalkerConfigurations::createWalkers(n, TotalNum);
  // no pre-existing walkers, need to initialized based on particleset.
  if (old_nw == 0)
    for (auto& awalker : WalkerList)
    {
      awalker->R     = R;
      awalker->spins = spins;
    }
  resizeWalkerHistories();
}


void MCWalkerConfiguration::resize(int numWalkers, int numPtcls)
{
  if (TotalNum && WalkerList.size())
    app_warning() << "MCWalkerConfiguration::resize cleans up the walker list." << std::endl;
  const int old_nw = getActiveWalkers();
  ParticleSet::resize(unsigned(numPtcls));
  WalkerConfigurations::resize(numWalkers, TotalNum);
  // no pre-existing walkers, need to initialized based on particleset.
  if (old_nw == 0)
    for (auto& awalker : WalkerList)
    {
      awalker->R     = R;
      awalker->spins = spins;
    }
}

/** Make Metropolis move to the walkers and save in a temporary array.
 * @param it the iterator of the first walker to work on
 * @param tauinv  inverse of the time step
 *
 * R + D + X
 */
void MCWalkerConfiguration::sample(iterator it, RealType tauinv)
{
  throw std::runtime_error("MCWalkerConfiguration::sample obsolete");
  //  makeGaussRandom(R);
  //  R *= tauinv;
  //  R += (*it)->R + (*it)->Drift;
}

//void MCWalkerConfiguration::clearAuxDataSet() {
//  UpdateMode=Update_Particle;
//  int nbytes=128*TotalNum*sizeof(RealType);//could be pagesize
//  if(WalkerList.size())//check if capacity is bigger than the estimated one
//    nbytes = (WalkerList[0]->DataSet.capacity()>nbytes)?WalkerList[0]->DataSet.capacity():nbytes;
//  iterator it(WalkerList.begin());
//  iterator it_end(WalkerList.end());
//  while(it!=it_end) {
//    (*it)->DataSet.clear();
//    //CHECK THIS WITH INTEL 10.1
//    //(*it)->DataSet.reserve(nbytes);
//    ++it;
//  }
//  ReadyForPbyP = true;
//}
//
//bool MCWalkerConfiguration::createAuxDataSet(int nfield) {
//
//  if(ReadyForPbyP) return false;
//
//  ReadyForPbyP=true;
//  UpdateMode=Update_Particle;
//  iterator it(WalkerList.begin());
//  iterator it_end(WalkerList.end());
//  while(it!=it_end) {
//    (*it)->DataSet.reserve(nfield); ++it;
//  }
//
//  return true;
//}


/** reset the Property container of all the walkers
 */
void MCWalkerConfiguration::resetWalkerProperty(int ncopy)
{
  int m(PropertyList.size());
  app_log() << "  Resetting Properties of the walkers " << ncopy << " x " << m << std::endl;
  try
  {
    Properties.resize(ncopy, m);
  }
  catch (std::domain_error& de)
  {
    app_error() << de.what() << '\n'
                << "This is likely because some object has attempted to add walker properties\n"
                << " in excess of WALKER_MAX_PROPERTIES.\n"
                << "build with cmake ... -DWALKER_MAX_PROPERTIES=at_least_properties_required" << std::endl;
    APP_ABORT("Fatal Exception");
  }

  iterator it(WalkerList.begin()), it_end(WalkerList.end());
  while (it != it_end)
  {
    (*it)->resizeProperty(ncopy, m);
    (*it)->Weight = 1;
    ++it;
  }
  resizeWalkerHistories();
}

void MCWalkerConfiguration::resizeWalkerHistories()
{
  //using std::vector<std::vector<RealType> > is too costly.
  int np = PropertyHistory.size();
  if (np)
    for (int iw = 0; iw < WalkerList.size(); ++iw)
      WalkerList[iw]->PropertyHistory = PropertyHistory;
  np = PHindex.size();
  if (np)
    for (int iw = 0; iw < WalkerList.size(); ++iw)
      WalkerList[iw]->PHindex = PHindex;
  ;
}

/** allocate the SampleStack
 * @param n number of samples per thread
 */
void MCWalkerConfiguration::setNumSamples(int n)
{
  samples.clearEnsemble();
  samples.setMaxSamples(n);
}

/** save the current walkers to SampleStack
 */
void MCWalkerConfiguration::saveEnsemble() { saveEnsemble(WalkerList.begin(), WalkerList.end()); }

/** save the [first,last) walkers to SampleStack
 */
void MCWalkerConfiguration::saveEnsemble(iterator first, iterator last)
{
  for (; first != last; first++)
  {
    samples.appendSample(MCSample(**first));
  }
}
/** load a single sample from SampleStack
 */
void MCWalkerConfiguration::loadSample(ParticleSet& pset, size_t iw) const { samples.loadSample(pset, iw); }

/** load SampleStack to WalkerList
 */
void MCWalkerConfiguration::loadEnsemble()
{
  using WP     = WalkerProperties::Indexes;
  int nsamples = std::min(samples.getMaxSamples(), samples.getNumSamples());
  if (samples.empty() || nsamples == 0)
    return;
  Walker_t::PropertyContainer_t prop(1, PropertyList.size(), 1, WP::MAXPROPERTIES);
  WalkerList.resize(nsamples);
  for (int i = 0; i < nsamples; ++i)
  {
    auto awalker = std::make_unique<Walker_t>(TotalNum);
    awalker->Properties.copy(prop);
    samples.getSample(i).convertToWalker(*awalker);
    WalkerList[i] = std::move(awalker);
  }
  resizeWalkerHistories();
  samples.clearEnsemble();
}

bool MCWalkerConfiguration::dumpEnsemble(std::vector<MCWalkerConfiguration*>& others,
                                         HDFWalkerOutput& out,
                                         int np,
                                         int nBlock)
{
  WalkerConfigurations wctemp;
  for (auto* mcwc : others)
  {
    const auto& astack(mcwc->getSampleStack());
    const size_t sample_size = std::min(mcwc->getMaxSamples(), mcwc->numSamples());
    for (int j = 0; j < sample_size; ++j)
    {
      const auto& sample     = astack.getSample(j);
      const size_t num_ptcls = sample.getNumPtcls();
      auto awalker           = std::make_unique<Walker_t>(num_ptcls);
      sample.convertToWalker(*awalker);
      wctemp.push_back(std::move(awalker));
    }
  }
  const int w = wctemp.getActiveWalkers();
  if (w == 0)
    return false;

  // The following code assumes the same amount of active walkers on all the MPI ranks
  std::vector<int> nwoff(np + 1, 0);
  for (int ip = 0; ip < np; ++ip)
    nwoff[ip + 1] = nwoff[ip] + w;
  wctemp.setGlobalNumWalkers(nwoff[np]);
  wctemp.setWalkerOffsets(nwoff);
  out.dump(wctemp, nBlock);
  return true;
}

int MCWalkerConfiguration::getMaxSamples() const { return samples.getMaxSamples(); }

void MCWalkerConfiguration::loadEnsemble(std::vector<MCWalkerConfiguration*>& others, bool doclean)
{
  using WP = WalkerProperties::Indexes;
  std::vector<int> off(others.size() + 1, 0);
  for (int i = 0; i < others.size(); ++i)
  {
    off[i + 1] = off[i] + std::min(others[i]->getMaxSamples(), others[i]->numSamples());
  }
  int nw_tot = off.back();
  if (nw_tot)
  {
    Walker_t::PropertyContainer_t prop(1, PropertyList.size(), 1, WP::MAXPROPERTIES);
    while (WalkerList.size())
      pop_back();
    WalkerList.resize(nw_tot);
    for (int i = 0; i < others.size(); ++i)
    {
      SampleStack& astack(others[i]->getSampleStack());
      for (int j = 0, iw = off[i]; iw < off[i + 1]; ++j, ++iw)
      {
        auto awalker = std::make_unique<Walker_t>(TotalNum);
        awalker->Properties.copy(prop);
        astack.getSample(j).convertToWalker(*awalker);
        WalkerList[iw] = std::move(awalker);
      }
      if (doclean)
        others[i]->clearEnsemble();
    }
  }
  if (doclean)
    resizeWalkerHistories();
}

void MCWalkerConfiguration::clearEnsemble() { samples.clearEnsemble(); }

} // namespace qmcplusplus

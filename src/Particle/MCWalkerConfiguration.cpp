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

#ifdef QMC_CUDA
#include "Particle/accept_kernel.h"
#endif

namespace qmcplusplus
{
MCWalkerConfiguration::MCWalkerConfiguration(const DynamicCoordinateKind kind)
    : ParticleSet(kind),
#ifdef QMC_CUDA
      RList_GPU("MCWalkerConfiguration::RList_GPU"),
      GradList_GPU("MCWalkerConfiguration::GradList_GPU"),
      LapList_GPU("MCWalkerConfiguration::LapList_GPU"),
      Rnew_GPU("MCWalkerConfiguration::Rnew_GPU"),
      NLlist_GPU("MCWalkerConfiguration::NLlist_GPU"),
      iatList_GPU("iatList_GPU"),
      AcceptList_GPU("MCWalkerConfiguration::AcceptList_GPU"),
#endif
      ReadyForPbyP(false),
      UpdateMode(Update_Walker),
      reptile(0),
      Polymer(0)
{
  //move to ParticleSet
  //initPropertyList();
}

MCWalkerConfiguration::MCWalkerConfiguration(const MCWalkerConfiguration& mcw)
    : ParticleSet(mcw),
#ifdef QMC_CUDA
      RList_GPU("MCWalkerConfiguration::RList_GPU"),
      GradList_GPU("MCWalkerConfiguration::GradList_GPU"),
      LapList_GPU("MCWalkerConfiguration::LapList_GPU"),
      Rnew_GPU("MCWalkerConfiguration::Rnew_GPU"),
      NLlist_GPU("MCWalkerConfiguration::NLlist_GPU"),
      iatList_GPU("iatList_GPU"),
      AcceptList_GPU("MCWalkerConfiguration::AcceptList_GPU"),
#endif
      ReadyForPbyP(false),
      UpdateMode(Update_Walker),
      Polymer(0)
{
  samples.clearEnsemble();
  samples.setMaxSamples(mcw.getMaxSamples());
  GlobalNumWalkers = mcw.GlobalNumWalkers;
  WalkerOffsets    = mcw.WalkerOffsets;
  Properties       = mcw.Properties;
  //initPropertyList();
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
  APP_ABORT("MCWalkerConfiguration::sample obsolete");
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

///** load SampleStack to WalkerList
// */
//void MCWalkerConfiguration::loadEnsemble(const Walker_t& wcopy)
//{
//  int nsamples=std::min(MaxSamples,CurSampleCount);
//  if(SampleStack.empty() || nsamples==0) return;
//
//  Walker_t::PropertyContainer_t prop(1,PropertyList.size());
//
//  while(WalkerList.size()) pop_back();
//  WalkerList.resize(nsamples);
//
//  for(int i=0; i<nsamples; ++i)
//  {
//    Walker_t* awalker=new Walker_t(TotalNum);
//    awalker->Properties.copy(prop);
//    SampleStack[i]->get(*awalker);
//    WalkerList[i]=awalker;
//  }
//  resizeWalkerHistories();
//  clearEnsemble();
//}
//
//void MCWalkerConfiguration::loadEnsemble(MCWalkerConfiguration& other)
//{
//  if(SampleStack.empty()) return;
//
//  Walker_t twalker(*WalkerList[0]);
//  for(int i=0; i<MaxSamples; ++i)
//  {
//    Walker_t* awalker=new Walker_t(twalker);
//    SampleStack[i]->get(*awalker);
//    other.WalkerList.push_back(awalker);
//  }
//
//  clearEnsemble();
//}

bool MCWalkerConfiguration::dumpEnsemble(std::vector<MCWalkerConfiguration*>& others,
                                         HDFWalkerOutput& out,
                                         int np,
                                         int nBlock)
{
  return samples.dumpEnsemble(others, out, np, nBlock);
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

#ifdef QMC_CUDA
void MCWalkerConfiguration::updateLists_GPU()
{
  int nw         = WalkerList.size();
  int NumSpecies = getSpeciesSet().TotalNum;
  if (Rnew_GPU.size() != nw * kblocksize)
  {
    Rnew_GPU.resize(nw * kblocksize);
    RhokLists_GPU.resize(NumSpecies);
    for (int isp = 0; isp < NumSpecies; isp++)
      RhokLists_GPU[isp].resize(nw);
    Rnew_host.resize(nw * kblocksize);
    Rnew.resize(nw * kblocksize);
    AcceptList_GPU.resize(nw);
    AcceptList_host.resize(nw);
    RList_GPU.resize(nw);
    GradList_GPU.resize(nw);
    LapList_GPU.resize(nw);
    DataList_GPU.resize(nw);
  }
  hostlist.resize(nw);
  hostlist_valueType.resize(nw);
  hostlist_AA.resize(nw);

  for (int iw = 0; iw < nw; iw++)
  {
    if (WalkerList[iw]->R_GPU.size() != R.size())
      std::cerr << "Error in R_GPU size for iw = " << iw << "!\n";
    hostlist[iw] = (CTS::RealType*)WalkerList[iw]->R_GPU.data();
  }
  RList_GPU = hostlist;

  for (int iw = 0; iw < nw; iw++)
  {
    if (WalkerList[iw]->Grad_GPU.size() != R.size())
      std::cerr << "Error in Grad_GPU size for iw = " << iw << "!\n";
    hostlist_valueType[iw] = (CTS::ValueType*)WalkerList[iw]->Grad_GPU.data();
  }
  GradList_GPU = hostlist_valueType;

  for (int iw = 0; iw < nw; iw++)
  {
    if (WalkerList[iw]->Lap_GPU.size() != R.size())
      std::cerr << "Error in Lap_GPU size for iw = " << iw << "!\n";
    hostlist_valueType[iw] = (CTS::ValueType*)WalkerList[iw]->Lap_GPU.data();
  }
  LapList_GPU = hostlist_valueType;

  for (int iw = 0; iw < nw; iw++)
    hostlist_valueType[iw] = WalkerList[iw]->cuda_DataSet.data();
  DataList_GPU = hostlist_valueType;

  for (int isp = 0; isp < NumSpecies; isp++)
  {
    for (int iw = 0; iw < nw; iw++)
      hostlist_AA[iw] = WalkerList[iw]->get_rhok_ptr(isp);
    RhokLists_GPU[isp] = hostlist_AA;
  }
}

void MCWalkerConfiguration::allocateGPU(size_t buffersize)
{
  int N    = WalkerList[0]->R.size();
  int Numk = 0;
  if (SK)
    Numk = SK->getKLists().numk;
  int NumSpecies = getSpeciesSet().TotalNum;
  for (int iw = 0; iw < WalkerList.size(); iw++)
  {
    Walker_t& walker = *(WalkerList[iw]);
    walker.resizeCuda(buffersize, NumSpecies, Numk);
  }
}


void MCWalkerConfiguration::copyWalkersToGPU(bool copyGrad)
{
  R_host.resize(WalkerList[0]->R.size());
  for (int iw = 0; iw < WalkerList.size(); iw++)
  {
    for (int i = 0; i < WalkerList[iw]->size(); i++)
      for (int dim = 0; dim < OHMMS_DIM; dim++)
        R_host[i][dim] = WalkerList[iw]->R[i][dim];
    WalkerList[iw]->R_GPU = R_host;
  }
  if (copyGrad)
    copyWalkerGradToGPU();
}

void MCWalkerConfiguration::copyWalkerGradToGPU()
{
  Grad_host.resize(WalkerList[0]->G.size());
  for (int iw = 0; iw < WalkerList.size(); iw++)
  {
    for (int i = 0; i < WalkerList[iw]->size(); i++)
      for (int dim = 0; dim < OHMMS_DIM; dim++)
        Grad_host[i][dim] = WalkerList[iw]->G[i][dim];
    WalkerList[iw]->Grad_GPU = Grad_host;
  }
}

void MCWalkerConfiguration::proposeMove_GPU(std::vector<PosType>& newPos, int iat)
{
  int nw = newPos.size();
  if (Rnew_host.size() < nw * kblocksize)
  {
    Rnew.resize(nw * kblocksize);
    Rnew_host.resize(nw * kblocksize);
  }
  // store things sequentially with k to make evaluation more straight-forward:
  //           k=0     k=1     k=kblocksize
  // Rnew = [0,..,nw|0,..,nw|...|0,..,nw]
  int offset = kcurr * nw;
  for (int i = 0; i < nw; i++)
  {
    for (int dim = 0; dim < OHMMS_DIM; dim++)
    {
      Rnew[i + offset][dim]      = newPos[i][dim];
      Rnew_host[i + offset][dim] = newPos[i][dim];
    }
  }
  if (kDelay)
  {
    kcurr  = (kcurr + 1) % kblocksize; // loop kcurr around every k blocks
    kstart = kblock * kblocksize;
    if (kcurr == 0)
      kblock++; // keep increasing kblock (even beyond available matrix blocks) - the update check takes care of self-consistency
    // only copy new position matrix when needed (when update is imminent)
    if (klinear)
    {
      Rnew_GPU.asyncCopy(&(Rnew_host[offset]), nw * kblocksize, offset, nw);
    }
    else if (kcurr == 0 || (kcurr + kblock * kblocksize >= getnat(iat)))
      Rnew_GPU.asyncCopy(Rnew_host);
  }
  else
    Rnew_GPU.asyncCopy(Rnew_host);
  CurrentParticle = iat;
}


void MCWalkerConfiguration::acceptMove_GPU(std::vector<bool>& toAccept, int k)
{
  if (AcceptList_host.size() < toAccept.size())
    AcceptList_host.resize(toAccept.size());
  for (int i = 0; i < toAccept.size(); i++)
    AcceptList_host[i] = (int)toAccept[i];
  AcceptList_GPU.asyncCopy(AcceptList_host);
  //   app_log() << "toAccept.size()        = " << toAccept.size() << std::endl;
  //   app_log() << "AcceptList_host.size() = " << AcceptList_host.size() << std::endl;
  //   app_log() << "AcceptList_GPU.size()  = " << AcceptList_GPU.size() << std::endl;
  //   app_log() << "WalkerList.size()      = " << WalkerList.size() << std::endl;
  //   app_log() << "Rnew_GPU.size()        = " << Rnew_GPU.size() << std::endl;
  //   app_log() << "RList_GPU.size()       = " << RList_GPU.size() << std::endl;
  if (RList_GPU.size() != WalkerList.size())
    std::cerr << "Error in RList_GPU size.\n";
  if (Rnew_GPU.size() != WalkerList.size() * kblocksize)
    std::cerr << "Error in Rnew_GPU size.\n";
  if (AcceptList_GPU.size() != WalkerList.size())
    std::cerr << "Error in AcceptList_GPU size.\n";
  accept_move_GPU_cuda(RList_GPU.data(), (CUDA_PRECISION*)Rnew_GPU.data(), AcceptList_GPU.data(), CurrentParticle++,
                       WalkerList.size(), k);
}


void MCWalkerConfiguration::NLMove_GPU(std::vector<Walker_t*>& walkers,
                                       std::vector<PosType>& newpos,
                                       std::vector<int>& iat)
{
  int N = walkers.size();
  if (NLlist_GPU.size() < N)
  {
    NLlist_GPU.resize(N);
    NLlist_host.resize(N);
  }
  if (Rnew_GPU.size() < N)
  {
    Rnew_host.resize(N);
    Rnew_GPU.resize(N);
  }
  for (int iw = 0; iw < N; iw++)
  {
    Rnew_host[iw]   = newpos[iw];
    NLlist_host[iw] = (CUDA_PRECISION*)(walkers[iw]->R_GPU.data()) + OHMMS_DIM * iat[iw];
  }
  Rnew_GPU   = Rnew_host;
  NLlist_GPU = NLlist_host;
  NL_move_cuda(NLlist_GPU.data(), (CUDA_PRECISION*)Rnew_GPU.data(), N);
}
#endif

} // namespace qmcplusplus

//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020  QMCPACK developers.
//
// File developed by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//
// File created by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "Particle/SampleStack.h"
#include "Particle/HDFWalkerOutput.h"
#include "Utilities/IteratorUtility.h"

namespace qmcplusplus
{
/** store minimum Walker data for the next section
 */
struct MCSample
{
  using WP       = WalkerProperties::Indexes;
  using Walker_t = ParticleSet::Walker_t;

  ParticleSet::ParticlePos_t R;
  ParticleSet::ParticleGradient_t G;
  ParticleSet::ParticleLaplacian_t L;
  ParticleSet::RealType LogPsi, Sign, PE, KE;

  inline MCSample(const Walker_t& w) : R(w.R), G(w.G), L(w.L)
  {
    LogPsi = w.Properties(WP::LOGPSI);
    Sign   = w.Properties(WP::SIGN);
    PE     = w.Properties(WP::LOCALPOTENTIAL);
    KE     = w.Properties(WP::LOCALENERGY) - PE;
  }

  inline MCSample(int n)
  {
    R.resize(n);
    G.resize(n);
    L.resize(n);
  }

  inline void put(const Walker_t& w)
  {
    R      = w.R;
    G      = w.G;
    L      = w.L;
    LogPsi = w.Properties(WP::LOGPSI);
    Sign   = w.Properties(WP::SIGN);
    PE     = w.Properties(WP::LOCALPOTENTIAL);
    KE     = w.Properties(WP::LOCALENERGY) - PE;
  }

  inline void get(Walker_t& w) const
  {
    w.R                              = R;
    w.G                              = G;
    w.L                              = L;
    w.Properties(WP::LOGPSI)         = LogPsi;
    w.Properties(WP::SIGN)           = Sign;
    w.Properties(WP::LOCALPOTENTIAL) = PE;
    w.Properties(WP::LOCALENERGY)    = PE + KE;
  }
};


SampleStack::SampleStack() : total_num_(0), max_samples_(10), current_sample_count_(0) {}

/** allocate the SampleStack
 * @param n number of samples per thread
 */
void SampleStack::setMaxSamples(int n)
{
  clearEnsemble();
  max_samples_ = n;
  //do not add anything
  if (n == 0)
    return;
  sample_vector_.reserve(n);
  int nadd = n - sample_vector_.size();
  while (nadd > 0)
  {
    sample_vector_.push_back(new MCSample(total_num_));
    --nadd;
  }
}

/** save the [first,last) walkers to SampleStack
 */
void SampleStack::saveEnsemble(walker_iterator first, walker_iterator last)
{
  //safety check
  if (max_samples_ == 0)
    return;
  while ((first != last) && (current_sample_count_ < max_samples_))
  {
    sample_vector_[current_sample_count_]->put(**first);
    ++first;
    ++current_sample_count_;
  }
}


/** load a single sample from SampleStack
 */
void SampleStack::loadSample(ParticleSet::ParticlePos_t& Pos, size_t iw) const { Pos = sample_vector_[iw]->R; }

void SampleStack::getSample(unsigned int i, Walker_t& w) const { sample_vector_[i]->get(w); }

void SampleStack::putSample(unsigned int i, const Walker_t& w) { sample_vector_[i]->put(w); }

bool SampleStack::dumpEnsemble(std::vector<MCWalkerConfiguration*>& others, HDFWalkerOutput* out, int np, int nBlock)
{
  MCWalkerConfiguration wtemp;
  wtemp.resize(0, total_num_);
  wtemp.loadEnsemble(others, false);
  int w = wtemp.getActiveWalkers();
  if (w == 0)
    return false;
  std::vector<int> nwoff(np + 1, 0);
  for (int ip = 0; ip < np; ++ip)
    nwoff[ip + 1] = nwoff[ip] + w;
  wtemp.setGlobalNumWalkers(nwoff[np]);
  wtemp.setWalkerOffsets(nwoff);
  out->dump(wtemp, nBlock);
  return true;
}

void SampleStack::clearEnsemble()
{
  //delete_iter(SampleStack.begin(),SampleStack.end());
  for (int i = 0; i < sample_vector_.size(); ++i)
    if (sample_vector_[i])
      delete sample_vector_[i];
  sample_vector_.clear();
  max_samples_     = 0;
  current_sample_count_ = 0;
}


} // namespace qmcplusplus

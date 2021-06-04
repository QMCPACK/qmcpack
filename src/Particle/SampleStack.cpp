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


#include "SampleStack.h"
#include "Particle/HDFWalkerOutput.h"
#include "Particle/MCSample.h"
#include "Utilities/IteratorUtility.h"

namespace qmcplusplus
{
SampleStack::SampleStack() : total_num_(0), max_samples_(10), current_sample_count_(0) {}

/** allocate the SampleStack
 * @param n number of samples per rank
 * @param num_ranks number of ranks. Used to set global number of samples.
 */
void SampleStack::setMaxSamples(int n, int num_ranks)
{
  max_samples_        = n;
  global_num_samples_ = n * num_ranks;
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

MCSample& SampleStack::getSample(unsigned int i) const { return *sample_vector_[i]; }

void SampleStack::saveEnsemble(std::vector<MCSample>& walker_list)
{
  //safety check
  if (max_samples_ == 0)
    return;
  auto first = walker_list.begin();
  auto last  = walker_list.end();
  while ((first != last) && (current_sample_count_ < max_samples_))
  {
    *sample_vector_[current_sample_count_] = *first;
    ++first;
    ++current_sample_count_;
  }
}

void SampleStack::appendSample(MCSample&& sample)
{
  // Ignore samples in excess of the expected number of samples
  if (current_sample_count_ < max_samples_)
  {
    *sample_vector_[current_sample_count_] = std::move(sample);
    current_sample_count_++;
  }
}


/** load a single sample from SampleStack
 */
void SampleStack::loadSample(ParticleSet& pset, size_t iw) const
{
  pset.R     = sample_vector_[iw]->R;
  pset.spins = sample_vector_[iw]->spins;
}

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
  max_samples_          = 0;
  current_sample_count_ = 0;
}

SampleStack::~SampleStack() { clearEnsemble(); }

void SampleStack::resetSampleCount() { current_sample_count_ = 0; }


} // namespace qmcplusplus

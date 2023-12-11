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
#include "Utilities/IteratorUtility.h"

namespace qmcplusplus
{

/** allocate the SampleStack
 * @param n number of samples per rank
 * @param num_ranks number of ranks. Used to set global number of samples.
 */
void SampleStack::setMaxSamples(size_t n, size_t num_ranks)
{
  max_samples_        = n;
  global_num_samples_ = n * num_ranks;
  current_sample_count_ = std::min(current_sample_count_, max_samples_);
  sample_vector_.resize(n, MCSample(0));
}

const MCSample& SampleStack::getSample(size_t i) const { return sample_vector_[i]; }

void SampleStack::appendSample(MCSample&& sample)
{
  // Ignore samples in excess of the expected number of samples
  if (current_sample_count_ < max_samples_)
  {
    sample_vector_[current_sample_count_] = std::move(sample);
    current_sample_count_++;
  }
}

/** load a single sample from SampleStack
 */
void SampleStack::loadSample(ParticleSet& pset, size_t iw) const
{
  pset.R     = sample_vector_[iw].R;
  pset.spins = sample_vector_[iw].spins;
}

void SampleStack::clearEnsemble()
{
  sample_vector_.clear();
  current_sample_count_ = 0;
}

void SampleStack::resetSampleCount() { current_sample_count_ = 0; }


} // namespace qmcplusplus

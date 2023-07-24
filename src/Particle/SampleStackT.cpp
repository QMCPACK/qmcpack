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

#include "SampleStackT.h"

#include "Utilities/IteratorUtility.h"

namespace qmcplusplus
{

/** allocate the SampleStack
 * @param n number of samples per rank
 * @param num_ranks number of ranks. Used to set global number of samples.
 */
template<typename T>
void SampleStackT<T>::setMaxSamples(size_t n, size_t num_ranks)
{
  max_samples_          = n;
  global_num_samples_   = n * num_ranks;
  current_sample_count_ = std::min(current_sample_count_, max_samples_);
  sample_vector_.resize(n, MCSample(0));
}

template<typename T>
const MCSample& SampleStackT<T>::getSample(size_t i) const
{
  return sample_vector_[i];
}

template<typename T>
void SampleStackT<T>::appendSample(MCSample&& sample)
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
template<typename T>
void SampleStackT<T>::loadSample(ParticleSetT<T>& pset, size_t iw) const
{
  pset.R     = sample_vector_[iw].R;
  pset.spins = sample_vector_[iw].spins;
}

template<typename T>
void SampleStackT<T>::clearEnsemble()
{
  sample_vector_.clear();
  current_sample_count_ = 0;
}

template<typename T>
void SampleStackT<T>::resetSampleCount()
{
  current_sample_count_ = 0;
}

template class SampleStackT<double>;
template class SampleStackT<float>;
template class SampleStackT<std::complex<double>>;
template class SampleStackT<std::complex<float>>;

} // namespace qmcplusplus

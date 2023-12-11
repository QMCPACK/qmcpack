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

/** @file SampleStack.h
 * @brief Stores particle configurations for later use in DMC and wavefunction optimization
 *
 * Stores less temporary data than the buffer object.
 */

#ifndef QMCPLUSPLUS_SAMPLE_STACK_H
#define QMCPLUSPLUS_SAMPLE_STACK_H

#include <vector>
#include "Particle/ParticleSet.h"
#include "Particle/MCSample.h"
#include "Particle/Walker.h"
#include "Particle/WalkerConfigurations.h"

namespace qmcplusplus
{
class SampleStack
{
public:
  using PropertySetType = QMCTraits::PropertySetType;

  size_t getMaxSamples() const { return max_samples_; }

  bool empty() const { return sample_vector_.empty(); }

  const MCSample& getSample(size_t i) const;

  //@{save/load/clear function for optimization
  inline size_t getNumSamples() const { return current_sample_count_; }
  ///set the number of max samples per rank.
  void setMaxSamples(size_t n, size_t number_of_ranks = 1);
  /// Global number of samples is number of samples per rank * number of ranks
  size_t getGlobalNumSamples() const { return global_num_samples_; }
  /// load a single sample from SampleStack
  void loadSample(ParticleSet& pset, size_t iw) const;

  void appendSample(MCSample&& sample);

  ///clear the ensemble
  void clearEnsemble();
  //@}
  ///  Set the sample count to zero but preserve the storage
  void resetSampleCount();

private:
  size_t max_samples_{10};
  size_t current_sample_count_{0};
  size_t global_num_samples_{max_samples_};

  std::vector<MCSample> sample_vector_;
};


} // namespace qmcplusplus
#endif

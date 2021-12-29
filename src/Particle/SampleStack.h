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
#include "Particle/Walker.h"
#include "Particle/WalkerConfigurations.h"

namespace qmcplusplus
{
struct MCSample;

class SampleStack
{
public:
  using PropertySetType = QMCTraits::PropertySetType;

  SampleStack();

  void setTotalNum(int total_num) { total_num_ = total_num; }

  int getMaxSamples() const { return max_samples_; }

  bool empty() const { return sample_vector_.empty(); }

  MCSample& getSample(unsigned int i) const;

  //@{save/load/clear function for optimization
  inline int getNumSamples() const { return current_sample_count_; }
  ///set the number of max samples per rank.
  void setMaxSamples(int n, int number_of_ranks = 1);
  /// Global number of samples is number of samples per rank * number of ranks
  uint64_t getGlobalNumSamples() const { return global_num_samples_; }
  ///save the position of current walkers
  void saveEnsemble(std::vector<MCSample>& walker_list);
  /// load a single sample from SampleStack
  void loadSample(ParticleSet& pset, size_t iw) const;

  void appendSample(MCSample&& sample);

  ///clear the ensemble
  void clearEnsemble();
  //@}
  ///  Set the sample count to zero but preserve the storage
  void resetSampleCount();

  ~SampleStack();

private:
  int total_num_;
  int max_samples_;
  int current_sample_count_;
  uint64_t global_num_samples_;

  std::vector<MCSample*> sample_vector_;
};


} // namespace qmcplusplus
#endif

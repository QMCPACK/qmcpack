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

namespace qmcplusplus
{
class MCWalkerConfiguration;
class HDFWalkerOutput;
class MCSample;

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
  ///set the number of max samples
  void setMaxSamples(int n);
  ///save the position of current walkers
  void saveEnsemble(std::vector<MCSample>& walker_list);
  /// load a single sample from SampleStack
  void loadSample(ParticleSet::ParticlePos_t& Pos, size_t iw) const;

  void appendSample(MCSample&& sample);

  bool dumpEnsemble(std::vector<MCWalkerConfiguration*>& others, HDFWalkerOutput* out, int np, int nBlock);
  ///clear the ensemble
  void clearEnsemble();
  //@}

private:
  int total_num_;
  int max_samples_;
  int current_sample_count_;

  std::vector<MCSample*> sample_vector_;
};


} // namespace qmcplusplus
#endif

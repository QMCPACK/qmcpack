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
  using Walker_t        = ParticleSet::Walker_t;
  using WalkerList_t    = std::vector<Walker_t*>;
  using walker_iterator = WalkerList_t::iterator;
  using PropertySetType = QMCTraits::PropertySetType;

  SampleStack();

  int setTotalNum(int total_num) { total_num_ = total_num; }

  int getMaxSamples() const { return max_samples_; }

  void setMaxSamples(int max_samples) { max_samples_ = max_samples; }

  bool empty() const { return sample_vector_.empty(); }

  void putSample(unsigned int i, const Walker_t& w);

  void getSample(unsigned int i, Walker_t& w) const;

  //@{save/load/clear function for optimization
  inline int getNumSamples() const { return current_sample_count_; }
  ///set the number of max samples
  void setNumSamples(int n);
  ///save the position of current walkers to SampleStack
  //void saveEnsemble();
  ///save the position of current walkers
  void saveEnsemble(walker_iterator first, walker_iterator last);
  /// load a single sample from SampleStack
  void loadSample(ParticleSet::ParticlePos_t& Pos, size_t iw) const;
  /** load SampleStack data to current walkers
   */
  //void loadEnsemble(PropertySetType& PropertyList, WalkerList_t& WalkerList);
  //void loadEnsemble(const Walker_t& wcopy);
  /** load SampleStack from others
    */
  //void loadEnsemble(std::vector<MCWalkerConfiguration*>& others, bool doclean = true);
  /** dump Samples to a file
   * @param others MCWalkerConfigurations whose samples will be collected
   * @param out engine to write the samples to state_0/walkers
   * @param np number of processors
   * @return true with non-zero samples
   */
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

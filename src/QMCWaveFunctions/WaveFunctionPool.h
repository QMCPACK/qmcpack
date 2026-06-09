//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


/**@file WaveFunctionPool.h
 * @brief Declaration of WaveFunctionPool
 */
#ifndef QMCPLUSPLUS_WAVEFUNCTIONPOOL_H
#define QMCPLUSPLUS_WAVEFUNCTIONPOOL_H

#include "OhmmsData/OhmmsElementBase.h"
#include "Message/MPIObjectBase.h"
#include "QMCWaveFunctions/WaveFunctionFactory.h"
#include "Utilities/RuntimeOptions.h"
#include "Utilities/ObjectPool.h"
#include <map>
#include <string>

namespace qmcplusplus
{
class ParticleSetPool;
class ParticleSet;

extern template class ObjectPool<TrialWaveFunction>;

/** @ingroup qmcapp
 * @brief Manage a collection of TrialWaveFunction objects
 *
 * This object handles \<wavefunction\> elements and
 * functions as a builder class for TrialWaveFunction objects.
 */
class WaveFunctionPool : public MPIObjectBase, public ObjectPool<TrialWaveFunction>
{
public:
  using PoolType = typename ObjectPool<TrialWaveFunction>::Pool;

  WaveFunctionPool(const RuntimeOptions& runtime_options, ParticleSetPool& pset_pool, Communicate* c);
  WaveFunctionPool(const WaveFunctionPool&)            = delete;
  WaveFunctionPool& operator=(const WaveFunctionPool&) = delete;
  WaveFunctionPool(WaveFunctionPool&&)                 = default;
  WaveFunctionPool& operator=(WaveFunctionPool&&)      = delete;

  ~WaveFunctionPool();

  bool put(xmlNodePtr cur);

  /** look up wavefunction by name
   * @param pname wavefunction name to look up
   * if pname is empty and the pool contains one entry, return the only entry
   * if pname is not empty and not found in the pool, throw error
   */
  OptionalRef<TrialWaveFunction> getWaveFunction(const std::string& pname = "");

  /** return a xmlNode containing Jastrow
   * @param id name of the wave function
   *
   * If the wavefunction with id does not exist, return 0
   */
  xmlNodePtr getWaveFunctionNode(const std::string& id);

  /** get the Pool object
   */
  inline const PoolType& getPool() const { return myPool; }

private:
  /// @brief top-level runtime options from project data information
  const RuntimeOptions& runtime_options_;

  /** pointer to ParticleSetPool
   *
   * TrialWaveFunction needs to know which ParticleSet object
   * is used as an input object for the evaluations.
   */
  ParticleSetPool& ptcl_pool_;
};
} // namespace qmcplusplus
#endif

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
#include <map>
#include <string>

namespace qmcplusplus
{
class ParticleSetPool;
class ParticleSet;

/** @ingroup qmcapp
 * @brief Manage a collection of TrialWaveFunction objects
 *
 * This object handles \<wavefunction\> elements and
 * functions as a builder class for TrialWaveFunction objects.
 */
class WaveFunctionPool : public MPIObjectBase
{
public:
  using PoolType = std::map<std::string, const std::unique_ptr<TrialWaveFunction>>;

  WaveFunctionPool(ParticleSetPool& pset_pool, Communicate* c, const char* aname = "wavefunction");
  WaveFunctionPool(const WaveFunctionPool&)            = delete;
  WaveFunctionPool& operator=(const WaveFunctionPool&) = delete;
  WaveFunctionPool(WaveFunctionPool&&)                 = default;
  WaveFunctionPool& operator=(WaveFunctionPool&&)      = delete;

  ~WaveFunctionPool();

  bool put(xmlNodePtr cur);

  inline bool empty() const { return myPool.empty(); }

  TrialWaveFunction* getPrimary() { return primary_psi_; }

  TrialWaveFunction* getWaveFunction(const std::string& pname)
  {
    if (auto pit(myPool.find(pname)); pit == myPool.end())
    {
      if (myPool.empty())
        return nullptr;
      else
        return myPool.begin()->second.get();
    }
    else
      return pit->second.get();
  }

  /** return a xmlNode containing Jastrow
   * @param id name of the wave function
   *
   * If the wavefunction with id does not exist, return 0
   */
  xmlNodePtr getWaveFunctionNode(const std::string& id);

  /** get the Pool object
   */
  inline const PoolType& getPool() const { return myPool; }

  /** add a TrialWaveFunction* to myPool
   */
  void addFactory(std::unique_ptr<TrialWaveFunction> psi, bool primary);

private:
  /// pointer to the primary TrialWaveFunction
  TrialWaveFunction* primary_psi_;

  /// storage of WaveFunctionFactory
  PoolType myPool;

  /** pointer to ParticleSetPool
   *
   * TrialWaveFunction needs to know which ParticleSet object
   * is used as an input object for the evaluations.
   */
  ParticleSetPool& ptcl_pool_;
};
} // namespace qmcplusplus
#endif

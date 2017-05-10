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

  typedef std::map<std::string,WaveFunctionFactory*> PoolType;

  WaveFunctionPool(Communicate* c, const char* aname = "wavefunction");
  ~WaveFunctionPool();

  bool put(xmlNodePtr cur);

  inline bool empty() const
  {
    return myPool.empty();
  }

  TrialWaveFunction* getPrimary()
  {
    return primaryPsi;
  }

  TrialWaveFunction* getWaveFunction(const std::string& pname)
  {
    std::map<std::string,WaveFunctionFactory*>::iterator pit(myPool.find(pname));
    if(pit == myPool.end())
    {
      if(myPool.empty())
        return 0;
      else
        return (*(myPool.begin())).second->targetPsi;
    }
    else
      return (*pit).second->targetPsi;
  }

  WaveFunctionFactory* getWaveFunctionFactory(const std::string& pname)
  {
    std::map<std::string,WaveFunctionFactory*>::iterator pit(myPool.find(pname));
    if(pit == myPool.end())
    {
      if(myPool.empty())
        return 0;
      else
        return (*(myPool.begin())).second;
    }
    else
      return (*pit).second;
  }

  /** assign a pointer of ParticleSetPool
   */
  inline void
  setParticleSetPool(ParticleSetPool* pset)
  {
    ptclPool=pset;
  }

  /** return a xmlNode containing Jastrow
   * @param id name of the wave function
   *
   * If the wavefunction with id does not exist, return 0
   */
  xmlNodePtr getWaveFunctionNode(const std::string& id);

  /** get the Pool object
   */
  inline PoolType& getPool()
  {
    return myPool;
  }

  /** add a WaveFunctionFactory* to myPool
   */
  void addFactory(WaveFunctionFactory* psifac);

private:

  /// pointer to the primary TrialWaveFunction
  TrialWaveFunction* primaryPsi;

  /// storage of WaveFunctionFactory
  PoolType myPool;

  /** pointer to ParticleSetPool
   *
   * TrialWaveFunction needs to know which ParticleSet object
   * is used as an input object for the evaluations.
   */
  ParticleSetPool* ptclPool;

};
}
#endif

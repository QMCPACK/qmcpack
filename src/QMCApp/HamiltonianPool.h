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
    
    



/**@file HamiltonianPool.h
 * @brief Declaration of HamiltonianPool
 */
#ifndef QMCPLUSPLUS_QMCHAMILTONIANS_H
#define QMCPLUSPLUS_QMCHAMILTONIANS_H

#include "QMCHamiltonians/HamiltonianFactory.h"
#include "OhmmsData/OhmmsElementBase.h"
#include "Message/MPIObjectBase.h"
#include <map>

namespace qmcplusplus
{

class ParticleSet;
class MCWalkerConfiguration;
class ParticleSetPool;
class WaveFunctionPool;
class Libxml2Document;

/** @ingroup qmcapp
 * @brief Manage a collection of QMCHamiltonian objects
 *
 * This object handles \<hamiltonian\> elements and
 * functions as a builder class for QMCHamiltonian objects.
 */
class HamiltonianPool : public MPIObjectBase
{

public:

  typedef std::map<std::string,HamiltonianFactory*> PoolType;

  HamiltonianPool(Communicate* c, const char* aname = "hamiltonian");

  bool put(xmlNodePtr cur);
  bool get(std::ostream& os) const;
  void reset();

  inline bool empty() const
  {
    return myPool.empty();
  }

  /** return the pointer to the primary QMCHamiltonian
   *
   * The first QMCHamiltonian is assigned to the primaryH.
   * The last QMCHamiltonian with role="primary" will be the primaryH.
   */
  inline QMCHamiltonian* getPrimary()
  {
    return primaryH;
  }

  /** return the pointer to a QMCHamiltonian with the name
   * @param pname name of the QMCHamiltonian
   */
  inline QMCHamiltonian* getHamiltonian(const std::string& pname)
  {
    PoolType::iterator hit(myPool.find(pname));
    if(hit == myPool.end())
    {
      if(myPool.empty())
        return 0;
      else
        return (*(myPool.begin())).second->targetH;
    }
    else
      return (*hit).second->targetH;
  }

  void setDocument(Libxml2Document* doc)
  {
    curDoc=doc;
  }

  /** assign a pointer to a ParticleSetPool
   */
  inline void setParticleSetPool(ParticleSetPool* pset)
  {
    ptclPool=pset;
  }

  /** assign a pointer to a WaveFunctionPool
   */
  inline void setWaveFunctionPool(WaveFunctionPool* pset)
  {
    psiPool=pset;
  }

private:

  /** pointer to the primary QMCHamiltonian
   */
  QMCHamiltonian* primaryH;

  /** pointer to a current QMCHamiltonian to be built.
   */
  QMCHamiltonian* curH;

  /** pointer to ParticleSetPool
   *
   * QMCHamiltonian needs to know which ParticleSet object
   * is used as an input object for the evaluations.
   * Any number of ParticleSet can be used to describe
   * a QMCHamiltonian.
   */
  ParticleSetPool* ptclPool;

  /** pointer to WaveFunctionPool
   *
   * For those QMCHamiltonianBase that depends on TrialWaveFunction,
   * e.g., NonLocalPPotential.
   */
  WaveFunctionPool* psiPool;


  /** point to the working document */
  Libxml2Document* curDoc;

  /** storage for HamiltonianFactory */
  PoolType myPool;
};
}
#endif

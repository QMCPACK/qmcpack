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

struct Libxml2Document;

namespace qmcplusplus
{
class ParticleSet;
class MCWalkerConfiguration;
class ParticleSetPool;
class WaveFunctionPool;

/** @ingroup qmcapp
 * @brief Manage a collection of QMCHamiltonian objects
 *
 * This object handles \<hamiltonian\> elements and
 * functions as a builder class for QMCHamiltonian objects.
 */
class HamiltonianPool : public MPIObjectBase
{
public:
  using PoolType = std::map<std::string, HamiltonianFactory*>;

  HamiltonianPool(ParticleSetPool& pset_pool,
                  WaveFunctionPool& psi_pool,
                  Communicate* c,
                  const char* aname = "hamiltonian");
  HamiltonianPool(const HamiltonianPool&) = delete;
  HamiltonianPool& operator=(const HamiltonianPool&) = delete;
  HamiltonianPool(HamiltonianPool&&)                 = default;
  HamiltonianPool& operator=(HamiltonianPool&&) = delete;
  ~HamiltonianPool();

  bool put(xmlNodePtr cur);
  bool get(std::ostream& os) const;
  void reset();

  inline bool empty() const { return myPool.empty(); }

  /** return the pointer to the primary QMCHamiltonian
   *
   * The first QMCHamiltonian is assigned to the primaryH.
   * The last QMCHamiltonian with role="primary" will be the primaryH.
   */
  inline QMCHamiltonian* getPrimary() { return primaryH; }

  /** return the pointer to a QMCHamiltonian with the name
   * @param pname name of the QMCHamiltonian
   */
  inline QMCHamiltonian* getHamiltonian(const std::string& pname)
  {
    PoolType::iterator hit(myPool.find(pname));
    if (hit == myPool.end())
    {
      if (myPool.empty())
        return nullptr;
      else
        return (*(myPool.begin())).second->getH();
    }
    else
      return (*hit).second->getH();
  }

  void setDocument(Libxml2Document* doc) { curDoc = doc; }

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
  ParticleSetPool& ptcl_pool_;

  /** pointer to WaveFunctionPool
   *
   * For those OperatorBase that depends on TrialWaveFunction,
   * e.g., NonLocalPPotential.
   */
  WaveFunctionPool& psi_pool_;


  /** point to the working document */
  Libxml2Document* curDoc;

  /** storage for HamiltonianFactory */
  PoolType myPool;
};
} // namespace qmcplusplus
#endif

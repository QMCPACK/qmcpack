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

#include "QMCHamiltonian.h"
#include "OhmmsData/OhmmsElementBase.h"
#include "Message/MPIObjectBase.h"
#include "Utilities/ObjectPool.h"

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
class HamiltonianPool : public MPIObjectBase, public ObjectPool<QMCHamiltonian>
{
public:
  HamiltonianPool(ParticleSetPool& pset_pool,
                  WaveFunctionPool& psi_pool,
                  Communicate* c,
                  const char* aname = "hamiltonian");
  HamiltonianPool(const HamiltonianPool&)            = delete;
  HamiltonianPool& operator=(const HamiltonianPool&) = delete;
  HamiltonianPool(HamiltonianPool&&)                 = default;
  HamiltonianPool& operator=(HamiltonianPool&&)      = delete;
  ~HamiltonianPool();

  bool put(xmlNodePtr cur);
  bool get(std::ostream& os) const;

  /** look up hamiltonian by name
   * @param pname hamiltonian name to look up
   * if pname is empty and the pool contains one entry, return the only entry
   * if pname is not empty and not found in the pool, throw error
   */
  OptionalRef<QMCHamiltonian> getHamiltonian(const std::string& pname = "");

private:
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
};
} // namespace qmcplusplus
#endif

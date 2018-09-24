//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2018 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_QMCHAMILTONIAN_POOL_INTERFACE_H
#define QMCPLUSPLUS_QMCHAMILTONIAN_POOL_INTERFACE_H
#include "QMCHamiltonians/HamiltonianFactory.h"
#include "OhmmsData/OhmmsElementBase.h"
#include "Message/MPIObjectBase.h"
#include <map>
#include "QMCApp/WaveFunctionPool.h"
#include "QMCWaveFunctions/Batching.h"
class Libxml2Document;

namespace qmcplusplus
{

class ParticleSet;
class MCWalkerConfiguration;
class ParticleSetPool;
  //class WaveFunctionPool;

/** @ingroup qmcapp
 * @brief Manage a collection of QMCHamiltonian objects
 *
 * This object handles \<hamiltonian\> elements and
 * functions as a builder class for QMCHamiltonian objects.
 */
class HamiltonianPoolInterface
{

public:

  virtual bool put(xmlNodePtr cur) = 0;
  virtual bool get(std::ostream& os) const = 0;
  virtual void reset() = 0;

  virtual bool empty() const = 0;

  /** return the pointer to the primary QMCHamiltonian
   *
   * The first QMCHamiltonian is assigned to the primaryH.
   * The last QMCHamiltonian with role="primary" will be the primaryH.
   */
  virtual QMCHamiltonian* getPrimary() = 0;

  /** return the pointer to a QMCHamiltonian with the name
   * @param pname name of the QMCHamiltonian
   */
  virtual QMCHamiltonian* getHamiltonian(const std::string& pname) = 0;

  virtual void setDocument(Libxml2Document* doc) = 0;

  /** assign a pointer to a ParticleSetPool
   */
  virtual void setParticleSetPool(ParticleSetPool* pset) = 0;

  /** assign a pointer to a WaveFunctionPool
   */
  virtual void setWaveFunctionPool(WaveFunctionPool* pset) = 0;

};
}
#endif

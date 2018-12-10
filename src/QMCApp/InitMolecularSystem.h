//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



/**@file InitMolecularSystem.h
 * @brief Declaration of InitMolecularSystem
 */
#ifndef QMCPLUSPLUS_INITMOLECULARSYSTEM_H
#define QMCPLUSPLUS_INITMOLECULARSYSTEM_H

#include "OhmmsData/OhmmsElementBase.h"
#include <map>

namespace qmcplusplus
{

class ParticleSet;
class ParticleSetPool;

/* Engine to initialize the initial electronic structure for a molecular system
 */
class InitMolecularSystem : public OhmmsElementBase
{

public:

  InitMolecularSystem(ParticleSetPool* pset, const char* aname = "mosystem");

  bool get(std::ostream& os) const;
  bool put(std::istream& is);
  bool put(xmlNodePtr cur);
  void reset();

  /** initialize els for an atom
   */
  void initAtom(ParticleSet* ions, ParticleSet* els);
  /** initialize els position for a molecule
   *
   * Use the valence of each ionic species on a sphere
   */
  void initMolecule(ParticleSet* ions, ParticleSet* els);
  /** initialize els for the systems with a mixed boundary
   *
   * Use the bound of the ionic systems and uniform random positions within a reduced box
   */
  void initWithVolume(ParticleSet* ions, ParticleSet* els);

private:

  /** pointer to ParticleSetPool
   *
   * QMCHamiltonian needs to know which ParticleSet object
   * is used as an input object for the evaluations.
   * Any number of ParticleSet can be used to describe
   * a QMCHamiltonian.
   */
  ParticleSetPool* ptclPool;

};
}
#endif

//////////////////////////////////////////////////////////////////
// (c) Copyright 2003- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//   Department of Physics, Ohio State University
//   Ohio Supercomputer Center
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
/**@file InitMolecularSystem.h
 * @brief Declaration of InitMolecularSystem
 */
#ifndef QMCPLUSPLUS_INITMOLECULARSYSTEM_H
#define QMCPLUSPLUS_INITMOLECULARSYSTEM_H

#include "OhmmsData/OhmmsElementBase.h"
#include <map>

namespace qmcplusplus {

  class ParticleSet;
  class ParticleSetPool;

  /* Engine to initialize the initial electronic structure for a molecular system
   */
  class InitMolecularSystem : public OhmmsElementBase {

  public:

    InitMolecularSystem(ParticleSetPool* pset, const char* aname = "mosystem");

    bool get(std::ostream& os) const;
    bool put(std::istream& is);
    bool put(xmlNodePtr cur);
    void reset();

    void initAtom(ParticleSet* ions, ParticleSet* els);
    void initMolecule(ParticleSet* ions, ParticleSet* els);

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
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

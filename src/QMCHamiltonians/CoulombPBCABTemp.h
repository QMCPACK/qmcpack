//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim and Kris Delaney
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
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
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef QMCPLUSPLUS_COULOMBPBCAB_TEMP_H
#define QMCPLUSPLUS_COULOMBPBCAB_TEMP_H
#include "QMCHamiltonians/QMCHamiltonianBase.h"
#include "LongRange/LRCoulombSingleton.h"

namespace qmcplusplus {

  /** @ingroup hamiltonian
   *\brief Calculates the AA Coulomb potential using PBCs
   *
   * Functionally identical to CoulombPBCAB but uses a templated version of
   * LRHandler.
   */
  struct CoulombPBCABTemp: public QMCHamiltonianBase {

    typedef LRCoulombSingleton::LRHandlerType LRHandlerType;
    ParticleSet* PtclA;
    ParticleSet* PtclB;
    LRHandlerType* AB;
    DistanceTableData* d_ab;

    bool FirstTime;
    int NumSpeciesA;
    int NumSpeciesB;
    int ChargeAttribIndxA;
    int ChargeAttribIndxB;
    int MemberAttribIndxA;
    int MemberAttribIndxB;
    int NptclA;
    int NptclB;
    RealType myConst;
    RealType myRcut;

    vector<RealType> Zat,Zspec; 
    vector<RealType> Qat,Qspec; 
    vector<int> NofSpeciesA;
    vector<int> NofSpeciesB;

    //This is set to true if the K_c of structure-factors are different
    bool kcdifferent; 
    RealType minkc;

    CoulombPBCABTemp(ParticleSet& ions, ParticleSet& elns);

    /// copy constructor
    CoulombPBCABTemp(const CoulombPBCABTemp& c);
    
    ~CoulombPBCABTemp();

    void resetTargetParticleSet(ParticleSet& P);

    Return_t evaluate(ParticleSet& P);

    inline Return_t evaluate(ParticleSet& P, vector<NonLocalData>& Txy) {
      return evaluate(P);
    }

    /** Do nothing */
    bool put(xmlNodePtr cur) {
      return true;
    }

    bool get(std::ostream& os) const {
      os << "CoulombPBCAB potential: " << PtclA->getName() << "-" << PtclB->getName();
      return true;
    }

    void initBreakup();

    Return_t evalSR();
    Return_t evalLR();
    Return_t evalConsts();

  };

}
#endif

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/


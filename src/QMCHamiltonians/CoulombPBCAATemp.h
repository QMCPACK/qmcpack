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
#ifndef QMCPLUSPLUS_COULOMBPBCAA_TEMP_H
#define QMCPLUSPLUS_COULOMBPBCAA_TEMP_H
#include "QMCHamiltonians/QMCHamiltonianBase.h"
#include "LongRange/LRCoulombSingleton.h"

namespace qmcplusplus {

  /** @ingroup hamiltonian
   *\brief Calculates the AA Coulomb potential using PBCs
   *
   * Functionally identical to CoulombPBCAA but uses a templated version of
   * LRHandler.
   */
  struct CoulombPBCAATemp: public QMCHamiltonianBase {

    typedef LRCoulombSingleton::LRHandlerType LRHandlerType;
    typedef LRCoulombSingleton::GridType       GridType;
    typedef LRCoulombSingleton::RadFunctorType RadFunctorType;
    ParticleSet* PtclRef;
    DistanceTableData* d_aa;
    LRHandlerType* AA;
    GridType* myGrid;
    RadFunctorType* rVs;

    bool is_active;
    bool FirstTime;
    int NumSpecies;
    int ChargeAttribIndx;
    int MemberAttribIndx;
    int NumCenters;
    RealType myConst;
    RealType myRcut;
    vector<RealType> Zat,Zspec; 
    vector<int> NofSpecies;

    Matrix<RealType> SR2;
    Vector<RealType> dSR;
    Vector<ComplexType> del_eikr;
    CoulombPBCAATemp(ParticleSet& ref, bool active);

    ///// copy constructor
    //CoulombPBCAATemp(const CoulombPBCAATemp& c);
    
    ~CoulombPBCAATemp();

    void resetTargetParticleSet(ParticleSet& P);

    Return_t evaluate(ParticleSet& P);

    inline Return_t evaluate(ParticleSet& P, vector<NonLocalData>& Txy) {
      return evaluate(P);
    }
    Return_t registerData(ParticleSet& P, BufferType& buffer);
    Return_t updateBuffer(ParticleSet& P, BufferType& buffer);
    void copyFromBuffer(ParticleSet& P, BufferType& buf);
    void copyToBuffer(ParticleSet& P, BufferType& buf);
    Return_t evaluatePbyP(ParticleSet& P, int iat);
    void acceptMove(int iat);

    /** Do nothing */
    bool put(xmlNodePtr cur) {
      return true;
    }

    bool get(std::ostream& os) const {
      os << "CoulombPBCAA potential: " << PtclRef->getName();
      return true;
    }

    QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi);

    void initBreakup();

    Return_t evalSR();
    Return_t evalLR();
    Return_t evalConsts();
    Return_t evaluateForPbyP(ParticleSet& P);

  };

}
#endif

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/


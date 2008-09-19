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
#include "Numerics/OneDimGridBase.h"
#include "Numerics/OneDimGridFunctor.h"
#include "Numerics/OneDimCubicSpline.h"

namespace qmcplusplus {

  /** @ingroup hamiltonian
   *\brief Calculates the AA Coulomb potential using PBCs
   *
   * Functionally identical to CoulombPBCAB but uses a templated version of
   * LRHandler.
   */
  struct CoulombPBCABTemp: public QMCHamiltonianBase {

    typedef LRCoulombSingleton::LRHandlerType LRHandlerType;
    typedef LRCoulombSingleton::GridType GridType;
    typedef LRCoulombSingleton::RadFunctorType RadFunctorType;

    ///source particle set
    ParticleSet& PtclA;
    ///long-range Handler
    LRHandlerType* AB;
    ///locator of the distance table 
    int myTableIndex;
    ///number of species of A particle set
    int NumSpeciesA;
    ///number of species of B particle set
    int NumSpeciesB;
    ///number of particles of A
    int NptclA;
    ///number of particles of B
    int NptclB;
    ///const energy after breakup
    RealType myConst;
    ///cutoff radius of the short-range part
    RealType myRcut;
    ///radial grid
    GridType* myGrid;
    ///Always mave a radial functor for the bare coulomb
    RadFunctorType* V0;

    ///number of particles per species of A
    vector<int> NofSpeciesA;
    ///number of particles per species of B
    vector<int> NofSpeciesB;
    ///Zat[iat] charge for the iat-th particle of A
    vector<RealType> Zat;
    ///Qat[iat] charge for the iat-th particle of B
    vector<RealType> Qat; 
    ///Zspec[spec] charge for the spec-th species of A
    vector<RealType> Zspec; 
    ///Qspec[spec] charge for the spec-th species of B
    vector<RealType> Qspec; 
    ///Short-range potential for each ion
    vector<RadFunctorType*> Vat;
    ///Short-range potential for each species
    vector<RadFunctorType*> Vspec;
    /*@{
     * @brief temporary data for pbyp evaluation
     */
    ///short-range part for the moved particle
    RealType SRtmp;
    ///long-range part for the moved particle
    RealType LRtmp;
    ///short-range per particle
    Vector<RealType> SRpart;
    ///long-range per particle
    Vector<RealType> LRpart;
    /*@}*/

    //This is set to true if the K_c of structure-factors are different
    bool kcdifferent; 
    RealType minkc;

    CoulombPBCABTemp(ParticleSet& ions, ParticleSet& elns);

    ///// copy constructor
    //CoulombPBCABTemp(const CoulombPBCABTemp& c);
    
    ~CoulombPBCABTemp();

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
      //os << "CoulombPBCAB potential: " << PtclA->getName() << "-" << PtclB->getName();
      return true;
    }

    QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi);

    void initBreakup(ParticleSet& P);

    Return_t evalSR(ParticleSet& P);
    Return_t evalLR(ParticleSet& P);
    Return_t evalConsts();
    Return_t evaluateForPyP(ParticleSet& P);
    void add(int groupID, RadFunctorType* ppot);
  };

}
#endif

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/


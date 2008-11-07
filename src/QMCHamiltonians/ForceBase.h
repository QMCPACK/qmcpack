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
#ifndef QMCPLUSPLUS_FORCE_BASE_HAMILTONIAN_H
#define QMCPLUSPLUS_FORCE_BASE_HAMILTONIAN_H
#include "QMCHamiltonians/QMCHamiltonianBase.h"

namespace qmcplusplus {

  struct ForceBase
  {
    int FirstForceIndex;
    int myTableIndex;
    int Nnuc;
    int Nel;
    int tries;
    bool FirstTime;
    ParticleSet& Ions;
    ParticleSet::ParticlePos_t forces;
    ParticleSet::ParticlePos_t forces_IonIon;
    vector<QMCTraits::RealType> Zat; 
    vector<QMCTraits::RealType> Qat;
    string prefix;
    string pairName;

    ForceBase(ParticleSet& ions, ParticleSet& elns);
    virtual ~ForceBase(){}

    void addObservablesF(QMCTraits::PropertySetType& plist);
    void setObservablesF(QMCTraits::PropertySetType& plist);
    void setParticleSetF(QMCTraits::PropertySetType& plist, int offset);
  };

  struct BareForce: public QMCHamiltonianBase, public ForceBase
  {
    BareForce(ParticleSet& ions, ParticleSet& elns);
    void resetTargetParticleSet(ParticleSet& P);

    Return_t evaluate(ParticleSet& P);

    inline Return_t evaluate(ParticleSet& P, vector<NonLocalData>& Txy) 
    {
      return evaluate(P);
    }

    void addObservables(PropertySetType& plist)
    {
      addObservablesF(plist);
      myIndex=FirstForceIndex;
    }
    void setObservables(PropertySetType& plist)
    {
      setObservablesF(plist);
    }

    void setParticlePropertyList(PropertySetType& plist, int offset)
    {
      setParticleSetF(plist, offset);
    }

    /** Do nothing */
    bool put(xmlNodePtr cur) {
      return true;
    }

    bool get(std::ostream& os) const {
      os << "Force Base Hamiltonian: " << pairName;
      return true;
    }

    QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi);
  };

}
#endif

/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 2945 $   $Date: 2008-08-05 10:21:33 -0500 (Tue, 05 Aug 2008) $
 * $Id: ForceBase.h 2945 2008-08-05 15:21:33Z jnkim $ 
 ***************************************************************************/

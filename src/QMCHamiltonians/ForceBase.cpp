//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim and Kris Delaney
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
#include "QMCHamiltonians/ForceBase.h"
#include "Particle/DistanceTable.h"
#include "Particle/DistanceTableData.h"
#include "Message/Communicate.h"
#include "Utilities/ProgressReportEngine.h"

namespace qmcplusplus {

  ForceBase::ForceBase(ParticleSet& ions, ParticleSet& elns)
    : FirstForceIndex(-1),tries(0), Ions(ions)
    {
      ReportEngine PRE("ForceBase","ForceBase");
      myTableIndex=elns.addTable(ions);
      FirstTime = true;

      Nnuc = ions.getTotalNum();
      Nel = elns.getTotalNum();
      SpeciesSet& tspeciesA(ions.getSpeciesSet());
      SpeciesSet& tspeciesB(elns.getSpeciesSet());
      int ChargeAttribIndxA = tspeciesA.addAttribute("charge");
      int MemberAttribIndxA = tspeciesA.addAttribute("membersize");
      int ChargeAttribIndxB = tspeciesB.addAttribute("charge");
      int MemberAttribIndxB = tspeciesB.addAttribute("membersize");
      int NumSpeciesA = tspeciesA.TotalNum;
      int NumSpeciesB = tspeciesB.TotalNum;
      //Store information about charges and number of each species
      Zat.resize(Nnuc);
      Qat.resize(Nel); 
      for(int iat=0; iat<Nnuc; iat++)
        Zat[iat] = tspeciesA(ChargeAttribIndxA,ions.GroupID[iat]);
      for(int iat=0; iat<Nel; iat++)
        Qat[iat] = tspeciesB(ChargeAttribIndxB,elns.GroupID[iat]);

      pairName=elns.getName()+"-"+ions.getName();

      forces.resize(Nnuc);
      forces = 0.0;
      forces_IonIon.resize(Nnuc);
      forces_IonIon = 0.0;
    }

  void ForceBase::addObservablesF(QMCTraits::PropertySetType& plist) {
    if(FirstForceIndex<0) 
      FirstForceIndex=plist.size();
    //    cerr << "FirstForceIndex = " << FirstForceIndex << endl;
    for(int iat=0; iat<Nnuc; iat++) {
      for(int x=0; x<OHMMS_DIM; x++) {
        ostringstream obsName;
        obsName << prefix << "_" << iat << "_" << x;
        plist.add(obsName.str());
      }
    }
  }

  void ForceBase::registerObservablesF(vector<observable_helper*>& h5list
      , hid_t gid) const
  {
    vector<int> ndim(2);
    ndim[0]=Nnuc;
    ndim[1]=OHMMS_DIM;
    observable_helper* h5o=new observable_helper(prefix);
    h5o->set_dimensions(ndim,FirstForceIndex);
    h5o->open(gid);
    h5list.push_back(h5o);
  }

  void ForceBase::setObservablesF(QMCTraits::PropertySetType& plist) {
    // constant ion-ion contribution
    if(FirstTime) {
      FirstTime = false;
      forces_IonIon = 0.0;
      DistanceTableData* d_aa=DistanceTable::add(Ions);
      for(int iat=0; iat<Nnuc; iat++) {
        for(int nn=d_aa->M[iat], jat=1; nn<d_aa->M[iat+1]; nn++,jat++) {
          int jid = d_aa->J[nn];
          QMCTraits::RealType rinv=d_aa->rinv(nn);
          QMCTraits::RealType r3zz=Zat[jid]*Zat[iat]*rinv*rinv*rinv;
          forces_IonIon[iat] -= r3zz*d_aa->dr(nn);
          forces_IonIon[jid] += r3zz*d_aa->dr(nn);
        }
      }
    }

    int index = FirstForceIndex;
    for(int iat=0; iat<Nnuc; iat++) {
      //      cerr << "Forces[iat] = " << forces[iat] << " index = " << index << endl;
      for(int x=0; x<OHMMS_DIM; x++) {
        plist[index] = forces[iat][x];// + forces_IonIon[iat][x];
        index++;
      }
    }
  }

  void ForceBase::setParticleSetF(QMCTraits::PropertySetType& plist, int offset)
  {
    int index = FirstForceIndex + offset;
    for(int iat=0; iat<Nnuc; iat++) {
      for(int x=0; x<OHMMS_DIM; x++) {
        plist[index] = forces[iat][x];// + forces_IonIon[iat][x];
        index++;
      }
    }
  }

  BareForce::BareForce(ParticleSet& ions, ParticleSet& elns): ForceBase(ions,elns)
  {
    myName = "HF_Force_Base";
    prefix="HFBase";
  }

  void BareForce::resetTargetParticleSet(ParticleSet& P) { }

  QMCHamiltonianBase* BareForce::makeClone(ParticleSet& qp, TrialWaveFunction& psi)
  {
    return new BareForce(*this);
  }


  BareForce::Return_t 
    BareForce::evaluate(ParticleSet& P) {
      forces = 0.0;
      const DistanceTableData* d_ab=P.DistTables[myTableIndex];

      //Loop over distinct eln-ion pairs
      for(int iat=0; iat<Nnuc; iat++)
      {
        for(int nn=d_ab->M[iat], jat=0; nn<d_ab->M[iat+1]; nn++,jat++) {
          RealType rinv=d_ab->rinv(nn);
          RealType r3zz=Qat[jat]*Zat[iat]*rinv*rinv*rinv;
          forces[iat] -= r3zz*d_ab->dr(nn);
        }

      }
      tries++;
      return 0.0;
    }

}

/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 3015 $   $Date: 2008-08-18 16:08:06 -0500 (Mon, 18 Aug 2008) $
 * $Id: ForceBase.cpp 3015 2008-08-18 21:08:06Z jnkim $ 
 ***************************************************************************/


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
#include "Estimators/HFDHE2PolymerEstimator.h"
#include "QMCHamiltonians/QMCHamiltonian.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "ParticleBase/ParticleAttribOps.h"
#include "Message/CommOperators.h"
#include "QMCDrivers/DriftOperators.h"
#include "QMCDrivers/MultiChain.h"

namespace qmcplusplus {

  /** constructor
   * @param h QMCHamiltonian to define the components
   * @param hcopy number of copies of QMCHamiltonians
   */
  HFDHE2PolymerEstimator::HFDHE2PolymerEstimator(QMCHamiltonian& h, int hcopy, MultiChain* polymer)
  {
    Reptile=polymer;
    NumCopies=hcopy;
//     Hpointer = &h;

    SizeOfHamiltonians = h.sizeOfObservables();
    FirstHamiltonian = h.startIndex();
//     cout<<"size of Hamiltonian "<<SizeOfHamiltonians<<" First one "<<FirstHamiltonian<<endl;

    elocal_name.push_back("LocalEnergy");
    for(int i=0; i < SizeOfHamiltonians; i++)
    {
      elocal_name.push_back(h.getObservableName(i));
      //elocal_name.push_back(h.getName(i));
    };
    
    elocal_name.push_back("SumPot");
    elocal_name.push_back("ElSumPot");
    elocal_name.push_back("Virial");
    elocal_name.push_back("MaxAge");

    scalars.resize(SizeOfHamiltonians+5);
    scalars_saved=scalars;
    pNorm=0.0;
    
    Pindex = h.getHamiltonian("HePress")->myIndex;
    HDFHE2index = h.getHamiltonian("HFDHE2")->myIndex;
  }


  HFDHE2PolymerEstimator::HFDHE2PolymerEstimator(const HFDHE2PolymerEstimator& mest): 
       PolymerEstimator(mest)
//       Reptile(mest.Reptile)
  {
    Reptile=mest.Reptile;
    NumCopies=mest.NumCopies;
    //d_data.resize(mest.d_data.size());
  }

  ScalarEstimatorBase* HFDHE2PolymerEstimator::clone()
  {
    return new HFDHE2PolymerEstimator(*this);
  }

  /**  add the local energy, variance and all the Hamiltonian components to the scalar record container
   *@param record storage of scalar records (name,value)
   */
  void 
  HFDHE2PolymerEstimator::add2Record(RecordNamedProperty<RealType>& record) {
    FirstIndex = record.add(elocal_name[0].c_str());
    for(int i=1; i<elocal_name.size(); i++) record.add(elocal_name[i].c_str());
    LastIndex=FirstIndex + elocal_name.size();
    clear();
    char aname[32];
    for(int i=0; i<NumCopies-1; i++) {
      for(int j=i+1; j<NumCopies; j++) {
        sprintf(aname,"DiffS%iS%i",i,j); 
        int dummy=record.add(aname);
      }
    }
  };

  void HFDHE2PolymerEstimator::accumulate(WalkerIterator first, WalkerIterator last, RealType wgt)
  {
    for(int i=0; i<NumCopies; i++) 
    {
      RealType uw(Reptile->UmbrellaWeight[i]);
      //Adding all parts from the Hamiltonian
      RealType* restrict HeadProp(Reptile->front()->getPropertyBase(i));
      RealType* restrict TailProp(Reptile->back()->getPropertyBase(i));
      RealType eloc = 0.5*( HeadProp[LOCALENERGY] + TailProp[LOCALENERGY]);
      scalars[0](eloc,uw);
      for(int obsi=0;obsi<SizeOfHamiltonians ;obsi++){
        scalars[obsi+1]( 0.5*( HeadProp[obsi+FirstHamiltonian] + TailProp[obsi+FirstHamiltonian]) , uw);
      };
      //Center Pressure
      RealType* restrict CenProp(Reptile->center()->getPropertyBase(i));
      scalars[SizeOfHamiltonians+3]( (2.0*eloc-2.0*CenProp[LOCALPOTENTIAL])*pNorm + CenProp[Pindex+1+FirstHamiltonian] ,uw);
      
      int Rage(Reptile->Age);
      int Bage=Rage;

      RealType tmpV=0.0;
      RealType tmpP=0.0;
      KEconst = (*(Reptile->begin()))->Drift.size()  * 1.5 * OneOverTau;
      for( MultiChain::iterator Bit = Reptile->begin();Bit != (Reptile->end());Bit++){
        tmpP+= (*Bit)->getPropertyBase(i)[Pindex+1+FirstHamiltonian];
        tmpV+=(*Bit)->getPropertyBase(i)[LOCALPOTENTIAL];
        Bage=min((*Bit)->stepmade,Bage);
      };

      tmpP-=0.5*Reptile->back()->getPropertyBase(i)[Pindex+1+FirstHamiltonian];
      tmpP-=0.5*Reptile->front()->getPropertyBase(i)[Pindex+1+FirstHamiltonian];
      tmpV-=0.5*Reptile->back()->getPropertyBase(i)[LOCALPOTENTIAL];
      tmpV-=0.5*Reptile->front()->getPropertyBase(i)[LOCALPOTENTIAL];
      tmpV *= -2.0*pNorm*Tau;
      tmpP *= Tau;

      scalars[SizeOfHamiltonians+1](tmpV+tmpP,uw);
      scalars[SizeOfHamiltonians+2](eloc*(tmpV+tmpP),uw);
      scalars[SizeOfHamiltonians+4](Rage-Bage,1.0);

    }
  }

  void 
  HFDHE2PolymerEstimator::evaluateDiff() 
  {

  }

}
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1926 $   $Date: 2007-04-20 12:30:26 -0500 (Fri, 20 Apr 2007) $
 * $Id: HFDHE2PolymerEstimator.cpp 1926 2007-04-20 17:30:26Z jnkim $
 ***************************************************************************/

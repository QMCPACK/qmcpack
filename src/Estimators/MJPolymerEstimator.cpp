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
#include "Estimators/MJPolymerEstimator.h"
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
  MJPolymerEstimator::MJPolymerEstimator(QMCHamiltonian& h, int hcopy, MultiChain* polymer)
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
//     elocal_name.push_back("ThermE");
    elocal_name.push_back("centerV");
    elocal_name.push_back("EcorrFun");
    elocal_name.push_back("Ehead");
    elocal_name.push_back("Etail");

    scalars.resize(SizeOfHamiltonians+9);
    scalars_saved=scalars;
    pNorm=0.0;
  };


  MJPolymerEstimator::MJPolymerEstimator(const MJPolymerEstimator& mest): 
       PolymerEstimator(mest)
//       Reptile(mest.Reptile)
  {
    Reptile=mest.Reptile;
    NumCopies=mest.NumCopies;
    //d_data.resize(mest.d_data.size());
  }

  ScalarEstimatorBase* MJPolymerEstimator::clone()
  {
    return new MJPolymerEstimator(*this);
  }

  /**  add the local energy, variance and all the Hamiltonian components to the scalar record container
   *@param record storage of scalar records (name,value)
   */
  void 
  MJPolymerEstimator::add2Record(RecordNamedProperty<RealType>& record) {
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

  void MJPolymerEstimator::accumulate(WalkerIterator first, WalkerIterator last, RealType wgt)
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
      scalars[SizeOfHamiltonians+3]( (2.0*eloc-CenProp[LOCALPOTENTIAL])*pNorm ,uw);
      scalars[SizeOfHamiltonians+5]( CenProp[LOCALPOTENTIAL] ,uw);
      scalars[SizeOfHamiltonians+6]( HeadProp[LOCALENERGY] * TailProp[LOCALENERGY] ,uw);
      scalars[SizeOfHamiltonians+7]( HeadProp[LOCALENERGY] ,uw);
      scalars[SizeOfHamiltonians+8]( TailProp[LOCALENERGY] ,uw);
      
      int Rage(Reptile->Age);
      int Bage=Rage;

//       RealType tmpE=0.0;
      RealType tmpV=0.0;
//       RealType tmpF=0.0;
//       RealType localE=0.0;
//       KEconst = (*(Reptile->begin()))->Drift.size()  * 1.5 * OneOverTau;
      for( MultiChain::iterator Bit = Reptile->begin();Bit != (Reptile->end());Bit++){
//         localE += 0.5*( (*Bit)->deltaRSquared[0] + (*Bit)->deltaRSquared[1]);
//         tmpF+= 0.5*(*Bit)->deltaRSquared[2];
//         tmpE+=(*Bit)->getPropertyBase(i)[LOCALENERGY];
        tmpV+=(*Bit)->getPropertyBase(i)[LOCALPOTENTIAL];
        Bage=min((*Bit)->stepmade,Bage);
      };

//       tmpF-=0.25*Reptile->back()->deltaRSquared[2];
//       tmpF-=0.25*Reptile->front()->deltaRSquared[2];
//       tmpE-=0.5*Reptile->back()->getPropertyBase(i)[LOCALENERGY];
//       tmpE-=0.5*Reptile->front()->getPropertyBase(i)[LOCALENERGY];
      tmpV-=0.5*Reptile->back()->getPropertyBase(i)[LOCALPOTENTIAL];
      tmpV-=0.5*Reptile->front()->getPropertyBase(i)[LOCALPOTENTIAL];
      tmpV *= -pNorm*Tau;

//       localE *= -0.5*OneOverTau*OneOverTau;
//       localE += KEconst* (Reptile->Last) + tmpF ;
//       localE += tmpE;

      scalars[SizeOfHamiltonians+1](tmpV,uw);
      scalars[SizeOfHamiltonians+2](eloc*tmpV,uw);

      scalars[SizeOfHamiltonians+4](Rage-Bage,1.0);

      ///This is the center bead energy using PIMC stuff.
//       localE *= 1.0/(Reptile->Last);
//       scalars[SizeOfHamiltonians+5]( localE ,uw);
    }
  }

  void 
  MJPolymerEstimator::evaluateDiff() 
  {

  }

}
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1926 $   $Date: 2007-04-20 12:30:26 -0500 (Fri, 20 Apr 2007) $
 * $Id: MJPolymerEstimator.cpp 1926 2007-04-20 17:30:26Z jnkim $
 ***************************************************************************/

//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
    
    
    
    
    
    
#include "Estimators/MJPolymerEstimator.h"
#include "QMCHamiltonians/QMCHamiltonian.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "ParticleBase/ParticleAttribOps.h"
#include "Message/CommOperators.h"
#include "QMCDrivers/DriftOperators.h"
#include "QMCDrivers/MultiChain.h"

namespace qmcplusplus
{

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
//     std::cout <<"size of Hamiltonian "<<SizeOfHamiltonians<<" First one "<<FirstHamiltonian<< std::endl;
  elocal_name.push_back("LocalEnergy");
  for(int i=0; i < SizeOfHamiltonians; i++)
  {
    elocal_name.push_back(h.getObservableName(i));
    //elocal_name.push_back(h.getName(i));
  };
  elocal_name.push_back("SumPot");
  elocal_name.push_back("ElSumPot");
  elocal_name.push_back("CenterTruncSumPot");
  elocal_name.push_back("Null");
  elocal_name.push_back("CenterTruncElSumPot");
  elocal_name.push_back("Virial");
  elocal_name.push_back("MaxAge");
  elocal_name.push_back("MaxTouch");
  elocal_name.push_back("centerV");
  elocal_name.push_back("EcorrFun");
  elocal_name.push_back("Ehead");
  elocal_name.push_back("Etail");
  elocal_name.push_back("Edecorr");
  elocal_name.push_back("Vdecorr");
//     elocal_name.push_back("RMC_HFCep_1_0");
//     elocal_name.push_back("RMC_HFCep_1_1");
// for(int i=0; i < elocal_name.size(); i++)
//   app_log()<<elocal_name[i]<<" ";
// app_log()<< std::endl;
  scalars.resize(SizeOfHamiltonians+14);
  scalars_saved=scalars;
  pNorm=0.0;
//     Findex = h.getObservable("HFCep_1_0");
//     app_log()<<"Force Index "<<Findex<< std::endl;
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
MJPolymerEstimator::add2Record(RecordNamedProperty<RealType>& record)
{
  FirstIndex = record.add(elocal_name[0].c_str());
  for(int i=1; i<elocal_name.size(); i++)
    record.add(elocal_name[i].c_str());
  LastIndex=FirstIndex + elocal_name.size();
  clear();
  char aname[32];
  for(int i=0; i<NumCopies-1; i++)
  {
    for(int j=i+1; j<NumCopies; j++)
    {
      sprintf(aname,"DiffS%iS%i",i,j);
      int dummy=record.add(aname);
    }
  }
};

void MJPolymerEstimator::accumulate(const MCWalkerConfiguration& W
                                    , WalkerIterator first, WalkerIterator last, RealType wgt)
{
  for(int i=0; i<NumCopies; i++)
  {
    RealType uw(Reptile->UmbrellaWeight[i]);
    //Adding all parts from the Hamiltonian
    RealType* restrict HeadProp(Reptile->front()->getPropertyBase(i));
    RealType* restrict TailProp(Reptile->back()->getPropertyBase(i));
    RealType eloc = 0.5*( HeadProp[LOCALENERGY] + TailProp[LOCALENERGY]);
    scalars[0](eloc,uw);
    for(int obsi=0; obsi<SizeOfHamiltonians ; obsi++)
    {
      scalars[obsi+1]( 0.5*( HeadProp[obsi+FirstHamiltonian] + TailProp[obsi+FirstHamiltonian]) , uw);
    };
    //Center Pressure
    RealType* restrict CenProp(Reptile->center()->getPropertyBase(i));
    scalars[SizeOfHamiltonians+6]( (2.0*eloc-CenProp[LOCALPOTENTIAL])*pNorm ,uw);
    scalars[SizeOfHamiltonians+9]( CenProp[LOCALPOTENTIAL] ,uw);
    scalars[SizeOfHamiltonians+10]( HeadProp[LOCALENERGY] * TailProp[LOCALENERGY] ,uw);
    RealType energy_head = HeadProp[LOCALENERGY];
    RealType energy_tail = TailProp[LOCALENERGY];
    scalars[SizeOfHamiltonians+11]( HeadProp[LOCALENERGY] ,uw);
    scalars[SizeOfHamiltonians+12]( TailProp[LOCALENERGY] ,uw);
//       int Rage(Reptile->Age);
//       int Bage=Rage;
//
// //       RealType tmpE=0.0;
//       RealType tmpV=0.0;
// //       RealType tmpF=0.0;
// //       RealType localE=0.0;
// //       KEconst = (*(Reptile->begin()))->Drift.size()  * 1.5 * OneOverTau;
//       for( MultiChain::iterator Bit = Reptile->begin();Bit != (Reptile->end());Bit++){
// //         localE += 0.5*( (*Bit)->deltaRSquared[0] + (*Bit)->deltaRSquared[1]);
// //         tmpF+= 0.5*(*Bit)->deltaRSquared[2];
// //         tmpE+=(*Bit)->getPropertyBase(i)[LOCALENERGY];
//         tmpV+=(*Bit)->getPropertyBase(i)[LOCALPOTENTIAL];
//         Bage= std::min((*Bit)->stepmade,Bage);
//       };
//
// //       tmpF-=0.25*Reptile->back()->deltaRSquared[2];
// //       tmpF-=0.25*Reptile->front()->deltaRSquared[2];
// //       tmpE-=0.5*Reptile->back()->getPropertyBase(i)[LOCALENERGY];
// //       tmpE-=0.5*Reptile->front()->getPropertyBase(i)[LOCALENERGY];
//       tmpV-=0.5*Reptile->back()->getPropertyBase(i)[LOCALPOTENTIAL];
//       tmpV-=0.5*Reptile->front()->getPropertyBase(i)[LOCALPOTENTIAL];
//       tmpV *= -pNorm*Tau;
//
// //       localE *= -0.5*OneOverTau*OneOverTau;
// //       localE += KEconst* (Reptile->Last) + tmpF ;
// //       localE += tmpE;
//
//       scalars[SizeOfHamiltonians+1](tmpV,uw);
//       scalars[SizeOfHamiltonians+2](eloc*tmpV,uw);
//
//       scalars[SizeOfHamiltonians+4](Rage-Bage,1.0);
    ///This is the center bead energy using PIMC stuff.
//       localE *= 1.0/(Reptile->Last);
//       scalars[SizeOfHamiltonians+5]( localE ,uw);
    int Rage(Reptile->Age);
    int Bage=Rage;
    RealType tmpV_head=0.0;
    RealType tmpV_tail=0.0;
    RealType tmpF=0.0;
    int maxtouch=0;
    MultiChain::iterator Bit = Reptile->begin();
    MultiChain::iterator Lit = Reptile->end();
    MultiChain::iterator Endit = Reptile->end();
    Endit--;
    //truncated version of estimator
    Lit--;
    for(int j=0 ; ( (j<truncLength[0]) && (Bit != Endit) ); Bit++,Lit--,j++)
    {
      tmpV_head+= (*Bit)->getPropertyBase(i)[LOCALPOTENTIAL];
      tmpV_tail+= (*Lit)->getPropertyBase(i)[LOCALPOTENTIAL];
      if ((*Bit)->timesTouched>maxtouch)
        maxtouch = (*Bit)->timesTouched;
      if ((*Bit)->stepmade < Bage)
        Bage = (*Bit)->stepmade;
    };
    RealType Ppref = (-pNorm*Tau);
    RealType tmpVsum = (energy_tail*tmpV_head + energy_head*tmpV_tail)*Ppref;
    scalars[SizeOfHamiltonians+3]((tmpV_head+tmpV_tail)*Ppref,uw);
    scalars[SizeOfHamiltonians+4](0.0,uw);
    scalars[SizeOfHamiltonians+5]( tmpVsum,uw);
    //continue sum for comparison to truncated version
    for( ; Bit != Endit; Bit++ )
    {
      tmpV_head+= (*Bit)->getPropertyBase(i)[LOCALPOTENTIAL];
//         tmpV_tail+= (*Lit)->getPropertyBase(i)[LOCALPOTENTIAL];
      if ((*Bit)->timesTouched>maxtouch)
        maxtouch = (*Bit)->timesTouched;
      if ((*Bit)->stepmade < Bage)
        Bage = (*Bit)->stepmade;
    };
    RealType tmpV = tmpV_head*Ppref;
    RealType tmpEVsum = 0.5*(energy_head  + energy_tail)*tmpV;
    scalars[SizeOfHamiltonians+1](tmpV,uw);
    scalars[SizeOfHamiltonians+2](tmpEVsum,uw);
    scalars[SizeOfHamiltonians+7](Rage-Bage,1.0);
    scalars[SizeOfHamiltonians+8](maxtouch,1.0);
    //calculate correlation between head energy and energy at some point in chain to truncate correctly
    Bit= Reptile->begin();
    ObsEvals+=1.0;
    RealType Opf=1.0/ObsEvals;
    for(int j=0; Bit != (Reptile->end()); Bit++,j++)
    {
      ObsCont[j]+=(*Bit)->getPropertyBase(i)[LOCALENERGY];
      ObsContAvg[j]=ObsCont[j]*Opf;
      ObsCont2[j]+=(*Bit)->getPropertyBase(i)[LOCALPOTENTIAL];
      ObsContAvg2[j]=ObsCont[j]*Opf;
    };
    int AC(-1);
    RealType tmpE=1.0;
    RealType Hac=HeadProp[LOCALENERGY]-ObsContAvg[0];
    Bit = Reptile->begin();
    Bit++;
    for(int j=1; (Bit!=(Reptile->end()))&&(tmpE>=0.0) ; Bit++,j++)
    {
      tmpE=Hac*( (*Bit)->getPropertyBase(i)[LOCALENERGY]-ObsContAvg[j] );
      AC+=1;
    };
    scalars[SizeOfHamiltonians+13](AC,uw);
    int AC2(-1);
    tmpE= -1.0;
    Bit = Reptile->begin();
    if (Hac*((*Bit)->getPropertyBase(i)[LOCALPOTENTIAL]-ObsContAvg2[0]) >0)
      tmpE=1.0;
    Bit++;
    for(int j=1; (Bit!=(Reptile->end()))&&(tmpE>=0.0) ; Bit++,j++)
    {
      tmpE=Hac*( (*Bit)->getPropertyBase(i)[LOCALPOTENTIAL]-ObsContAvg2[j] );
      AC2+=1;
    }
    scalars[SizeOfHamiltonians+14](AC2,uw);
//       scalars[SizeOfHamiltonians+15](CenProp[Findex+FirstHamiltonian],uw);
//       scalars[SizeOfHamiltonians+16](CenProp[Findex+1+FirstHamiltonian],uw);
  }
}

void MJPolymerEstimator::registerObservables(std::vector<observable_helper*>& h5dec, hid_t gid)
{
  //IMPLEMENT for hdf5
}

}

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
    
    
    
    
    
    
    
    
#include "Estimators/HFDHE2PolymerEstimator.h"
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
HFDHE2PolymerEstimator::HFDHE2PolymerEstimator(QMCHamiltonian& h, int hcopy, MultiChain* polymer)
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
  elocal_name.push_back("TruncSumPot");
  elocal_name.push_back("TruncElSumPot");
  elocal_name.push_back("Virial");
  elocal_name.push_back("MaxAge");
  elocal_name.push_back("MaxTouch");
  elocal_name.push_back("centerV");
  elocal_name.push_back("Edecorr");
  elocal_name.push_back("Vdecorr");
  scalars.resize(SizeOfHamiltonians+11);
  scalars_saved=scalars;
  pNorm=0.0;
  Pindex = h.getHamiltonian("HePress")->myIndex;
  HDFHE2index = h.getHamiltonian("HFDHE2")->myIndex;
  ObsEvals=0;
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
HFDHE2PolymerEstimator::add2Record(RecordNamedProperty<RealType>& record)
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

void HFDHE2PolymerEstimator::accumulate(const MCWalkerConfiguration& W
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
    scalars[SizeOfHamiltonians+5]( (2.0*eloc-2.0*CenProp[LOCALPOTENTIAL])*pNorm + CenProp[Pindex+1+FirstHamiltonian] ,uw);
    int Rage(Reptile->Age);
    int Bage=Rage;
    RealType tmpV=0.0;
    RealType tmpF=0.0;
//       if (truncLength<0){
//         for( MultiChain::iterator Bit = Reptile->begin();Bit != (Reptile->end());Bit++){
//           tmpF+= (*Bit)->getPropertyBase(i)[Pindex+2+FirstHamiltonian];
//           tmpV+=(*Bit)->getPropertyBase(i)[LOCALPOTENTIAL];
//           Bage= std::min((*Bit)->stepmade,Bage);
//         };
//
//       tmpF-=0.5*Reptile->back()->getPropertyBase(i)[Pindex+2+FirstHamiltonian];
//       tmpF-=0.5*Reptile->front()->getPropertyBase(i)[Pindex+2+FirstHamiltonian];
//       tmpV-=0.5*Reptile->back()->getPropertyBase(i)[LOCALPOTENTIAL];
//       tmpV-=0.5*Reptile->front()->getPropertyBase(i)[LOCALPOTENTIAL];
//       tmpV *= -2.0*pNorm*Tau;
//       tmpF *= Tau;
//
//       scalars[SizeOfHamiltonians+1](tmpV+tmpF,uw);
//       scalars[SizeOfHamiltonians+2](eloc*(tmpV+tmpF),uw);
//       scalars[SizeOfHamiltonians+3](tmpV+tmpF,uw);
//       scalars[SizeOfHamiltonians+4](eloc*(tmpV+tmpF),uw);
//       } else {
    int maxtouch=0;
    MultiChain::iterator Bit = Reptile->begin();
    MultiChain::iterator Lit = Reptile->end();
    MultiChain::iterator Endit = Reptile->end();
    Endit--;
    //truncated version of estimator
    Lit--;
    for(int j=0 ; ( (j<truncLength) && (Bit != Endit) ); Bit++,Lit--,j++)
    {
      tmpF+= 0.5*(*Bit)->getPropertyBase(i)[Pindex+2+FirstHamiltonian];
      tmpV+= 0.5*(*Bit)->getPropertyBase(i)[LOCALPOTENTIAL];
      tmpF+= 0.5*(*Lit)->getPropertyBase(i)[Pindex+2+FirstHamiltonian];
      tmpV+= 0.5*(*Lit)->getPropertyBase(i)[LOCALPOTENTIAL];
      if ((*Bit)->timesTouched>maxtouch)
        maxtouch = (*Bit)->timesTouched;
      if ((*Bit)->stepmade < Bage)
        Bage = (*Bit)->stepmade;
    };
    RealType tmpval = (tmpV*(-2.0*pNorm) + tmpF)*Tau;
    scalars[SizeOfHamiltonians+3](tmpval,uw);
    scalars[SizeOfHamiltonians+4](eloc*tmpval,uw);
    //continue sum for comparison to truncated version
    for( ; Bit != Endit; Bit++,Lit--)
    {
      tmpF+= 0.5*(*Bit)->getPropertyBase(i)[Pindex+2+FirstHamiltonian];
      tmpV+= 0.5*(*Bit)->getPropertyBase(i)[LOCALPOTENTIAL];
      tmpF+= 0.5*(*Lit)->getPropertyBase(i)[Pindex+2+FirstHamiltonian];
      tmpV+= 0.5*(*Lit)->getPropertyBase(i)[LOCALPOTENTIAL];
      if ((*Bit)->timesTouched>maxtouch)
        maxtouch = (*Bit)->timesTouched;
      if ((*Bit)->stepmade < Bage)
        Bage = (*Bit)->stepmade;
    };
    tmpV *= -2.0*pNorm*Tau;
    tmpF *= Tau;
    scalars[SizeOfHamiltonians+1](tmpV+tmpF,uw);
    scalars[SizeOfHamiltonians+2](eloc*(tmpV+tmpF),uw);
    scalars[SizeOfHamiltonians+6](Rage-Bage,1.0);
    scalars[SizeOfHamiltonians+7](maxtouch,1.0);
    scalars[SizeOfHamiltonians+8](CenProp[LOCALPOTENTIAL] ,uw);
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
    scalars[SizeOfHamiltonians+9](AC,uw);
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
    };
    scalars[SizeOfHamiltonians+10](AC2,uw);
  }
}

void HFDHE2PolymerEstimator::registerObservables(std::vector<observable_helper*>& h5dec, hid_t gid)
{
  //IMPLEMENT for hdf5
}

}

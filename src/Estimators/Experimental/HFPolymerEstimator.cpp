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
    
    
    
    
    
    
    
    



#include <sstream>

#include "Estimators/HFPolymerEstimator.h"
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
HFPolymerEstimator::HFPolymerEstimator(QMCHamiltonian& h, int hcopy, MultiChain* polymer): H(h)
{
  Reptile=polymer;
  NumCopies=hcopy;
//     Hpointer = &h;
  SizeOfHamiltonians = h.sizeOfObservables();
  FirstHamiltonian = h.startIndex();
//     std::cout <<"size of Hamiltonian "<<SizeOfHamiltonians<<" First one "<<FirstHamiltonian<< std::endl;
  elocal_name.push_back("LocalEnergy");
  elocal_name.push_back("HeadEnergy");
  elocal_name.push_back("TailEnergy");
  elocal_name.push_back(h.getObservableName(0));
//     for(int i=0; i < SizeOfHamiltonians; i++)
//     {
//       elocal_name.push_back(h.getObservableName(i));
//     };
};


void HFPolymerEstimator::add_HF_Observables(int nchains, int strd)
{
  stride=strd;
//     never FW on kinetic energy
  resize_HF(nchains,SizeOfHamiltonians-1);
  nobs=0;
  if(nchains%stride>0)
    nobs=1+(nchains-1)/stride;
  else
    nobs=nchains/stride;
  app_log()<<"Writing observables each "<<stride<<" steps"<< std::endl;
  app_log()<<"Writing "<<nobs<<" observables"<< std::endl;
  int nadded(elocal_name.size());
  for(int j=0; j < nobs; j++)
    for(int i=1; i < SizeOfHamiltonians; i++)
    {
      std::ostringstream ss;
      ss << "BR_"<<H.getObservableName(i)<<"_"<<stride*j;
      elocal_name.push_back(ss.str());
    }
  for(int j=0; j < nobs; j++)
    for(int i=1; i < SizeOfHamiltonians; i++)
    {
      std::ostringstream ss;
      ss << "ED_"<<H.getObservableName(i)<<"_"<<stride*j;
      elocal_name.push_back(ss.str());
      ss.str("");
      ss << "DH_"<<H.getObservableName(i)<<"_"<<stride*j;
      elocal_name.push_back(ss.str());
      ss.str("");
      ss << "DT_"<<H.getObservableName(i)<<"_"<<stride*j;
      elocal_name.push_back(ss.str());
    }
  scalars.resize(elocal_name.size()); //SizeOfHamiltonians+elocal_name.size()-nadded);
  scalars_saved=scalars;
};


HFPolymerEstimator::HFPolymerEstimator(const HFPolymerEstimator& mest):
  PolymerEstimator(mest), H(mest.H)
//       Reptile(mest.Reptile)
{
  Reptile=mest.Reptile;
  NumCopies=mest.NumCopies;
  //d_data.resize(mest.d_data.size());
}

ScalarEstimatorBase* HFPolymerEstimator::clone()
{
  return new HFPolymerEstimator(*this);
}

/**  add the local energy, variance and all the Hamiltonian components to the scalar record container
 *@param record storage of scalar records (name,value)
 */
void
HFPolymerEstimator::add2Record(RecordNamedProperty<RealType>& record)
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

void HFPolymerEstimator::accumulate(const MCWalkerConfiguration& W
                                    , WalkerIterator first, WalkerIterator last, RealType wgt)
{
  for(int i=0; i<1; i++)
  {
//       RealType uw(Reptile->UmbrellaWeight[i]);
    RealType uw(1);
    //Adding all parts from the Hamiltonian
    RealType* restrict HeadProp(Reptile->front()->getPropertyBase(i));
    RealType* restrict TailProp(Reptile->back()->getPropertyBase(i));
    RealType eloc = 0.5*( HeadProp[LOCALENERGY] + TailProp[LOCALENERGY]);
    int oindx=0;
    scalars[oindx++](eloc,uw);
    scalars[oindx++](HeadProp[LOCALENERGY],uw);
    scalars[oindx++](TailProp[LOCALENERGY],uw);
//       for(int obsi=0;obsi<SizeOfHamiltonians ;obsi++){
    scalars[oindx++]( 0.5*( HeadProp[0+FirstHamiltonian] + TailProp[0+FirstHamiltonian]) , uw);
//       };
//    observables in reptile.
    ObsSumR=0;
    ObsSumL=0;
    int xi(0);
    MultiChain::iterator Bit_first = Reptile->begin();
    MultiChain::iterator Bit_last = Reptile->end();
    for( MultiChain::iterator Bit = Bit_first; Bit != Bit_last; Bit++,xi++)
    {
      for(int obsi=1; obsi<SizeOfHamiltonians ; obsi++)
      {
        if(xi%stride==0)
          scalars[oindx++]((*Bit)->getPropertyBase(i)[obsi+FirstHamiltonian], uw);
        ObsSumL(xi+1,obsi-1)=(*Bit)->getPropertyBase(i)[obsi+FirstHamiltonian] + ObsSumL(xi,obsi-1);
      }
    }
    Bit_last--;
    Bit_first--;
    xi--;
    for( MultiChain::iterator Bit = Bit_last; Bit != Bit_first; Bit--,xi--)
      for(int obsi=1; obsi<SizeOfHamiltonians ; obsi++)
        ObsSumR(xi-1,obsi-1)=ObsSumR(xi,obsi-1)+(*Bit)->getPropertyBase(i)[obsi+FirstHamiltonian];
    xi=0;
    Bit_last++;
    Bit_first++;
    for( MultiChain::iterator Bit = Bit_first; Bit != Bit_last; Bit++,xi++)
    {
      if(xi%stride==0)
      {
        for(int obsi=1; obsi<SizeOfHamiltonians ; obsi++)
        {
          //           E*D
          scalars[oindx++](Tau*(TailProp[LOCALENERGY]*ObsSumL(xi,obsi-1)+HeadProp[LOCALENERGY]*ObsSumR(xi,obsi-1)), uw);
          //           L and R
          scalars[oindx++](Tau*ObsSumL(xi,obsi-1), uw);
          scalars[oindx++](Tau*ObsSumR(xi,obsi-1), uw);
        }
      }
    }
//       //Center Pressure
//       RealType* restrict CenProp(Reptile->center()->getPropertyBase(i));
//       scalars[SizeOfHamiltonians+6]( (2.0*eloc-CenProp[LOCALPOTENTIAL])*pNorm ,uw);
//       scalars[SizeOfHamiltonians+9]( CenProp[LOCALPOTENTIAL] ,uw);
//       scalars[SizeOfHamiltonians+10]( HeadProp[LOCALENERGY] * TailProp[LOCALENERGY] ,uw);
//       RealType energy_head = HeadProp[LOCALENERGY];
//       RealType energy_tail = TailProp[LOCALENERGY];
//       scalars[SizeOfHamiltonians+11]( HeadProp[LOCALENERGY] ,uw);
//       scalars[SizeOfHamiltonians+12]( TailProp[LOCALENERGY] ,uw);
//
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
//       int Rage(Reptile->Age);
//       int Bage=Rage;
//
//       RealType tmpV_head=0.0;
//       RealType tmpV_tail=0.0;
//       RealType tmpF=0.0;
//       int maxtouch=0;
//       MultiChain::iterator Bit = Reptile->begin();
//       MultiChain::iterator Lit = Reptile->end();
//       MultiChain::iterator Endit = Reptile->end();
//       Endit--;
//
//         //truncated version of estimator
//       Lit--;
//       for(int j=0 ;( (j<truncLength[0]) && (Bit != Endit) );Bit++,Lit--,j++){
//         tmpV_head+= (*Bit)->getPropertyBase(i)[LOCALPOTENTIAL];
//         tmpV_tail+= (*Lit)->getPropertyBase(i)[LOCALPOTENTIAL];
//         if ((*Bit)->timesTouched>maxtouch) maxtouch = (*Bit)->timesTouched;
//         if ((*Bit)->stepmade < Bage) Bage = (*Bit)->stepmade;
//       };
//       RealType Ppref = (-pNorm*Tau);
//
//
//
//
//
//       RealType tmpVsum = (energy_tail*tmpV_head + energy_head*tmpV_tail)*Ppref;
//       scalars[SizeOfHamiltonians+3]((tmpV_head+tmpV_tail)*Ppref,uw);
//       scalars[SizeOfHamiltonians+4](0.0,uw);
//       scalars[SizeOfHamiltonians+5]( tmpVsum,uw);
//
//         //continue sum for comparison to truncated version
//       for( ;Bit != Endit;Bit++ ){
//         tmpV_head+= (*Bit)->getPropertyBase(i)[LOCALPOTENTIAL];
// //         tmpV_tail+= (*Lit)->getPropertyBase(i)[LOCALPOTENTIAL];
//         if ((*Bit)->timesTouched>maxtouch) maxtouch = (*Bit)->timesTouched;
//         if ((*Bit)->stepmade < Bage) Bage = (*Bit)->stepmade;
//       };
//
//       RealType tmpV = tmpV_head*Ppref;
//       RealType tmpEVsum = 0.5*(energy_head  + energy_tail)*tmpV;
//       scalars[SizeOfHamiltonians+1](tmpV,uw);
//       scalars[SizeOfHamiltonians+2](tmpEVsum,uw);
//
//
//       scalars[SizeOfHamiltonians+7](Rage-Bage,1.0);
//       scalars[SizeOfHamiltonians+8](maxtouch,1.0);
//
//       //calculate correlation between head energy and energy at some point in chain to truncate correctly
//       Bit= Reptile->begin();
//
//       ObsEvals+=1.0;
//       RealType Opf=1.0/ObsEvals;
//       for(int j=0;Bit != (Reptile->end());Bit++,j++){
//         ObsCont[j]+=(*Bit)->getPropertyBase(i)[LOCALENERGY];
//         ObsContAvg[j]=ObsCont[j]*Opf;
//         ObsCont2[j]+=(*Bit)->getPropertyBase(i)[LOCALPOTENTIAL];
//         ObsContAvg2[j]=ObsCont[j]*Opf;
//       };
//       int AC(-1);
//       RealType tmpE=1.0;
//       RealType Hac=HeadProp[LOCALENERGY]-ObsContAvg[0];
//
//       Bit = Reptile->begin();
//       Bit++;
//       for(int j=1; (Bit!=(Reptile->end()))&&(tmpE>=0.0) ;Bit++,j++){
//         tmpE=Hac*( (*Bit)->getPropertyBase(i)[LOCALENERGY]-ObsContAvg[j] );
//         AC+=1;
//       };
//       scalars[SizeOfHamiltonians+13](AC,uw);
//
//       int AC2(-1);
//       tmpE= -1.0;
//       Bit = Reptile->begin();
//       if (Hac*((*Bit)->getPropertyBase(i)[LOCALPOTENTIAL]-ObsContAvg2[0]) >0) tmpE=1.0;
//       Bit++;
//       for(int j=1; (Bit!=(Reptile->end()))&&(tmpE>=0.0) ;Bit++,j++)
//       {
//         tmpE=Hac*( (*Bit)->getPropertyBase(i)[LOCALPOTENTIAL]-ObsContAvg2[j] );
//         AC2+=1;
//       }
//
//       scalars[SizeOfHamiltonians+14](AC2,uw);
//       scalars[SizeOfHamiltonians+15](CenProp[Findex+FirstHamiltonian],uw);
//       scalars[SizeOfHamiltonians+16](CenProp[Findex+1+FirstHamiltonian],uw);
  }
}

void HFPolymerEstimator::registerObservables(std::vector<observable_helper*>& h5dec, hid_t gid)
{
  //IMPLEMENT for hdf5
}

}

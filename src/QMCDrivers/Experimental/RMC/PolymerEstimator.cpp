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
#include "QMCDrivers/PolymerEstimator.h"

namespace qmcplusplus
{

PolymerEstimator::PolymerEstimator(PolymerChain& reptile, int npsi):
  Reptile(reptile),
  Middle(0), Counter(0), nPsi(npsi),
  fout(0), OutLocEn(0),OutPotEn(0)
{
  Middle=Reptile.Middle;
  PEavg.resize(Middle+1,0.0);
  PE2.resize(Middle+1,0.0);
}

void PolymerEstimator::clean()
{
  if(fout)
  {
    delete fout;
    delete OutLocEn;
    delete OutPotEn;
    fout=0;
    OutLocEn=0;
    OutPotEn=0;
  }
}

void PolymerEstimator::resetReportSettings(const string& aname)
{
  clean();
  std::string filename(aname);
  filename.append(".pe.dat");
  fout = new ofstream(filename.c_str());
  //filename=aname+".Eloc.bin";
  filename="Eloc.bin";
  OutLocEn=new ofstream(filename.c_str(),ios::binary);
  //filename=aname+".Epot.bin";
  filename="Epot.bin";
  OutPotEn=new ofstream(filename.c_str(),ios::binary);
  int uno=1;
  int isize=sizeof(int);
  int dsize=sizeof(double);
  ReptileLength=Reptile.size();
  EpotLength=ReptileLength*nPsi;
  EpotSize=EpotLength*dsize;
  ElocSize=nPsi*dsize;
  OutPotEn->write(reinterpret_cast<char*>(&nPsi),isize);
  OutPotEn->write(reinterpret_cast<char*>(&ReptileLength),isize);
  OutLocEn->write(reinterpret_cast<char*>(&nPsi),isize);
  OutLocEn->write(reinterpret_cast<char*>(&uno),isize);
  AvgLocalEnergy.resize(nPsi);
  AvgPotentialEnergy.resize(EpotLength);
  TotalWeight.resize(nPsi);
  reset();
}

void PolymerEstimator::report(int iter)
{
  (*fout) << iter;
  double wgtinv = 1.0/static_cast<double>(Counter);
  for(int i=0; i<PEavg.size(); i++)
  {
    double eavg = PEavg[i]*wgtinv;
    (*fout)  <<  setw(15) << eavg;
  }
  (*fout) << endl;
  //double oneoversteps=1.0/static_cast<double>(Counter);
  for(int i=0; i<nPsi; i++)
  {
    double wgtinv = 1.0/TotalWeight[i];
    AvgLocalEnergy[i] *= wgtinv;
    AvgPotentialEnergy[i] *= wgtinv;
  }
  OutLocEn->write(reinterpret_cast<char*>(AvgLocalEnergy.data()),ElocSize);
  OutLocEn->flush();
  OutPotEn->write(reinterpret_cast<char*>(AvgPotentialEnergy.data()),EpotSize);
  OutPotEn->flush();
  reset();
}

void PolymerEstimator::accumulate()
{
  Counter++;
  register double e;
  for(int k=0,j=Reptile.size()-1; k<Middle; k++,j--)
  {
    e =0.5*(Reptile[k]->Properties(LOCALPOTENTIAL)+Reptile[j]->Properties(LOCALPOTENTIAL));
    PEavg[k] += e;
    PE2[k] += e*e;
  }
  e= Reptile[Middle]->Properties(LOCALPOTENTIAL);
  PEavg[Middle] += e;
  PE2[Middle] += e*e;
  // Accumulate Block Averages
  for(int ipsi=0; ipsi<nPsi; ipsi++)
  {
    TotalWeight[ipsi]+=Reptile.UmbrellaWeight[ipsi];
    double w=0.5*Reptile.UmbrellaWeight[ipsi];
    AvgLocalEnergy[ipsi]
    += (w*(Reptile.front()->Properties(ipsi,LOCALENERGY)+Reptile.back()->Properties(ipsi,LOCALENERGY)));
  }
  //Loop over beads. Unnecessary, i.e. sum every 5/10 beads.
  PolymerChain::iterator beadIter(Reptile.begin());
  int ibead=0;
  while(ibead<ReptileLength)
  {
    for(int ipsi=0; ipsi<nPsi; ipsi++)
    {
      AvgPotentialEnergy[ibead]
      +=Reptile.UmbrellaWeight[ipsi]*(*beadIter)->Properties(ipsi,LOCALPOTENTIAL);
    }
    ++ibead;
    ++beadIter;
  }
}
}
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 759 $   $Date: 2005-10-31 10:10:28 -0600 (Mon, 31 Oct 2005) $
 * $Id: PolymerEstimator.cpp 759 2005-10-31 16:10:28Z jnkim $
 ***************************************************************************/

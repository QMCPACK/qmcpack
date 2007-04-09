//////////////////////////////////////////////////////////////////
// (c) Copyright 2004- by Jeongnim Kim 
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#include "Estimators/PairCorrEstimator.h"
#include "Particle/DistanceTable.h"
#include "Particle/DistanceTableData.h"

namespace qmcplusplus {
  PairCorrEstimator::PairCorrEstimator(ParticleSet& source): 
    Symmetric(true),sourcePtcl(source)
  {
    myTable = DistanceTable::add(source);
    int ns=sourcePtcl.groups();

    vector<int> mask(ns*ns,-1);
    int ij=0;
    for(int i=0; i<ns; i++)
      for(int j=i; j<ns; j++,ij++) 
      {
        mask[j+i*ns]=ij;
        char fname[32];
        sprintf(fname,"gofr.%s_%d_%d.dat",myTable->Name.c_str(),i,j);
        fout.push_back(new ofstream(fname));
        fout[ij]->setf(ios::scientific, ios::floatfield);
        fout[ij]->precision(5);
      }

    NumPairTypes=ij;
    Centers=sourcePtcl.getTotalNum();

    PairID.resize(myTable->getTotNadj());
    for(int iat=0; iat<Centers; iat++) {
      for(int nn=myTable->M[iat]; nn<myTable->M[iat+1]; nn++)
      {
        PairID[nn]=mask[myTable->PairID[nn]];
      }
    }
    setBound(10,0.1);
  }

  PairCorrEstimator::PairCorrEstimator(const ParticleSet& source, ParticleSet& target):
    Symmetric(false),sourcePtcl(source)
  {
    myTable = DistanceTable::add(source,target);
    NumPairTypes=sourcePtcl.getSpeciesSet().getTotalNum(); 
    for(int i=0; i<NumPairTypes; i++) 
    {
      char fname[32];
      sprintf(fname,"gofr.%s_%s.dat",myTable->Name.c_str(),
          sourcePtcl.getSpeciesSet().speciesName[i].c_str());
      fout.push_back(new ofstream(fname));
      fout[i]->setf(ios::scientific, ios::floatfield);
      fout[i]->precision(5);
    }

    Centers=sourcePtcl.getTotalNum();
    PairID.resize(myTable->getTotNadj());
    for(int iat=0; iat<Centers; iat++) {
      for(int nn=myTable->M[iat]; nn<myTable->M[iat+1]; nn++)
      {
        PairID[nn]=sourcePtcl.GroupID[iat];
      }
    }
    setBound(10,0.1);
  }

  PairCorrEstimator::~PairCorrEstimator()
  {
  }

  void PairCorrEstimator::resetTargetParticleSet(ParticleSet& p)
  {
    if(Symmetric)
      myTable=DistanceTable::add(p);
    else
      myTable=DistanceTable::add(sourcePtcl,p);
  }

  void PairCorrEstimator::startAccumulate()
  {
    dCInst=0.0;
  }

  /** accumulate the observables */
  void PairCorrEstimator::accumulate(ParticleSet& p, RealType wgt)
  {
    for(int iat=0; iat<Centers; iat++) {
      for(int nn=myTable->M[iat]; nn<myTable->M[iat+1]; nn++)
      {
        if(myTable->r(nn)>=Dmax) continue;
        dCInst(PairID[nn],DeltaInv*myTable->r(nn))+=wgt;
      }
    }
  }

  /** reweight of the current cummulative  values */
  void PairCorrEstimator::stopAccumulate(RealType wgtinv)
  {
    //RealType norm=wgtinv/static_cast<RealType>(myTable->getTotNadj());
    //dCBlock += dCInst*norm;
    dCBlock += dCInst*wgtinv;
  }

  void PairCorrEstimator::stopBlock(RealType wgtinv)
  {
    for(int i=0; i<dCBlock.size1(); ++i)
    {
      RealType r=0.0;
      for(int j=0; j<dCBlock.size2(); ++j, r+=Delta)
        *fout[i] << setw(15) << r << setw(15) << wgtinv*dCBlock(i,j) << endl;
      *fout[i] << endl;
    }
  }

  void PairCorrEstimator::startBlock(int steps)
  {
    dCInst=0.0;
    dCBlock=0.0;
  }

  void PairCorrEstimator::setBound(RealType rmax, RealType dr)
  {
    Dmax=rmax;
    Delta=dr;
    DeltaInv=1.0/dr;
    int n=(Dmax)/Delta+1;
    dCInst.resize(NumPairTypes,n);
    dCBlock.resize(NumPairTypes,n);

  }

//  /////////////////////////////////////////////////////
//  PairCorrEstimator::PairCorrEstimator(ParticleSet& source): 
//    Symmetric(true),sourcePtcl(source), fout(0)
//  {
//    myTable = DistanceTable::add(source);
//    setBound(10,0.1);
//  }
//
//  PairCorrEstimator::PairCorrEstimator(const ParticleSet& source, ParticleSet& target):
//    Symmetric(false),sourcePtcl(source), fout(0)
//  {
//    myTable = DistanceTable::add(source,target);
//    setBound(10,0.1);
//  }
//
//  PairCorrEstimator::~PairCorrEstimator()
//  {
//  }
//
//  void PairCorrEstimator::resetTargetParticleSet(ParticleSet& p)
//  {
//    if(Symmetric)
//      myTable=DistanceTable::add(p);
//    else
//      myTable=DistanceTable::add(sourcePtcl,p);
//  }
//
//  void PairCorrEstimator::startAccumulate()
//  {
//    dCInst=0.0;
//  }
//
//  /** accumulate the observables */
//  void PairCorrEstimator::accumulate(ParticleSet& p, RealType wgt)
//  {
//    for(int i=0; i<myTable->getTotNadj(); i++)
//    {
//      if(myTable->r(i)<Dmax) dCInst[DeltaInv*myTable->r(i)]+=wgt;
//    }
//  }
//
//  /** reweight of the current cummulative  values */
//  void PairCorrEstimator::stopAccumulate(RealType wgtinv)
//  {
//    //RealType norm=wgtinv/static_cast<RealType>(myTable->getTotNadj());
//    //dCBlock += dCInst*norm;
//    dCBlock += dCInst*wgtinv;
//  }
//
//  void PairCorrEstimator::stopBlock(RealType wgtinv)
//  {
//    if(!fout)
//    {
//      char fname[32];
//      sprintf(fname,"gofr.%s.dat",myTable->Name.c_str());
//      fout = new ofstream(fname);
//      fout->setf(ios::scientific, ios::floatfield);
//      fout->precision(5);
//    }
//    RealType r=0.0;
//    for(int i=0; i<dCBlock.size(); i++, r+=Delta) 
//      *fout << setw(15) << r << setw(15) << wgtinv*dCBlock[i] << endl;
//    *fout << endl;
//  }
//
//  void PairCorrEstimator::startBlock(int steps)
//  {
//    dCInst=0.0;
//    dCBlock=0.0;
//  }
//
//  void PairCorrEstimator::setBound(RealType rmax, RealType dr)
//  {
//    Dmax=rmax;
//    Delta=dr;
//    DeltaInv=1.0/dr;
//    int n=(Dmax)/Delta+1;
//    dCInst.resize(n);
//    dCBlock.resize(n);
//  }
}

/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1415 $   $Date: 2006-10-23 11:51:53 -0500 (Mon, 23 Oct 2006) $
 * $Id: CompositeEstimatorBase.h 1415 2006-10-23 16:51:53Z jnkim $ 
 ***************************************************************************/

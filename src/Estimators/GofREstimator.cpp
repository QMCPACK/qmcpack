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
#include "Estimators/GofREstimator.h"
#include "Particle/DistanceTable.h"
#include "Particle/DistanceTableData.h"
#include "Utilities/IteratorUtility.h"
#define PRINT_DEBUG

namespace qmcplusplus {


  GofREstimator::GofREstimator(ParticleSet& source): 
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
        sprintf(fname,"%s_%d_%d",myTable->Name.c_str(),i,j);
        PairName.push_back(fname);
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

    setBound(0.1);
  }

  GofREstimator::GofREstimator(const ParticleSet& source, ParticleSet& target):
    Symmetric(false),sourcePtcl(source)
  {
    myTable = DistanceTable::add(source,target);
    NumPairTypes=sourcePtcl.getSpeciesSet().getTotalNum(); 
    for(int i=0; i<NumPairTypes; i++) 
    {
      char fname[32];
      sprintf(fname,"%s_%s",myTable->Name.c_str(),
          sourcePtcl.getSpeciesSet().speciesName[i].c_str());
      PairName.push_back(fname);
    }
    Centers=sourcePtcl.getTotalNum();
    PairID.resize(myTable->getTotNadj());
    for(int iat=0; iat<Centers; iat++) {
      for(int nn=myTable->M[iat]; nn<myTable->M[iat+1]; nn++)
      {
        PairID[nn]=sourcePtcl.GroupID[iat];
      }
    }

    setBound(0.1);
  }

  GofREstimator::~GofREstimator()
  {
    delete_iter(gofr.begin(),gofr.end());
  }

  void GofREstimator::resetTargetParticleSet(ParticleSet& p)
  {
    if(Symmetric)
      myTable=DistanceTable::add(p);
    else
      myTable=DistanceTable::add(sourcePtcl,p);
  }

  void GofREstimator::open(hid_t hroot)
  {
  //  if(GroupID<0)
  //  {
  //    Title="pc_"+myTable->Name;
  //    GroupID = H5Gcreate(hroot,Title.c_str(),0);
  //    v_h = new AppendData(gofr);
  //    v2_h = new AppendData(gofr2);
  //  }
  }

  void GofREstimator::close()
  {
//    if(GroupID>-1)
//    {
//      delete v_h;
//      delete v2_h;
//      H5Gclose(GroupID);
//      GroupID=-1;
//    }
  }

  /** ready to accumulate the measurements over the walkers
   */
  void GofREstimator::startAccumulate()
  {
    gofrInst=0.0;
  }

  /** accumulate the observables for a walker image*/
  void GofREstimator::accumulate(ParticleSet& p)
  {
    for(int iat=0; iat<Centers; iat++) {
      for(int nn=myTable->M[iat]; nn<myTable->M[iat+1]; nn++)
      {
        if(myTable->r(nn)>=Dmax) continue;
        //need a better average
        int ig=static_cast<int>(DeltaInv*myTable->r(nn));
        gofrInst(PairID[nn],ig)+=1;
      }
    }
  }

  /** add gofrInst which contains sum over walkers */
  void GofREstimator::stopAccumulate(RealType wgtinv)
  {
    for(int p=0; p<NumPairTypes; p++) 
    {
      //gofr[p]->accumulate(gofrInst[p],wgtinv);
      gofr[p]->accumulate(gofrInst[p],normFactor.begin());
    }
  }

  void GofREstimator::startBlock(int steps)
  {
    for(int p=0; p<NumPairTypes; p++) gofr[p]->reset();
  }

  /** save the block average */
  void GofREstimator::stopBlock(RealType wgtnorm, RealType errnorm)
  {
#if defined(PRINT_DEBUG)
    for(int p=0; p<NumPairTypes; p++) 
    {
      gofr[p]->takeBlockAverage(wgtnorm);
      ofstream fout(PairName[p].c_str());
      RealType r=0;
      for(int i=0; i< NumBins; i++, r+=Delta)
      {
        RealType a=gofr[p]->d_sum[i]*wgtnorm;
        fout << r << " " 
          << gofr[p]->d_data[2*i] << " " 
          << a << " " 
          << (gofr[p]->d_sum2[i]*errnorm-a*a) << endl;;
      }
    }
#endif
//    v_h->write(GroupID,"v");
//    v2_h->write(GroupID,"v2");
  }

  void GofREstimator::setBound(RealType dr)
  {
    RealType vnorm=1.0;
    if(sourcePtcl.Lattice.SuperCellEnum) 
    {
      Dmax=sourcePtcl.Lattice.LR_rc;
      /** normalizaton factor */
      vnorm=4.0*M_PI*myTable->size(DistanceTableData::SourceIndex)*myTable->size(DistanceTableData::VisitorIndex);
      vnorm=sourcePtcl.Lattice.Volume/vnorm;
    }
    else
    {
      Dmax=10.0; //choose a sensible number
    }

    //Dmax=rmax;
    Delta=dr;
    DeltaInv=1.0/dr;
    NumBins=static_cast<int>((Dmax)*DeltaInv+1);

    normFactor.resize(NumBins,0.0);
    RealType r=Delta;
    for(int i=1; i<NumBins; i++, r+=Delta) normFactor[i]=vnorm/(r*r); 

    gofrInst.resize(NumPairTypes,NumBins);
    delete_iter(gofr.begin(),gofr.end());
    gofr.resize(NumPairTypes,0);
    for(int i=0; i<NumPairTypes; i++)
    {
      gofr[i]=new VectorEstimatorType(NumBins);
      gofr[i]->Name=PairName[i];
    }
  }
}

/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1415 $   $Date: 2006-10-23 11:51:53 -0500 (Mon, 23 Oct 2006) $
 * $Id: CompositeEstimatorBase.h 1415 2006-10-23 16:51:53Z jnkim $ 
 ***************************************************************************/

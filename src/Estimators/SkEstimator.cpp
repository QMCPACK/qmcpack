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
#include "Estimators/SkEstimator.h"
#include "Particle/DistanceTable.h"
#include "Particle/DistanceTableData.h"
#include "Utilities/IteratorUtility.h"
#include "LongRange/StructFact.h"
#define PRINT_DEBUG

namespace qmcplusplus {


  SkEstimator::SkEstimator(ParticleSet& source) 
  {
    NumSpecies=source.getSpeciesSet().getTotalNum();
    NumK=source.SK->KLists.numk;
    OneOverN=1.0/static_cast<RealType>(source.getTotalNum());
    Kshell=source.SK->KLists.kshell;
    MaxKshell=Kshell.size()-1;
    SkInst.resize(NumK);
    RhokTot.resize(NumK);
    Sk.resize(NumK);

    Kmag.resize(MaxKshell);
    OneOverDnk.resize(MaxKshell);
    for(int ks=0, k=0; ks<MaxKshell; ks++)
    {
      Kmag[ks]=std::sqrt(source.SK->KLists.ksq[Kshell[ks]]);
      OneOverDnk[ks]=1.0/static_cast<RealType>(Kshell[ks+1]-Kshell[ks]);
    }

#if defined(PRINT_DEBUG)
    ofstream fout("Sk.dat");
#endif
  }


  SkEstimator::~SkEstimator()
  {
  }

  void SkEstimator::resetTargetParticleSet(ParticleSet& p)
  { }

  void SkEstimator::open(hid_t hroot)
  {
  //  if(GroupID<0)
  //  {
  //    Title="pc_"+myTable->Name;
  //    GroupID = H5Gcreate(hroot,Title.c_str(),0);
  //    v_h = new AppendData(gofr);
  //    v2_h = new AppendData(gofr2);
  //  }
  }

  void SkEstimator::close()
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
  void SkEstimator::startAccumulate()
  {
    SkInst=0.0;
  }

  /** accumulate the observables for a walker image*/
  void SkEstimator::accumulate(ParticleSet& p)
  {
    //sum over species
    std::copy(p.SK->rhok[0],p.SK->rhok[0]+NumK,RhokTot.begin());
    for(int i=1; i<NumSpecies; i++)
      accumulate_elements(p.SK->rhok[i],RhokTot.begin(),NumK);

    Vector<ComplexType>::const_iterator iit(RhokTot.begin());
    Vector<RealType>::iterator oit(SkInst.begin());
    for(int k=0; k<NumK; k++,iit++,oit++)
    {
      (*oit)+=(*iit).real()*(*iit).real()+(*iit).imag()*(*iit).imag();
    }
  }

  /** add gofrInst which contains sum over walkers */
  void SkEstimator::stopAccumulate(RealType wgtinv)
  {
    Sk.accumulate(SkInst.begin(),wgtinv*OneOverN);
  }

  void SkEstimator::startBlock(int steps)
  {
    Sk.reset();
  }

  /** save the block average */
  void SkEstimator::stopBlock(RealType wgtnorm, RealType errnorm)
  {
#if defined(PRINT_DEBUG)
    Sk.takeBlockAverage(wgtnorm);
    ofstream fout("Sk.dat",ios::app);
    for(int ks=0, k=0; ks<MaxKshell; ks++)
    {
      RealType s=0.0;
      for(;k<Kshell[ks+1]; k++) s += Sk.d_sum[k]; 
      fout << Kmag[ks] << " " << s*OneOverDnk[ks] << endl;
    }
    fout << endl;
#endif
//    v_h->write(GroupID,"v");
//    v2_h->write(GroupID,"v2");
  }
}

/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1415 $   $Date: 2006-10-23 11:51:53 -0500 (Mon, 23 Oct 2006) $
 * $Id: CompositeEstimatorBase.h 1415 2006-10-23 16:51:53Z jnkim $ 
 ***************************************************************************/

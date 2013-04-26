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
//#define PRINT_DEBUG

namespace qmcplusplus
{


SkEstimator::SkEstimator(ParticleSet& source)
{
  Title="sk_"+source.getName();
  NumSpecies=source.getSpeciesSet().getTotalNum();
  NumK=source.SK->KLists.numk;
  OneOverN=1.0/static_cast<RealType>(source.getTotalNum());
  Kshell=source.SK->KLists.kshell;
  MaxKshell=Kshell.size()-1;
  SkInst.resize(NumK);
  RhokTot.resize(NumK);
  Kmag.resize(MaxKshell);
  OneOverDnk.resize(MaxKshell);
  for(int ks=0, k=0; ks<MaxKshell; ks++)
  {
    Kmag[ks]=std::sqrt(source.SK->KLists.ksq[Kshell[ks]]);
    OneOverDnk[ks]=1.0/static_cast<RealType>(Kshell[ks+1]-Kshell[ks]);
  }
  nList.push_back(Title);
  dList.push_back(new VectorEstimatorType(Title,NumK));
#if defined(PRINT_DEBUG)
  ofstream fout("Sk.dat");
#endif
}

SkEstimator::SkEstimator(const SkEstimator& a):
  NumSpecies(a.NumSpecies),NumK(a.NumK),MaxKshell(a.MaxKshell),
  OneOverN(a.OneOverN),Kshell(a.Kshell),Kmag(a.Kmag),
  OneOverDnk(a.OneOverDnk),RhokTot(a.RhokTot),SkInst(a.SkInst)
{
  nList.push_back(Title);
  dList.push_back(new VectorEstimatorType(Title,NumK));
}

CompositeEstimatorBase* SkEstimator::clone()
{
  return new SkEstimator(*this);
}

SkEstimator::~SkEstimator()
{
}

void SkEstimator::resetTargetParticleSet(ParticleSet& p)
{
  Title="sk_"+p.getName();
}

// void SkEstimator::open(hid_t hroot)
// {
//   if(GroupID<0)
//   {
//     GroupID = H5Gcreate(hroot,Title.c_str(),0);
//     Sk_h = new HDFAttribIO<VectorEstimatorType>(Sk);
//     Sk_h->reserve(GroupID);
//   }
// }

// void SkEstimator::close()
// {
//   if(GroupID>-1)
//   {
//     delete Sk_h;
//     Sk_h=0;
//     H5Gclose(GroupID);
//     GroupID=-1;
//   }
// }

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
    accumulate_elements(p.SK->rhok[i],p.SK->rhok[i]+NumK,RhokTot.begin());
  //accumulate_elements(p.SK->rhok[i],RhokTot.begin(),NumK);
  Vector<ComplexType>::const_iterator iit(RhokTot.begin());
  Vector<RealType>::iterator oit(SkInst.begin());
  for(int k=0; k<NumK; k++,iit++,oit++)
  {
    (*oit)+=(*iit).real()*(*iit).real()+(*iit).imag()*(*iit).imag();
  }
}

/** add gofrInst which contains sum over walkers */
void SkEstimator::stopAccumulate(RealType wgtnorm)
{
  dList[0]->accumulate(SkInst.begin(),wgtnorm*OneOverN);
}

hid_t SkEstimator::createGroup(hid_t gid)
{
  hid_t h1 = H5Gcreate(gid,this->Title.c_str(),0);
  //save k vector and kshell
  return h1;
}

// void SkEstimator::writeHeaders(hid_t gid)
// {
//   hid_t h1 = H5Gcreate(gid,"sk",0);
//   //save k vector and kshell
//   H5Gclose(h1);
// }
}

/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1415 $   $Date: 2006-10-23 11:51:53 -0500 (Mon, 23 Oct 2006) $
 * $Id: CompositeEstimatorBase.h 1415 2006-10-23 16:51:53Z jnkim $
 ***************************************************************************/

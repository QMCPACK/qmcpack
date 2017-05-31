//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Bryan Clark, bclark@Princeton.edu, Princeton University
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Bryan Clark, bclark@Princeton.edu, Princeton University
//////////////////////////////////////////////////////////////////////////////////////
    
    
    
    
    
    
    
    
#include "Estimators/GofREstimator.h"
#include "Particle/DistanceTable.h"
#include "Particle/DistanceTableData.h"
#include "Utilities/IteratorUtility.h"
#include "Numerics/HDFNumericAttrib.h"
#include "OhmmsData/HDFStringAttrib.h"
//#define PRINT_DEBUG_GOFR

namespace qmcplusplus
{

GofREstimator::GofREstimator(ParticleSet& source):
  Symmetric(true),sourcePtcl(&source),targetPtcl(0)
{
  myTable = DistanceTable::add(source);
  int ns=sourcePtcl->groups();
  nList.clear();
  std::vector<int> mask(ns*ns,-1);
  int ij=0;
  for(int i=0; i<ns; i++)
    for(int j=i; j<ns; j++,ij++)
    {
      mask[j+i*ns]=ij;
      char fname[32];
      sprintf(fname,"%s_%d_%d",myTable->Name.c_str(),i,j);
      nList.push_back(fname);
    }
  NumPairTypes=ij;
  Centers=sourcePtcl->getTotalNum();
  PairID.resize(myTable->getTotNadj());
  for(int iat=0; iat<Centers; iat++)
  {
    for(int nn=myTable->M[iat]; nn<myTable->M[iat+1]; nn++)
    {
      PairID[nn]=mask[myTable->PairID[nn]];
    }
  }
  setBound(0.1);
}

GofREstimator::GofREstimator(ParticleSet& source,ParticleSet& target):
  Symmetric(false),sourcePtcl(&source),targetPtcl(&target)
{
  myTable = DistanceTable::add(source,target);
  NumPairTypes=sourcePtcl->getSpeciesSet().getTotalNum();
  for(int i=0; i<NumPairTypes; i++)
  {
    char fname[32];
    sprintf(fname,"%s_%s",myTable->Name.c_str(),
            sourcePtcl->getSpeciesSet().speciesName[i].c_str());
    nList.push_back(fname);
  }
  Centers=sourcePtcl->getTotalNum();
  PairID.resize(myTable->getTotNadj());
  for(int iat=0; iat<Centers; iat++)
  {
    for(int nn=myTable->M[iat]; nn<myTable->M[iat+1]; nn++)
    {
      PairID[nn]=sourcePtcl->GroupID[iat];
    }
  }
  setBound(0.1);
}

GofREstimator::~GofREstimator()
{
}

CompositeEstimatorBase* GofREstimator::clone()
{
  CompositeEstimatorBase* est=0;
  if(Symmetric)
    est= new GofREstimator(*sourcePtcl);
  else
    est= new GofREstimator(*sourcePtcl,*targetPtcl);
  est->Title=this->Title;
  return est;
}


void GofREstimator::resetTargetParticleSet(ParticleSet& p)
{
  if(Symmetric)
    myTable=DistanceTable::add(p);
  else
    myTable=DistanceTable::add(*sourcePtcl,p);
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
  for(int iat=0; iat<Centers; iat++)
  {
    for(int nn=myTable->M[iat]; nn<myTable->M[iat+1]; nn++)
    {
      if(myTable->r(nn)>=Dmax)
        continue;
      //need a better average
      int ig=static_cast<int>(DeltaInv*myTable->r(nn));
      gofrInst(PairID[nn],ig)+=1;
    }
  }
}

/** add gofrInst which contains sum over walkers */
void GofREstimator::stopAccumulate(RealType wgtnorm)
{
  for(int p=0; p<NumPairTypes; p++)
  {
    dList[p]->accumulate(gofrInst[p],gofrInst[p]+NumBins,normFactor.begin(),wgtnorm);
  }
}

void GofREstimator::setBound(RealType dr)
{
  RealType vnorm=1.0;
  if(sourcePtcl->Lattice.SuperCellEnum)
  {
    Dmax=sourcePtcl->Lattice.LR_rc;
    /** normalizaton factor */
    vnorm=4.0*M_PI*myTable->size(DistanceTableData::SourceIndex)*myTable->size(DistanceTableData::VisitorIndex);
    vnorm=sourcePtcl->Lattice.Volume/vnorm;
  }
  else
  {
    Dmax=10.0; //choose a sensible number
  }
  //Dmax=rmax;
  Delta=dr;
  DeltaInv=1.0/dr;
  NumBins=static_cast<int>((Dmax)*DeltaInv+1);
  normFactor.resize(NumBins);
  RealType r=Delta;
  for(int i=1; i<NumBins; i++, r+=Delta)
    normFactor[i]=vnorm/(r*r);
  gofrInst.resize(NumPairTypes,NumBins);
  //clean up the data before using
  delete_iter(dList.begin(),dList.end());
  for(int p=0; p<NumPairTypes; p++)
  {
    dList.push_back(new VectorEstimatorType(nList[p],NumBins));
  }
}

hid_t GofREstimator::createGroup(hid_t gid)
{
  hid_t h1=H5Gcreate(gid,this->Title.c_str(),0);
  Vector<RealType> rv(NumBins);
  RealType r=0;
  for(int i=0; i<NumBins; i++,r+=Delta)
    rv[i]=r;
  HDFAttribIO<Vector<RealType> > ro(rv);
  ro.write(h1,"distance_table");
  return h1;
}

//  void GofREstimator::writeHeaders(hid_t gid)
//  {
//    //hid_t h1 = H5Gcreate(gid,"gofr",0);
//    Vector<RealType> rv(NumBins);
//    RealType r=0;
//    for(int i=0; i<NumBins;i++,r+=Delta) rv[i]=r;
//    HDFAttribIO<Vector<RealType> > ro(rv);
//    ro.write(gid,"distances");
//    std::ostringstream o;
//    for(int i=0; i<NumPairTypes-1; i++) o << nList[i] <<":";
//    o<<nList.back();
//    std::string banner(o.str());
//    HDFAttribIO<std::string> so(banner);
//    so.write(h1,"pair_ids");
//
//    H5Gclose(h1);
//  }
}


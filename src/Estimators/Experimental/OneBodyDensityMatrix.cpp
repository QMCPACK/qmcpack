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
    
    
    
    
    
    
    
    



#include "Estimators/OneBodyDensityMatrix.h"
#include "Particle/DistanceTable.h"
#include "Particle/DistanceTableData.h"
#include "Utilities/IteratorUtility.h"
#include "Numerics/HDFNumericAttrib.h"
#include "OhmmsData/HDFStringAttrib.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "ParticleBase/RandomSeqGenerator.h"
//#define PRINT_DEBUG_GOFR

namespace qmcplusplus
{
nofrEstimator::nofrEstimator(TrialWaveFunction &psi,MCWalkerConfiguration& w):
  Symmetric(true),sourcePtcl(0),targetPtcl(0), Psi(psi), NumSamples(0),W(w)
{
  nList.clear();
  nList.push_back("nofrvals");
  setBound(0.1);
//     myTable = DistanceTable::add(source);
//     int ns=sourcePtcl->groups();
//
//     std::vector<int> mask(ns*ns,-1);
//     int ij=0;
//     for(int i=0; i<ns; i++)
//       for(int j=i; j<ns; j++,ij++)
//       {
//         mask[j+i*ns]=ij;
//         char fname[32];
//         sprintf(fname,"%s_%d_%d",myTable->Name.c_str(),i,j);
//         nList.push_back(fname);
//       }
//     NumPairTypes=ij;
//     Centers=sourcePtcl->getTotalNum();
//     PairID.resize(myTable->getTotNadj());
//     for(int iat=0; iat<Centers; iat++) {
//       for(int nn=myTable->M[iat]; nn<myTable->M[iat+1]; nn++)
//       {
//         PairID[nn]=mask[myTable->PairID[nn]];
//       }
//     }
//     setBound(0.1);
}

//  nofrEstimator(const nofrEstimator& pc,TrialWaveFunction &psi): sourcePtcl(pc.sourcePtcl),Psi(psi) {;}

//   nofrEstimator::nofrEstimator(ParticleSet& source,ParticleSet& target,TrialWaveFunction &psi):
//     Symmetric(false),sourcePtcl(&source),targetPtcl(&target),Psi(psi)
//   {
//     myTable = DistanceTable::add(source,target);
//     NumPairTypes=sourcePtcl->getSpeciesSet().getTotalNum();
//     for(int i=0; i<NumPairTypes; i++)
//     {
//       char fname[32];
//       sprintf(fname,"%s_%s",myTable->Name.c_str(),
//           sourcePtcl->getSpeciesSet().speciesName[i].c_str());
//       nList.push_back(fname);
//     }
//     Centers=sourcePtcl->getTotalNum();
//     PairID.resize(myTable->getTotNadj());
//     for(int iat=0; iat<Centers; iat++) {
//       for(int nn=myTable->M[iat]; nn<myTable->M[iat+1]; nn++)
//       {
//         PairID[nn]=sourcePtcl->GroupID[iat];
//       }
//     }

//     setBound(0.1);
//   }

nofrEstimator::~nofrEstimator()
{
}

CompositeEstimatorBase* nofrEstimator::clone()
{
  CompositeEstimatorBase* est=0;
  est= new nofrEstimator(Psi,W);
  est->Title=this->Title;
  return est;
}


void nofrEstimator::resetTargetParticleSet(ParticleSet& p)
{
  if(Symmetric)
    myTable=DistanceTable::add(p);
  else
    myTable=DistanceTable::add(*sourcePtcl,p);
}

/** ready to accumulate the measurements over the walkers
 */
void nofrEstimator::startAccumulate()
{
  gofrInst=0.0;
}

void
nofrEstimator::PutInBox (PosType &v)
{
  for (int i=0; i<DIM; i++)
  {
    double n = -std::floor(v(i)*(1.0/W.Lattice.R(i,i))+0.5);
    v(i) += n*W.Lattice.R(i,i);
  }
}


/** accumulate the observables for a walker image*/
void nofrEstimator::accumulate(ParticleSet& p)
{
  NumSamples++;
  std::cerr <<"Accumulating nofr"<< std::endl;
  ParticleSet::ParticlePos_t gRand;
  gRand.resize(1);
  ///Here we are going to compute the ODLRO stuff
  for (int ptcl=0; ptcl<p.R.size(); ptcl++)
  {
    double dist;
    //      do {
    makeUniformRandom(gRand);
    //      makeGaussRandom(gRand);
    ///      gRand=Box[0]*gRand;
    do
    {
      for (int dim=0; dim<gRand[0].size(); dim++)
        gRand[0][dim]=gRand[0][dim]*W.Lattice.R(dim,dim)/2.0;
      PutInBox(gRand[0]);
      dist=std::sqrt(dot(gRand[0],gRand[0]));
    }
    while (dist>W.Lattice.R(0,0)/2.0);
    p.makeMove(ptcl,gRand[0]);
    RealType ratio = Psi.ratio(p,ptcl);
    //      ratio=ratio*ratio;
    p.rejectMove(ptcl);
    //      std::cerr <<"Ratio is "<<dist<<" "<<ratio<< std::endl;
    if (!(dist>=Dmax))
    {
      int ig=static_cast<int>(DeltaInv*dist);
      gofrInst(0,ig)+=ratio;
    }
  }
}






//     for (int ptcl=0;ptcl<p.R.size();ptcl++){
//       //      for (int count=0;count<p.R.size();count++){
// 	makeGaussRandom(gRand);

// 	PosType oldPos=p.R[ptcl];

// 	gRand=0.1*gRand;
// 	//	p.R[ptcl]=p.R[ptcl]+gRand[0];
// 	/////	PutInBox(p.R[ptcl]);

// 	PosType diff=p.R[ptcl]-oldPos;
// 	//	std::cerr <<"RAND: "<<gRand[0]<<" "<<diff<<" "<< std::endl;
// 	PutInBox(diff);
// 	//	std::cerr <<"post-diff: "<<diff<< std::endl;
// 	//      p.Lattice.inBox(diff);//applyBC(diff);
// 	double dist=std::sqrt(dot(diff,diff));
// 	//	std::cerr <<"dist: "<<dist<< std::endl;
// 	RealType ratio = Psi.ratio(p,ptcl);
// 	p.R[ptcl]=oldPos;
// 	std::cerr <<dist<<" "<<ratio<< std::endl;
// 	if (!(dist>=Dmax)){
// 	  int ig=static_cast<int>(DeltaInv*dist);
// 	  gofrInst(0,ig)+=ratio*ratio;
// 	}
// 	//      }
//     }
//     std::cerr <<"Done accumulating nofr"<< std::endl;

//     for(int iat=0; iat<Centers; iat++) {
//       for(int nn=myTable->M[iat]; nn<myTable->M[iat+1]; nn++)
//       {
//         if(myTable->r(nn)>=Dmax) continue;
//         //need a better average
//         int ig=static_cast<int>(DeltaInv*myTable->r(nn));
//         gofrInst(PairID[nn],ig)+=1;
//       }
//     }

//  }

/** add gofrInst which contains sum over walkers */
void nofrEstimator::stopAccumulate(RealType wgtnorm)
{
  for(int p=0; p<NumPairTypes; p++)
  {
    dList[p]->accumulate(gofrInst[p],gofrInst[p]+NumBins,normFactor.begin(),wgtnorm);
  }
}

void nofrEstimator::setBound(RealType dr)
{
  ////    RealType vnorm=1.0;
  RealType vnorm=1.0/(double)NumSamples;
  NumSamples=0;
//     if(sourcePtcl->Lattice.SuperCellEnum)
//     {
//       Dmax=sourcePtcl->Lattice.LR_rc;
//       /** normalizaton factor */
//       //vnorm=4.0*M_PI*myTable->size(DistanceTableData::SourceIndex)*myTable->size(DistanceTableData::VisitorIndex);
//       //      vnorm=sourcePtcl->Lattice.Volume/vnorm;
//       vnorm=1.0;
//     }
//     else
//     {
//       Dmax=10.0; //choose a sensible number
//     }
  vnorm=1.0;
  Dmax=10.0;
  NumPairTypes=1;
  //Dmax=rmax;
  Delta=dr;
  DeltaInv=1.0/dr;
  NumBins=static_cast<int>((Dmax)*DeltaInv+1);
  normFactor.resize(NumBins);
  RealType r=Delta;
  for(int i=1; i<NumBins; i++, r+=Delta)
    normFactor[i]=(vnorm*Delta)/
                  ((r+Delta/2)*(r+Delta/2)*(r+Delta/2)-((r-Delta/2)*(r-Delta/2)*(r-Delta/2)));  //normFactor[i]=vnorm/(r*r);
  gofrInst.resize(NumPairTypes,NumBins);
  //clean up the data before using
  delete_iter(dList.begin(),dList.end());
  for(int p=0; p<NumPairTypes; p++)
  {
    dList.push_back(new VectorEstimatorType(nList[p],NumBins));
  }
}

hid_t nofrEstimator::createGroup(hid_t gid)
{
  std::cout <<"nofr group creating"<< std::endl;
  hid_t h1=H5Gcreate(gid,this->Title.c_str(),0);
  Vector<RealType> rv(NumBins);
  RealType r=0;
  for(int i=0; i<NumBins; i++,r+=Delta)
    rv[i]=r;
  HDFAttribIO<Vector<RealType> > ro(rv);
  ro.write(h1,"distance_table");
  return h1;
}

//  void nofrEstimator::writeHeaders(hid_t gid)
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


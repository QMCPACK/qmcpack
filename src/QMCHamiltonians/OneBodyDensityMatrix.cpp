//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim and Bryan Clark
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
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
#include "QMCHamiltonians/OneBodyDensityMatrix.h"

namespace qmcplusplus
{

OneBodyDensityMatrix::OneBodyDensityMatrix(TrialWaveFunction&  psi
    , ParticleSet& ions, ParticleSet& elns )
  :Psi(psi)
{
}

OneBodyDensityMatrix::OneBodyDensityMatrix(TrialWaveFunction&  psi
    , ParticleSet& elns)
  :Psi(psi)
{
}

void OneBodyDensityMatrix::resetTargetParticleSet(ParticleSet& P)
{
}

OneBodyDensityMatrix::Return_t OneBodyDensityMatrix::evaluate(ParticleSet& P)
{
  //ParticleSet::ParticlePos_t gRand;
  //gRand.resize(1);
  /////Here we are going to compute the ODLRO stuff
  //for (int ptcl=0;ptcl<p.R.size();ptcl++){
  //  double dist;
  //  //      do {
  //  makeUniformRandom(gRand);
  //  //      makeGaussRandom(gRand);
  //  ///      gRand=Box[0]*gRand;
  //  do {
  //    for (int dim=0;dim<gRand[0].size();dim++)
  //      gRand[0][dim]=gRand[0][dim]*W.Lattice.R(dim,dim)/2.0;
  //    PutInBox(gRand[0]);
  //    dist=std::sqrt(dot(gRand[0],gRand[0]));
  //  } while (dist>W.Lattice.R(0,0)/2.0);
  //  p.makeMove(ptcl,gRand[0]);
  //  RealType ratio = Psi.ratio(p,ptcl);
  //  //      ratio=ratio*ratio;
  //  p.rejectMove(ptcl);
  //  //      cerr<<"Ratio is "<<dist<<" "<<ratio<<endl;
  //  if (dist<Dmax)
  //  {
  //    gofr[static_cast<int>(DeltaInv*dist)]+=ratio;
  //  }
  //}
  return 0.0;
}

void OneBodyDensityMatrix::addObservables(PropertySetType& plist)
{
  for(int i=0; i<gofr.size(); ++i)
  {
    ostringstream obsName;
    osbName<<myName<<"_"<<i;
    myIndex=plist.add(obsName.str());
  }
  myIndex-=(gfor.size()-1);
}

void OneBodyDensityMatrix::setObservables(PropertySetType& plist)
{
  std::copy(gofr.begin(),gofr.end(),plist.begin()+myIndex);
}

void OneBodyDensityMatrix::setParticlePropertyList(PropertySetType& plist
    , int offset)
{
  std::copy(gofr.begin(),gofr.end(),plist.begin()+offset);
}

bool OneBodyDensityMatrix::put(xmlNodePtr cur)
{
  ///process attributes
  resize();
}

bool OneBodyDensityMatrix::get(std::ostream& os) const
{
  os << myName << " dmax=" << Dmax << endl;
}

QMCHamiltonianBase* OneBodyDensityMatrix::makeClone(ParticleSet& qp
    , TrialWaveFunction& psi)
{
  OneBodyDensityMatrix* myclone=new OneBodyDensityMatrix(psi,qp);
  myclone->resize();
  return myclone;
}

void  OneBodyDensityMatrix::resize()
{
  //vnorm=1.0;
  //Dmax=10.0;
  //NumPairTypes=1;
  ////Dmax=rmax;
  //Delta=dr;
  //DeltaInv=1.0/dr;
  //NumBins=static_cast<int>((Dmax)*DeltaInv+1);
  //normFactor.resize(NumBins);
  //RealType r=Delta;
  //for(int i=1; i<NumBins; i++, r+=Delta)
  //  normFactor[i]=(vnorm*Delta)/
  //    ((r+Delta/2)*(r+Delta/2)*(r+Delta/2)-((r-Delta/2)*(r-Delta/2)*(r-Delta/2)));  //normFactor[i]=vnorm/(r*r);
  //gofrInst.resize(NumPairTypes,NumBins);
}
}

/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 2945 $   $Date: 2008-08-05 10:21:33 -0500 (Tue, 05 Aug 2008) $
 * $Id: ForceBase.h 2945 2008-08-05 15:21:33Z jnkim $
 ***************************************************************************/

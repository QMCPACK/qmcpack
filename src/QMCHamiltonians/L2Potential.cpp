//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
    
    
#include "Particle/ParticleSet.h"
#include "Particle/DistanceTableData.h"
#include "QMCHamiltonians/QMCHamiltonianBase.h"
#include "QMCHamiltonians/L2Potential.h"
#include "Utilities/IteratorUtility.h"
#include "QMCWaveFunctions/OrbitalBase.h"

namespace qmcplusplus
{

typedef OrbitalBase::ValueType ValueType;
typedef OrbitalBase::HessVector_t HessVector_t;


L2Potential::L2Potential(const ParticleSet& ions, ParticleSet& els, TrialWaveFunction& psi):
  IonConfig(ions)
{
  set_energy_domain(potential);
  two_body_quantum_domain(ions,els);
  NumIons=ions.getTotalNum();
  myTableIndex=els.addTable(ions,DT_SOA_PREFERRED);
  size_t ns = ions.getSpeciesSet().getTotalNum();
  PPset.resize(ns,0);
  PP.resize(NumIons,nullptr);
  psi_ref = &psi;
}


///destructor
L2Potential::~L2Potential()
{
  delete_iter(PPset.begin(),PPset.end());
}


void L2Potential::resetTargetParticleSet(ParticleSet& P)
{
  int tid=P.addTable(IonConfig,DT_SOA_PREFERRED);
  if(tid != myTableIndex)
  {
    APP_ABORT("  L2Potential::resetTargetParticleSet found a different distance table index.");
  }
}


void L2Potential::add(int groupID, L2RadialPotential* ppot)
{
  PPset[groupID]=ppot;
  for(int iat=0; iat<PP.size(); iat++)
    if(IonConfig.GroupID[iat]==groupID)
      PP[iat]=ppot;
}


L2Potential::Return_t
L2Potential::evaluate(ParticleSet& P)
{
  // compute the Hessian
  HessVector_t D2;
  // evaluateHessian gives the Hessian of log(Psi)
  psi_ref->evaluateHessian(P,D2);
  // add gradient terms to get (Hessian(Psi))/Psi instead
  size_t N = P.getTotalNum();
  for(size_t n=0;n<N;n++)
    for(size_t i=0;i<DIM;i++)
      for(size_t j=0;j<DIM;j++)
        D2[n](i,j) += P.G[n][i]*P.G[n][j];
  // zero out diagonal elements for efficient L2 use
  for(int n=0;n<N;n++)
    for(int d=0;d<DIM;d++)
      D2[n](d,d)=0.0;

  const DistanceTableData& d_table(*P.DistTables[myTableIndex]);
  Value=0.0;
  if(d_table.DTType==DT_SOA)
  {
    const size_t Nelec = P.getTotalNum();
    for(size_t iel=0; iel<Nelec; ++iel)
    {
      const RealType* restrict dist=d_table.Distances[iel];
      Return_t esum = 0.0;
      for(size_t iat=0; iat<NumIons; ++iat)
        if(PP[iat]!=nullptr) esum+=PP[iat]->evaluate(dist[iat]);
      Value += esum;
    }
  }
  else
  {
    for(int iat=0; iat<NumIons; iat++)
    {
      L2RadialPotential* ppot = PP[iat];
      if(ppot==nullptr) continue;
      Return_t esum = 0.0;
      for(int nn=d_table.M[iat]; nn<d_table.M[iat+1]; ++nn)
        esum += ppot->evaluate(d_table.r(nn));
      Value += esum;
    }
  }
  return Value;
}


QMCHamiltonianBase* L2Potential::makeClone(ParticleSet& qp, TrialWaveFunction& psi)
{
  L2Potential* myclone=new L2Potential(IonConfig,qp,psi);
  for(int ig=0; ig<PPset.size(); ++ig)
    if(PPset[ig])
    {
      L2RadialPotential* ppot=PPset[ig]->makeClone();
      myclone->add(ig,ppot);
    }
  return myclone;
}

}


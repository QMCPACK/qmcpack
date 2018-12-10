//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#include "QMCHamiltonians/ZeroVarianceForce.h"
#include "Particle/DistanceTable.h"
#include "Particle/DistanceTableData.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "ParticleBase/ParticleAttribOps.h"

namespace qmcplusplus
{

ZeroVarianceForce::ZeroVarianceForce(ParticleSet& ions, ParticleSet& elns,
                                     TrialWaveFunction &psi)
  : ForceBase(ions, elns), Ions(ions), Electrons(elns),
    Psi(psi)
{
  for (int dim=0; dim<OHMMS_DIM; dim++)
  {
    grad_grad_psi[dim].resize(Nel);
    lapl_grad_psi[dim].resize(Nel);
  }
  F_ZV1.resize(Nnuc);
  F_ZV2.resize(Nnuc);
}

void
ZeroVarianceForce::resetTargetParticleSet(ParticleSet& P)
{
  int tid=P.addTable(Ions,DT_AOS);
  if(tid != myTableIndex)
    APP_ABORT("ZeroVarianceForce::resetTargetParticleSet found inconsistent table index");
}

void ZeroVarianceForce::addObservables(PropertySetType& plist
                                       , BufferType& collectables)
{
  if(FirstForceIndex<0)
    FirstForceIndex=plist.size();
    
 //   app_log() << "ionionforce = "<<ionionforce<< std::endl;
  app_log() << "addionion="<<addionion<< std::endl;
  app_log() << "FirstTime= "<<FirstTime<< std::endl;
  for(int iat=0; iat<Nnuc; iat++)
  {
    for(int x=0; x<OHMMS_DIM; x++)
    {
      std::ostringstream obsName1, obsName2;
      obsName1 << "F_ZV1" << "_" << iat << "_" << x;
      obsName2 << "F_ZV2" << "_" << iat << "_" << x;
      app_log() << "Adding " << obsName1.str() << " to observable list.\n";
      app_log() << "Adding " << obsName2.str() << " to observable list.\n";
      plist.add(obsName1.str());
      plist.add(obsName2.str());
    }
  }
}

void
ZeroVarianceForce::registerObservables(std::vector<observable_helper*>& h5list,
                                       hid_t gid) const
{
  QMCHamiltonianBase::registerObservables(h5list, gid);
  std::vector<int> ndim(2);
  ndim[0]=Nnuc;
  ndim[1]=OHMMS_DIM;
  observable_helper* h5o1 = new observable_helper("F_ZV1");
  observable_helper* h5o2 = new observable_helper("F_ZV2");
  h5o1->set_dimensions(ndim,FirstForceIndex);
  h5o1->open(gid);
  h5list.push_back(h5o1);
  h5o2->set_dimensions(ndim,FirstForceIndex);
  h5o2->open(gid);
  h5list.push_back(h5o2);
}

void
ZeroVarianceForce::setObservables(QMCTraits::PropertySetType& plist)
{
  //QMCHamiltonianBase::setObservables(plist);
  setObservablesF(plist);
  int index = FirstForceIndex;
  for(int iat=0; iat<Nnuc; iat++)
  {
    for(int x=0; x<OHMMS_DIM; x++)
    {
      double ZV = F_ZV1[iat][x] + F_ZV2[iat][x];
      // ZV = min (ZV, 50.0);
      // ZV = max (ZV, -50.0);
      plist[index++] = F_ZV1[iat][x];
      plist[index++] = F_ZV2[iat][x];
      //plist[index++] = ZV - F_ZV1[iat][x];
    }
  }
}


void
ZeroVarianceForce::setParticlePropertyList
(QMCTraits::PropertySetType& plist, int offset)
{
  QMCHamiltonianBase::setParticlePropertyList (plist, offset);
  int index = FirstForceIndex + offset;
  for(int iat=0; iat<Nnuc; iat++)
    for(int x=0; x<OHMMS_DIM; x++)
    {
      double ZV = F_ZV1[iat][x] + F_ZV2[iat][x];
      // ZV = min (ZV, 50.0);
      // ZV = max (ZV, -50.0);
      plist[index++] = F_ZV1[iat][x];
      plist[index++] = F_ZV2[iat][x];
    //  plist[index++] = ZV - F_ZV1[iat][x];
	//	plist[index++]= forces[iat][x];
	//	plist[index++]=ZV;
    }
}

ZeroVarianceForce::Return_t
ZeroVarianceForce::evaluate(ParticleSet& P)
{
  forces = forces_IonIon;
  const DistanceTableData* d_ab=P.DistTables[myTableIndex];
  const ParticleScalar_t* restrict Zat=Ions.Z.first_address();
  const ParticleScalar_t* restrict Qat=P.Z.first_address();
  //Loop over distinct eln-ion pairs
  for(int iat=0; iat<Nnuc; iat++)
  {
    for(int nn=d_ab->M[iat], jat=0; nn<d_ab->M[iat+1]; nn++,jat++)
    {
      real_type rinv=d_ab->rinv(nn);
      real_type r3zz=Qat[jat]*Zat[iat]*rinv*rinv*rinv;
      forces[iat] -= r3zz*d_ab->dr(nn);
    }
    
    F_ZV1[iat]=forces[iat];
  }
  
  for (int ion=0; ion < Nnuc; ion++)
  {
    GradType grad = Psi.evalGradSource(P, Ions, ion, grad_grad_psi, lapl_grad_psi);
    for (int dim=0; dim < OHMMS_DIM; dim++)
    {
    //  F_ZV1[ion][dim] = -0.5*(Sum(lapl_grad_psi[dim]));
      //F_ZV1[ion][dim]=0.0; //forces[ion][dim];
     // F_ZV2[ion][dim] = -Dot(grad_grad_psi[dim], P.G);
		double ZV = -0.5*(Sum(lapl_grad_psi[dim]))-Dot(grad_grad_psi[dim], P.G);
		F_ZV1[ion][dim]-=ZV;
		F_ZV2[ion][dim]=ZV;
    }
  }
  return Value=0.0;
}
}

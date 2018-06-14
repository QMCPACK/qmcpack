//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#include "QMCHamiltonians/PulayForce.h"
#include "Particle/DistanceTable.h"
#include "Particle/DistanceTableData.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"

namespace qmcplusplus
{

PulayForce::PulayForce(ParticleSet& ions, ParticleSet& elns,
                       TrialWaveFunction &psi)
  : ForceBase(ions, elns), Ions(ions), Electrons(elns),
    Psi(psi)
{
  GradLogPsi.resize(Nnuc);
  EGradLogPsi.resize(Nnuc);
  WarpNorm.resize(Nel, 0.0);
}

void
PulayForce::resetTargetParticleSet(ParticleSet& P)
{
  int tid=P.addTable(Ions,DT_AOS);
  if(tid != myTableIndex)
    APP_ABORT("PulayForce::resetTargetParticleSet found inconsistent table index");
}

void PulayForce::addObservables(PropertySetType& plist
                                , BufferType& collectables)
{
  if(FirstForceIndex<0)
    FirstForceIndex=plist.size();
  for(int iat=0; iat<Nnuc; iat++)
  {
    for(int x=0; x<OHMMS_DIM; x++)
    {
      std::ostringstream obsName1, obsName2;
      obsName1 << "grad_log_psi" << "_" << iat << "_" << x;
      plist.add(obsName1.str());
      obsName2 << "Egrad_log_psi" << "_" << iat << "_" << x;
      plist.add(obsName2.str());
    }
  }
}

void
PulayForce::registerObservables(std::vector<observable_helper*>& h5list,
                                hid_t gid) const
{
  QMCHamiltonianBase::registerObservables(h5list, gid);
  std::vector<int> ndim(2);
  ndim[0]=Nnuc;
  ndim[1]=OHMMS_DIM;
  observable_helper* h5o1 = new observable_helper("grad_log_psi");
  h5o1->set_dimensions(ndim,FirstForceIndex);
  h5o1->open(gid);
  h5list.push_back(h5o1);
  observable_helper* h5o2 = new observable_helper("Egrad_log_psi");
  h5o2->set_dimensions(ndim,FirstForceIndex+Nnuc*OHMMS_DIM);
  h5o2->open(gid);
  h5list.push_back(h5o2);
}

void
PulayForce::setObservables(QMCTraits::PropertySetType& plist)
{
  QMCHamiltonianBase::setObservables(plist);
  int index = FirstForceIndex;
  for(int iat=0; iat<Nnuc; iat++)
  {
    for(int x=0; x<OHMMS_DIM; x++)
    {
      plist[index++] = GradLogPsi[iat][x];
      plist[index++] = EGradLogPsi[iat][x];
    }
  }
}


void
PulayForce::setParticlePropertyList(QMCTraits::PropertySetType& plist,
                                    int offset)
{
  QMCHamiltonianBase::setParticlePropertyList (plist, offset);
  int index = FirstForceIndex + offset;
  for(int iat=0; iat<Nnuc; iat++)
  {
    for(int x=0; x<OHMMS_DIM; x++)
    {
      plist[index++] = GradLogPsi[iat][x];
      plist[index++] = EGradLogPsi[iat][x];
    }
  }
}

PulayForce::Return_t
PulayForce::evaluate(ParticleSet& P)
{
  const DistanceTableData& d_ab = *P.DistTables[myTableIndex];
  for (int ion=0; ion < Nnuc; ion++)
  {
    GradLogPsi[ion] = Psi.evalGradSource(P, Ions, ion);
    //      RealType E = tWalker->Properties(0,LOCALENERGY);
    RealType E = P.PropertyList[LOCALENERGY];
    EGradLogPsi[ion] = E * GradLogPsi[ion];
  }
  return Value=0.0;
  //    std::cerr << "In PulayForce::evaluate(ParticleSet& P).\n";
  // Compute normalization of the warp tranform for each electron
  for (int elec=0; elec<Nel; elec++)
    WarpNorm[elec] = 0.0;
  for (int ion=0; ion<Nnuc; ion++)
    for(int nn=d_ab.M[ion], elec=0; nn<d_ab.M[ion+1]; ++nn,++elec)
      WarpNorm[elec] += WarpFunction(d_ab.r(nn));
  for (int elec=0; elec<Nel; elec++)
    WarpNorm[elec] = 1.0/WarpNorm[elec];
  // Now, compute the approximation to the gradient of psi
  // w.r.t. the ions by using the gradient w.r.t. the electrons.
  for (int ion=0; ion < Nnuc; ion++)
  {
    GradLogPsi[ion] = EGradLogPsi[ion] = PosType();
    for(int nn=d_ab.M[ion], elec=0; nn<d_ab.M[ion+1]; ++nn,++elec)
      GradLogPsi[ion] -= P.G[elec] * static_cast<ParticleSet::ParticleValue_t>(WarpNorm[elec] * WarpFunction(d_ab.r(nn)));
      //GradLogPsi[ion] -= P.G[elec] * (ParticleSet::Scalar_t)(WarpNorm[elec] * WarpFunction(d_ab.r(nn)));
    RealType E = tWalker->Properties(0,LOCALENERGY);
    // EGradLogPsi[ion] = P.getPropertyBase()[LOCALENERGY] * GradLogPsi[ion];
    EGradLogPsi[ion] = E * GradLogPsi[ion];
  }
  //    std::cerr << "Finish PulayForce::evaluate(ParticleSet& P).\n";
  return Value = 0.0;
}


}

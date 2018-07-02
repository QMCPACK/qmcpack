//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    

#include "QMCHamiltonians/CoulombPBCAB_CUDA.h"
#include "Particle/MCWalkerConfiguration.h"

namespace qmcplusplus
{

CoulombPBCAB_CUDA::CoulombPBCAB_CUDA
(ParticleSet& ions, ParticleSet& elns, bool cloning) :
  CoulombPBCAB(ions,elns,cloning),
  ElecRef(elns), IonRef(ions),
  SumGPU("CoulombPBCABTemp::SumGPU"),
  IGPU("CoulombPBCABTemp::IGPU"),
  L("CoulombPBCABTemp::L"),
  Linv("CoulombPBCABTemp::Linv"),
  kpointsGPU("CoulombPBCABTemp::kpointsGPU"),
  kshellGPU("CoulombPBCABTemp::kshellGPU"),
  FkGPU("CoulombPBCABTemp::FkGPU"),
  RhoklistGPU("CoulombPBCABTemp::RhoklistGPU"),
  RhokElecGPU("CoulombPBCABTemp::RhokElecGPU")
{
  MaxGridPoints = 8191;
  SpeciesSet &sSet = ions.getSpeciesSet();
  NumIonSpecies = sSet.getTotalNum();
  NumIons  = ions.getTotalNum();
  NumElecs = elns.getTotalNum();
#ifdef QMC_CUDA
  gpu::host_vector<CUDA_PRECISION_FULL> LHost(9), LinvHost(9);
  for (int i=0; i<3; i++)
    for (int j=0; j<3; j++)
    {
      LHost[3*i+j]    = elns.Lattice.a(j)[i];
      LinvHost[3*i+j] = elns.Lattice.b(i)[j];
    }
  L    = LHost;
  Linv = LinvHost;
  // Copy center positions to GPU, sorting by GroupID
  gpu::host_vector<CUDA_PRECISION> I_host(OHMMS_DIM*NumIons);
  int index=0;
  for (int cgroup=0; cgroup<NumIonSpecies; cgroup++)
  {
    IonFirst.push_back(index);
    for (int i=0; i<NumIons; i++)
    {
      if (ions.GroupID[i] == cgroup)
      {
        for (int dim=0; dim<OHMMS_DIM; dim++)
          I_host[OHMMS_DIM*index+dim] = ions.R[i][dim];
        SortedIons.push_back(ions.R[i]);
        index++;
      }
    }
    IonLast.push_back(index-1);
  }
  IGPU = I_host;
  initBreakup(elns);
  setupLongRangeGPU();
#endif
}

void
CoulombPBCAB_CUDA::initBreakup(ParticleSet& P)
{
  CoulombPBCAB::initBreakup(P);
#ifdef QMC_CUDA
  V0Spline = new TextureSpline;
  V0Spline->set(V0->data(), V0->size(), V0->grid().rmin(),
                V0->grid().rmax());
  SRSplines.resize(NumIonSpecies, V0Spline);
#endif
}

void
CoulombPBCAB_CUDA::add(int groupID, RadFunctorType* ppot)
{
  RadFunctorType* savefunc = Vspec[groupID];
  CoulombPBCAB::add(groupID, ppot);
  RadFunctorType* rfunc = Vspec[groupID];
  if (rfunc != savefunc)
  {
    // Setup CUDA spline
    SRSplines[groupID] = new TextureSpline();
    SRSplines[groupID]->set(rfunc->data(), rfunc->size(),
                            rfunc->grid().rmin(), rfunc->grid().rmax());
  }
}

void
CoulombPBCAB_CUDA::setupLongRangeGPU()
{
  StructFact &SK = *(ElecRef.SK);
  Numk = SK.KLists.numk;
  gpu::host_vector<CUDA_PRECISION_FULL> kpointsHost(OHMMS_DIM*Numk);
  for (int ik=0; ik<Numk; ik++)
    for (int dim=0; dim<OHMMS_DIM; dim++)
      kpointsHost[ik*OHMMS_DIM+dim] = SK.KLists.kpts_cart[ik][dim];
  kpointsGPU = kpointsHost;
  gpu::host_vector<CUDA_PRECISION_FULL> FkHost(Numk);
  for (int ik=0; ik<Numk; ik++)
    FkHost[ik] = AB->Fk[ik];
  FkGPU = FkHost;
  // Now compute Rhok for the ions
  RhokIonsGPU.resize(NumIonSpecies);
  gpu::host_vector<CUDA_PRECISION_FULL> RhokIons_host(2*Numk);
  for (int sp=0; sp<NumIonSpecies; sp++)
  {
    for (int ik=0; ik < Numk; ik++)
    {
      PosType k = SK.KLists.kpts_cart[ik];
      RhokIons_host[2*ik+0] = 0.0;
      RhokIons_host[2*ik+1] = 0.0;
      for (int ion=IonFirst[sp]; ion<=IonLast[sp]; ion++)
      {
        PosType ipos = SortedIons[ion];
        RealType phase = dot(k,ipos);
        double s,c;
        sincos(phase, &s, &c);
        RhokIons_host[2*ik+0] += c;
        RhokIons_host[2*ik+1] += s;
      }
    }
    RhokIonsGPU[sp].set_name ("CoulombPBCABTemp::RhokIonsGPU");
    RhokIonsGPU[sp] = RhokIons_host;
  }
}


void
CoulombPBCAB_CUDA::addEnergy(MCWalkerConfiguration &W,
                             std::vector<RealType> &LocalEnergy)
{
  std::vector<Walker_t*> &walkers = W.WalkerList;
  // Short-circuit for constant contribution (e.g. fixed ions)
  // if (!is_active) {
  //   for (int iw=0; iw<walkers.size(); iw++) {
  // 	walkers[iw]->getPropertyBase()[NUMPROPERTIES+myIndex] = Value;
  // 	LocalEnergy[iw] += Value;
  //   }
  //   return;
  // }
  int nw = walkers.size();
  int N = NumElecs;
  if (SumGPU.size() < nw)
  {
    SumGPU.resize(nw);
    SumHost.resize(nw);
    RhokElecGPU.resize(2*nw*Numk);
    RhokIonsGPU.resize(NumIonSpecies);
    for (int sp=0; sp<NumIonSpecies; sp++)
      RhokIonsGPU[sp].resize(2*Numk);
    SumHost.resize(nw);
    RhoklistGPU.resize(nw);
    RhoklistHost.resize(nw);
  }
  for (int iw=0; iw<nw; iw++)
  {
    RhoklistHost[iw] = &(RhokElecGPU.data()[2*Numk*iw]);
    SumHost[iw] = 0.0;
  }
  RhoklistGPU = RhoklistHost;
  SumGPU = SumHost;
  // First, do short-range part
  std::vector<double> esum(nw, 0.0);
  for (int sp=0; sp<NumIonSpecies; sp++)
  {
    if (SRSplines[sp])
    {
      CoulombAB_SR_Sum
      (W.RList_GPU.data(), N, IGPU.data(), IonFirst[sp], IonLast[sp],
       SRSplines[sp]->rMax, SRSplines[sp]->NumPoints,
       SRSplines[sp]->MyTexture, L.data(), Linv.data(), SumGPU.data(), nw);
      SumHost = SumGPU;
      for (int iw=0; iw<nw; iw++)
        esum[iw] += Zspec[sp]*Qspec[0]* SumHost[iw];
    }
  }
  // Now, do long-range part:
  int first = 0;
  int last  = N-1;
  eval_rhok_cuda(W.RList_GPU.data(), first, last, kpointsGPU.data(), Numk,
                 RhoklistGPU.data(), nw);
  for (int sp=0; sp<NumIonSpecies; sp++)
  {
    for (int iw=0; iw<nw; iw++)
      SumHost[iw] = 0.0;
    SumGPU = SumHost;
    eval_vk_sum_cuda(RhoklistGPU.data(), RhokIonsGPU[sp].data(),
                     FkGPU.data(), Numk, SumGPU.data(), nw);
    SumHost = SumGPU;
    for (int iw=0; iw<nw; iw++)
      esum[iw] += Zspec[sp]*Qspec[0]* SumHost[iw];
  }
// #ifdef DEBUG_CUDA_RHOK
//     gpu::host_vector<CUDA_PRECISION_FULL> RhokHost;
//     RhokHost = RhokGPU;
//     for (int ik=0; ik<Numk; ik++) {
//       std::complex<double> rhok(0.0, 0.0);
//       PosType k = PtclRef.SK->KLists.kpts_cart[ik];
//       for (int ir=0; ir<N; ir++) {
//     	PosType r = walkers[0]->R[ir];
//     	double s, c;
//     	double phase = dot(k,r);
//     	sincos(phase, &s, &c);
//     	rhok += std::complex<double>(c,s);
//       }
//       fprintf (stderr, "GPU:   %d   %14.6f  %14.6f\n",
//     	       ik, RhokHost[2*ik+0], RhokHost[2*ik+1]);
//       fprintf (stderr, "CPU:   %d   %14.6f  %14.6f\n",
//     	       ik, rhok.real(), rhok.imag());
//     }
// #endif
//     for (int sp1=0; sp1<NumSpecies; sp1++)
//       for (int sp2=sp1; sp2<NumSpecies; sp2++)
//     	eval_vk_sum_cuda(RhoklistsGPU[sp1].data(), RhoklistsGPU[sp2].data(),
//     			 FkGPU.data(), Numk, SumGPU.data(), nw);
//    SumHost = SumGPU;
  for (int iw=0; iw<walkers.size(); iw++)
  {
    // fprintf (stderr, "Energy = %18.6f\n", SumHost[iw]);
    walkers[iw]->getPropertyBase()[NUMPROPERTIES+myIndex] =
      esum[iw] + myConst;
    LocalEnergy[iw] += esum[iw] + myConst;
  }
}


QMCHamiltonianBase*
CoulombPBCAB_CUDA::makeClone(ParticleSet& qp, TrialWaveFunction& psi)
{
  CoulombPBCAB_CUDA* myclone=new CoulombPBCAB_CUDA(PtclA, qp, true);
  if(myGrid) myclone->myGrid=new GridType(*myGrid);
  for(int ig=0; ig<Vspec.size(); ++ig)
  {
    if(Vspec[ig])
    {
      RadFunctorType* apot=Vspec[ig]->makeClone();
      myclone->Vspec[ig]=apot;
      for(int iat=0; iat<PtclA.getTotalNum(); ++iat)
      {
        if(PtclA.GroupID[iat]==ig) myclone->Vat[iat]=apot;
      }
    }
    myclone->V0Spline = V0Spline;
    myclone->SRSplines[ig] = SRSplines[ig];
  }
  return myclone;
}

}

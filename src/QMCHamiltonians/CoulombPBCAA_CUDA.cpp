//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    

#include "QMCHamiltonians/CoulombPBCAA_CUDA.h"
#include "Particle/MCWalkerConfiguration.h"

namespace qmcplusplus
{

CoulombPBCAA_CUDA::CoulombPBCAA_CUDA
(ParticleSet& ref, bool active, bool cloning) :
  CoulombPBCAA(ref,active,cloning),
  PtclRef(ref),
  SumGPU("CoulombPBCAATemp::SumGPU"),
  kpointsGPU("CoulombPBCAATemp::kpointsGPU"),
  kshellGPU("CoulombPBCAATemp::kshellGPU"),
  FkGPU("CoulombPBCAATemp::FkGPU"),
  RhokGPU("CoulombPBCAATemp::RhokGPU"),
  L("CoulombPBCAATemp::L"),
  Linv("CoulombPBCAATemp::Linv")
{
#ifdef QMC_CUDA
  gpu::host_vector<CUDA_PRECISION_FULL> LHost(9), LinvHost(9);
  for (int i=0; i<3; i++)
    for (int j=0; j<3; j++)
    {
      LHost[3*i+j]    = ref.Lattice.a(j)[i];
      LinvHost[3*i+j] = ref.Lattice.b(i)[j];
    }
  L    = LHost;
  Linv = LinvHost;
  initBreakup(ref,cloning);
#endif
}

void CoulombPBCAA_CUDA::initBreakup(ParticleSet& P, bool cloning)
{
#ifdef QMC_CUDA
  if (!cloning)
  {
    SRSpline = new TextureSpline;
    SRSpline->set(rVs->data(), rVs->size(), rVs->grid().rmin(),
                  rVs->grid().rmax());
  }
  setupLongRangeGPU(P);
#endif
}

QMCHamiltonianBase* CoulombPBCAA_CUDA::makeClone(ParticleSet& qp, TrialWaveFunction& psi)
{
  CoulombPBCAA_CUDA *myclone;
  if(is_active)
    myclone = new CoulombPBCAA_CUDA(qp,is_active, true);
  else
    myclone = new CoulombPBCAA_CUDA(*this);//nothing needs to be re-evaluated
  myclone->SRSpline = SRSpline;
  return myclone;
}

void
CoulombPBCAA_CUDA::setupLongRangeGPU(ParticleSet &P)
{
  if (is_active)
  {
    StructFact &SK = *(P.SK);
    Numk = SK.KLists.numk;
    gpu::host_vector<CUDA_PRECISION_FULL> kpointsHost(OHMMS_DIM*Numk);
    for (int ik=0; ik<Numk; ik++)
      for (int dim=0; dim<OHMMS_DIM; dim++)
        kpointsHost[ik*OHMMS_DIM+dim] = SK.KLists.kpts_cart[ik][dim];
    kpointsGPU = kpointsHost;
    gpu::host_vector<CUDA_PRECISION_FULL> FkHost(Numk);
    for (int ik=0; ik<Numk; ik++)
      FkHost[ik] = AA->Fk[ik];
    FkGPU = FkHost;
  }
}

void
CoulombPBCAA_CUDA::addEnergy(MCWalkerConfiguration &W,
                             std::vector<RealType> &LocalEnergy)
{
  std::vector<Walker_t*> &walkers = W.WalkerList;
  // Short-circuit for constant contribution (e.g. fixed ions)
  if (!is_active)
  {
    for (int iw=0; iw<walkers.size(); iw++)
    {
      walkers[iw]->getPropertyBase()[NUMPROPERTIES+myIndex] = Value;
      LocalEnergy[iw] += Value;
    }
    return;
  }
  int nw = walkers.size();
  int N = NumCenters;
  if (SumGPU.size() < nw)
  {
    SumGPU.resize(nw);
    RhokGPU.resize(2*nw*Numk*NumSpecies);
    RhoklistsGPU.resize(NumSpecies);
    RhoklistsHost.resize(NumSpecies);
    for (int sp=0; sp<NumSpecies; sp++)
    {
      RhoklistsGPU[sp].set_name("CoulombPBCAATemp::RhoklistsGPU");
      RhoklistsGPU[sp].resize(nw);
      RhoklistsHost[sp].resize(nw);
    }
  }
  for (int iw=0; iw<nw; iw++)
  {
    for (int sp=0; sp<NumSpecies; sp++)
      //RhoklistsHost[sp][iw] = &(RhokGPU.data()[2*nw*Numk*sp + 2*Numk*iw]);
      RhoklistsHost[sp][iw] = W.WalkerList[iw]->get_rhok_ptr(sp);
  }
  for (int sp=0; sp<NumSpecies; sp++)
    RhoklistsGPU[sp] = RhoklistsHost[sp];
  // First, do short-range part
  CoulombAA_SR_Sum(W.RList_GPU.data(), N, SRSpline->rMax, SRSpline->NumPoints,
                   SRSpline->MyTexture, L.data(), Linv.data(),
                   SumGPU.data(), nw);
  // Now, do long-range part:
  for (int sp=0; sp<NumSpecies; sp++)
  {
    int first = PtclRef.first(sp);
    int last  = PtclRef.last(sp)-1;
    eval_rhok_cuda(W.RList_GPU.data(), first, last,
                   kpointsGPU.data(), Numk,
                   W.RhokLists_GPU[sp].data(), walkers.size());
  }
#ifdef DEBUG_CUDA_RHOK
  gpu::host_vector<CUDA_PRECISION_FULL> RhokHost;
  RhokHost = RhokGPU;
  for (int ik=0; ik<Numk; ik++)
  {
    std::complex<double> rhok(0.0, 0.0);
    PosType k = PtclRef.SK->KLists.kpts_cart[ik];
    for (int ir=0; ir<N; ir++)
    {
      PosType r = walkers[0]->R[ir];
      double s, c;
      double phase = dot(k,r);
      sincos(phase, &s, &c);
      rhok += std::complex<double>(c,s);
    }
    fprintf (stderr, "GPU:   %d   %14.6f  %14.6f\n",
             ik, RhokHost[2*ik+0], RhokHost[2*ik+1]);
    fprintf (stderr, "CPU:   %d   %14.6f  %14.6f\n",
             ik, rhok.real(), rhok.imag());
  }
#endif
  for (int sp1=0; sp1<NumSpecies; sp1++)
    for (int sp2=sp1; sp2<NumSpecies; sp2++)
      eval_vk_sum_cuda(W.RhokLists_GPU[sp1].data(), W.RhokLists_GPU[sp2].data(),
                       FkGPU.data(), Numk, SumGPU.data(), nw);
  SumHost = SumGPU;
  for (int iw=0; iw<walkers.size(); iw++)
  {
    // fprintf (stderr, "Energy = %18.6f\n", SumHost[iw]);
    walkers[iw]->getPropertyBase()[NUMPROPERTIES+myIndex] =
      SumHost[iw] + myConst;
    LocalEnergy[iw] += SumHost[iw] + myConst;
  }
}


}

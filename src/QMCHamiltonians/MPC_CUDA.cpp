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
    
    
#include "QMCHamiltonians/MPC_CUDA.h"
#include "QMCHamiltonians/CudaCoulomb.h"
#include "Particle/MCWalkerConfiguration.h"
// #include "Lattice/ParticleBConds.h"
// #include "OhmmsPETE/OhmmsArray.h"
// #include "OhmmsData/AttributeSet.h"
// #include "Particle/DistanceTable.h"
// #include "Particle/DistanceTableData.h"


namespace qmcplusplus
{

MPC_CUDA::MPC_CUDA(ParticleSet& ptcl, double cutoff) :
  MPC(ptcl, cutoff), SumGPU("MPC::SumGPU"),
  L("MPC::L"), Linv("MPC::Linv")
{
  initBreakup();

  //create virtual particle
  myPtcl.resize(omp_get_max_threads());
  app_log() << "MPC_CUDA::evalLR uses " << myPtcl.size() << " threads." << std::endl;
  for(int i=0; i<myPtcl.size(); ++i)
    myPtcl[i]=new ParticleSet(ptcl);
}

void
MPC_CUDA::initBreakup()
{
  gpu::host_vector<CUDA_PRECISION> LHost(9), LinvHost(9);
  for (int i=0; i<3; i++)
    for (int j=0; j<3; j++)
    {
      LHost[3*i+j]    = PtclRef->Lattice.a(j)[i];
      LinvHost[3*i+j] = PtclRef->Lattice.b(i)[j];
    }
  L    = LHost;
  Linv = LinvHost;
//  app_log() << "    Starting to copy MPC spline to GPU memory.\n";
//  CudaSpline = create_UBspline_3d_s_cuda_conv (VlongSpline);
//  app_log() << "    Finished copying MPC spline to GPU memory.\n";
}

QMCHamiltonianBase*
MPC_CUDA::makeClone(ParticleSet& qp, TrialWaveFunction& psi)
{
  // return new MPC(qp, Ecut);
  MPC_CUDA* newMPC = new MPC_CUDA(*this);
  newMPC->resetTargetParticleSet(qp);
  return newMPC;
}

void
MPC_CUDA::addEnergy(MCWalkerConfiguration &W,
                    std::vector<RealType> &LocalEnergy)
{
  init_Acuda();
  std::vector<Walker_t*> &walkers = W.WalkerList;
  const int nw = walkers.size();
  const int N = NParticles;
  if (SumGPU.size() < nw)
  {
    SumGPU.resize(nw, 1.25);
    SumHost.resize(nw);
  }
  std::vector<double> esum(nw, 0.0);

  //for (int iw=0; iw<nw; iw++)
  //  SumHost[iw] = 0.0;
  //SumGPU = SumHost;
  // First, do short-range part
  //MPC_SR_Sum (W.RList_GPU.data(), N,
  //            L.data(), Linv.data(), SumGPU.data(), nw);
  //SumHost = SumGPU;
  //for (int iw=0; iw<nw; iw++)  esum[iw] += SumHost[iw];
  //// Now, do long-range part:
  //MPC_LR_Sum (W.RList_GPU.data(), N, CudaSpline,
  //            Linv.data(), SumGPU.data(), nw);
  //SumHost = SumGPU;
  //for (int iw=0; iw<nw; iw++)   esum[iw] += SumHost[iw];

  if(myPtcl.size()==1)
  {//serial
    for(int iw=0; iw<nw; ++iw)
    {
      ParticleSet& p(*myPtcl[0]);
      p.R=W[iw]->R;
      esum[iw]=MPC::evalLR(p);
    }
    SumGPU = SumHost;
    MPC_SR_Sum (W.RList_GPU.data(), N, L.data(), Linv.data(), SumGPU.data(), nw);
    SumHost = SumGPU;
  }
  else
  {
#pragma omp parallel
    {
      int ip=omp_get_thread_num();
      if(ip)
      {
        int np=omp_get_max_threads()-1;
        int nw_thread=(nw+np-1)/np;
        int first=nw_thread*(ip-1);
        int last=(ip<np)? nw_thread*ip:nw;
        if(last>nw) last=nw;
        ParticleSet& p(*myPtcl[ip]);
        for(int iw=first; iw<last; ++iw)
        {
          p.R=W[iw]->R;
          esum[iw]=MPC::evalLR(p);
        }
      }
      else
      {
        SumGPU = SumHost;
        MPC_SR_Sum (W.RList_GPU.data(), N, L.data(), Linv.data(), SumGPU.data(), nw);
        SumHost = SumGPU;
      }
    }
  }

  for (int iw=0; iw<nw; iw++)
  {
    double e=esum[iw]+SumHost[iw]+Vconst;
    walkers[iw]->getPropertyBase()[NUMPROPERTIES+myIndex] = e;
    LocalEnergy[iw] += e;
  }
}
}

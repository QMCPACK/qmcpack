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
}

void
MPC_CUDA::initBreakup()
{
#ifdef QMC_CUDA
  gpu::host_vector<CUDA_PRECISION> LHost(9), LinvHost(9);
  for (int i=0; i<3; i++)
    for (int j=0; j<3; j++)
    {
      LHost[3*i+j]    = PtclRef->Lattice.a(j)[i];
      LinvHost[3*i+j] = PtclRef->Lattice.b(i)[j];
    }
  L    = LHost;
  Linv = LinvHost;
  app_log() << "    Starting to copy MPC spline to GPU memory.\n";
  CudaSpline = create_UBspline_3d_s_cuda_conv (VlongSpline);
  app_log() << "    Finished copying MPC spline to GPU memory.\n";
#endif
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
                    vector<RealType> &LocalEnergy)
{
  init_Acuda();
  vector<Walker_t*> &walkers = W.WalkerList;
  int nw = walkers.size();
  int N = NParticles;
  if (SumGPU.size() < nw)
  {
    SumGPU.resize(nw, 1.25);
    SumHost.resize(nw);
  }
  for (int iw=0; iw<nw; iw++)
    SumHost[iw] = 0.0;
  SumGPU = SumHost;
  // First, do short-range part
  vector<double> esum(nw, 0.0);
  MPC_SR_Sum (W.RList_GPU.data(), N,
              L.data(), Linv.data(), SumGPU.data(), nw);
  SumHost = SumGPU;
  for (int iw=0; iw<nw; iw++)  esum[iw] += SumHost[iw];
  // Now, do long-range part:
  MPC_LR_Sum (W.RList_GPU.data(), N, CudaSpline,
              Linv.data(), SumGPU.data(), nw);
  SumHost = SumGPU;
  for (int iw=0; iw<nw; iw++)   esum[iw] += SumHost[iw];
  for (int iw=0; iw<nw; iw++)
  {
    walkers[iw]->getPropertyBase()[NUMPROPERTIES+myIndex] = esum[iw] + Vconst;
    LocalEnergy[iw] += esum[iw];
  }
}
}

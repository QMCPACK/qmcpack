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
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "OneBodyJastrowOrbitalBspline.h"
#include "CudaSpline.h"
#include "Lattice/ParticleBConds.h"
#include "QMCWaveFunctions/detail/CUDA_legacy/BsplineJastrowCuda.h"
#include "QMCWaveFunctions/detail/CUDA_legacy/BsplineJastrowCudaPBC.h"
#include "Configuration.h"

namespace qmcplusplus
{
template<class FT>
void OneBodyJastrowOrbitalBspline<FT>::recompute(MCWalkerConfiguration& W, bool firstTime)
{}

template<class FT>
void OneBodyJastrowOrbitalBspline<FT>::reserve(PointerPool<gpu::device_vector<CTS::RealType>>& pool)
{}

template<class FT>
void OneBodyJastrowOrbitalBspline<FT>::checkInVariables(opt_variables_type& active)
{
  J1OrbitalSoA<BsplineFunctor<WaveFunctionComponent::RealType>>::checkInVariables(active);
  for (int i = 0; i < NumCenterGroups; i++)
    if (JBase::J1UniqueFunctors[i] != nullptr)
      GPUSplines[i]->set(*JBase::J1UniqueFunctors[i]);
}

template<class FT>
void OneBodyJastrowOrbitalBspline<FT>::addFunc(int ig, std::unique_ptr<FT> j, int jg)
{
  auto newSpline = std::make_unique<CudaSpline<CTS::RealType>>(*j);
  GPUSplines[ig] = newSpline.get();
  UniqueSplines.push_back(std::move(newSpline));
  J1OrbitalSoA<BsplineFunctor<WaveFunctionComponent::RealType>>::addFunc(ig, std::move(j), jg);
}

template<class FT>
void OneBodyJastrowOrbitalBspline<FT>::addLog(MCWalkerConfiguration& W, std::vector<RealType>& logPsi)
{
  auto& walkers = W.WalkerList;
  if (SumHost.size() < 4 * walkers.size())
  {
    SumGPU.resize(4 * walkers.size());
    SumHost.resize(4 * walkers.size());
    UpdateListHost.resize(walkers.size());
    UpdateListGPU.resize(walkers.size());
  }
  int numGL = 4 * N * walkers.size();
  if (this->GradLaplGPU.size() < numGL)
  {
    this->GradLaplGPU.resize(numGL);
    this->GradLaplHost.resize(numGL);
  }
  std::vector<CTS::RealType> RHost(OHMMS_DIM * N * walkers.size());
  for (int iw = 0; iw < walkers.size(); iw++)
  {
    Walker_t& walker = *(walkers[iw]);
    SumHost[iw]      = 0.0;
  }
  SumGPU     = SumHost;
  int efirst = 0;
  int elast  = N - 1;
  for (int cgroup = 0; cgroup < NumCenterGroups; cgroup++)
  {
    int cfirst                        = CenterFirst[cgroup];
    int clast                         = CenterLast[cgroup];
    CudaSpline<CTS::RealType>& spline = *(GPUSplines[cgroup]);
    if (GPUSplines[cgroup])
    {
      if (UsePBC)
        one_body_sum_PBC(C.data(), W.RList_GPU.data(), cfirst, clast, efirst, elast, spline.coefs.data(),
                         spline.coefs.size(), spline.rMax, L.data(), Linv.data(), SumGPU.data(), walkers.size());
      else
        one_body_sum(C.data(), W.RList_GPU.data(), cfirst, clast, efirst, elast, spline.coefs.data(),
                     spline.coefs.size(), spline.rMax, SumGPU.data(), walkers.size());
    }
    // Copy data back to CPU memory
  }
  SumHost = SumGPU;
  for (int iw = 0; iw < walkers.size(); iw++)
    logPsi[iw] -= SumHost[iw];
#ifdef CUDA_DEBUG
  DTD_BConds<double, 3, SUPERCELL_BULK> bconds;
  double host_sum = 0.0;
  for (int cptcl = 0; cptcl < CenterRef.getTotalNum(); cptcl++)
  {
    PosType c = CenterRef.R[cptcl];
    for (int eptcl = 0; eptcl < N; eptcl++)
    {
      PosType r    = walkers[0]->R[eptcl];
      PosType disp = r - c;
      double dist  = std::sqrt(bconds.apply(ElecRef.Lattice, disp));
      host_sum -= J1Unique[cptcl]->evaluate(dist);
    }
  }
  fprintf(stderr, "host = %25.16f\n", host_sum);
  fprintf(stderr, "cuda = %25.16f\n", logPsi[0]);
#endif
}

template<class FT>
void OneBodyJastrowOrbitalBspline<FT>::update(MCWalkerConfiguration* W,
                                              std::vector<Walker_t*>& walkers,
                                              int iat,
                                              std::vector<bool>* acc,
                                              int k)
{
  // for (int iw=0; iw<walkers.size(); iw++)
  //   UpdateListHost[iw] = (CTS::RealType*)walkers[iw]->R_GPU.data();
  // UpdateListGPU = UpdateListHost;
  // one_body_update(UpdateListGPU.data(), N, iat, N);
}

// #define DEBUG_DELAYED

template<class FT>
void OneBodyJastrowOrbitalBspline<FT>::ratio(MCWalkerConfiguration& W,
                                             int iat,
                                             std::vector<ValueType>& psi_ratios,
                                             std::vector<GradType>& grad,
                                             std::vector<ValueType>& lapl)
{
  auto& walkers = W.WalkerList;
  int N         = W.Rnew_GPU.size();
  int nw        = walkers.size();
  bool zero     = true;
  if (SumGPU.size() < 4 * N)
    SumGPU.resize(4 * N);
  for (int group = 0; group < NumCenterGroups; group++)
  {
    int first = CenterFirst[group];
    int last  = CenterLast[group];
    if (GPUSplines[group])
    {
      CudaSpline<CTS::RealType>& spline = *(GPUSplines[group]);
      // 	one_body_ratio_PBC (C.data(), W.RList_GPU.data(), first, last,
      //  			    (CTS::RealType*)W.Rnew_GPU.data(), iat,
      // 			    spline.coefs.data(), spline.coefs.size(),
      // 			    spline.rMax, L.data(), Linv.data(),
      // 			    SumGPU.data(), N);
      if (UsePBC)
      {
        bool use_fast_image = W.getLattice().SimulationCellRadius >= spline.rMax;
        one_body_ratio_grad_PBC(C.data(), W.RList_GPU.data(), first, last, (CTS::RealType*)W.Rnew_GPU.data(), iat,
                                spline.coefs.data(), spline.coefs.size(), nw, spline.rMax, L.data(), Linv.data(), zero,
                                SumGPU.data(), N, use_fast_image);
      }
      else
        one_body_ratio_grad(C.data(), W.RList_GPU.data(), first, last, (CTS::RealType*)W.Rnew_GPU.data(), iat,
                            spline.coefs.data(), spline.coefs.size(), nw, spline.rMax, zero, SumGPU.data(), N);
      zero = false;
    }
  }
  // Copy data back to CPU memory
  SumHost = SumGPU;
  for (int iw = 0; iw < N; iw++)
  {
#ifdef DEBUG_DELAYED
    if (iw % nw == 0)
      fprintf(stderr, "k = %i:\n", iw / nw);
    fprintf(stderr, "walker %i ratio 1B Jastrow: %f|(%f,%f,%f) -> ", iw % nw, psi_ratios[iw], grad[iw][0], grad[iw][1],
            grad[iw][2]);
#endif
    psi_ratios[iw] *= std::exp(-SumHost[4 * iw + 0]);
    grad[iw][0] -= SumHost[4 * iw + 1];
    grad[iw][1] -= SumHost[4 * iw + 2];
    grad[iw][2] -= SumHost[4 * iw + 3];
#ifdef DEBUG_DELAYED
    fprintf(stderr, "%f|(%f,%f,%f)\n", psi_ratios[iw], grad[iw][0], grad[iw][1], grad[iw][2]);
#endif
  }
#ifdef CUDA_DEBUG
  DTD_BConds<double, 3, SUPERCELL_BULK> bconds;
  int iw           = 0;
  Walker_t& walker = *(walkers[iw]);
  double host_sum  = 0.0;
  for (int cptcl = 0; cptcl < CenterRef.getTotalNum(); cptcl++)
  {
    FT* func     = J1Unique[cptcl];
    PosType disp = new_pos[iw] - CenterRef.R[cptcl];
    double dist  = std::sqrt(bconds.apply(ElecRef.Lattice, disp));
    host_sum += func->evaluate(dist);
    disp = walkers[iw]->R[iat] - CenterRef.R[cptcl];
    dist = std::sqrt(bconds.apply(ElecRef.Lattice, disp));
    host_sum -= func->evaluate(dist);
  }
  fprintf(stderr, "Host sum = %18.12e\n", host_sum);
  fprintf(stderr, "CUDA sum = %18.12e\n", SumHost[0]);
#endif
}

template<class FT>
void OneBodyJastrowOrbitalBspline<FT>::calcRatio(MCWalkerConfiguration& W,
                                                 int iat,
                                                 std::vector<ValueType>& psi_ratios,
                                                 std::vector<GradType>& grad,
                                                 std::vector<ValueType>& lapl)
{
  int N         = W.Rnew_GPU.size();
  auto& walkers = W.WalkerList;
  int nw        = walkers.size();
  int kd        = W.getkDelay();
  int k         = W.getkcurr() - (kd > 1);
  if (k < 0)
    k += W.getkupdate();
  int offset = 0;
  if (W.getklinear())
    offset = k * nw;
  bool zero = true;
  if (SumGPU.size() < 4 * nw)
    SumGPU.resize(4 * nw);
  for (int group = 0; group < NumCenterGroups; group++)
  {
    int first = CenterFirst[group];
    int last  = CenterLast[group];
    if (GPUSplines[group])
    {
      CudaSpline<CTS::RealType>& spline = *(GPUSplines[group]);
      if (UsePBC)
      {
        bool use_fast_image = W.getLattice().SimulationCellRadius >= spline.rMax;
        one_body_ratio_grad_PBC(C.data(), W.RList_GPU.data(), first, last,
                                &(((CTS::RealType*)W.Rnew_GPU.data())[3 * offset]), iat, spline.coefs.data(),
                                spline.coefs.size(), walkers.size(), spline.rMax, L.data(), Linv.data(), zero,
                                SumGPU.data(), nw, use_fast_image);
      }
      else
        one_body_ratio_grad(C.data(), W.RList_GPU.data(), first, last,
                            &(((CTS::RealType*)W.Rnew_GPU.data())[3 * offset]), iat, spline.coefs.data(),
                            spline.coefs.size(), walkers.size(), spline.rMax, zero, SumGPU.data(), nw);
      zero = false;
    }
  }
  // Copy data back to CPU memory
  gpu::streamsSynchronize();
  SumHost.asyncCopy(SumGPU);
  cudaEventRecord(gpu::ratioSyncOneBodyEvent, gpu::memoryStream);
}

template<class FT>
void OneBodyJastrowOrbitalBspline<FT>::addRatio(MCWalkerConfiguration& W,
                                                int iat,
                                                int k,
                                                std::vector<ValueType>& psi_ratios,
                                                std::vector<GradType>& grad,
                                                std::vector<ValueType>& lapl)
{
  int N         = W.Rnew_GPU.size();
  auto& walkers = W.WalkerList;
  int nw        = walkers.size();
  cudaEventSynchronize(gpu::ratioSyncOneBodyEvent);
  for (int iw = 0; iw < nw; iw++)
  {
#ifdef DEBUG_DELAYED
    fprintf(stderr, "-> 1B Jastrow walker %i: ratio = %f ; grad = (%f,%f,%f) -> ", iw, psi_ratios[iw], grad[iw][0],
            grad[iw][1], grad[iw][2]);
#endif
    psi_ratios[iw] *= std::exp(-SumHost[4 * iw + 0]);
    grad[iw][0] -= SumHost[4 * iw + 1];
    grad[iw][1] -= SumHost[4 * iw + 2];
    grad[iw][2] -= SumHost[4 * iw + 3];
#ifdef DEBUG_DELAYED
    fprintf(stderr, "ratio = %f ; grad = (%f,%f,%f)\n", psi_ratios[iw], grad[iw][0], grad[iw][1], grad[iw][2]);
#endif
  }
#ifdef CUDA_DEBUG
  DTD_BConds<double, 3, SUPERCELL_BULK> bconds;
  int iw           = 0;
  Walker_t& walker = *(walkers[iw]);
  double host_sum  = 0.0;
  for (int cptcl = 0; cptcl < CenterRef.getTotalNum(); cptcl++)
  {
    FT* func     = J1Unique[cptcl];
    PosType disp = new_pos[iw] - CenterRef.R[cptcl];
    double dist  = std::sqrt(bconds.apply(ElecRef.Lattice, disp));
    host_sum += func->evaluate(dist);
    disp = walkers[iw]->R[iat] - CenterRef.R[cptcl];
    dist = std::sqrt(bconds.apply(ElecRef.Lattice, disp));
    host_sum -= func->evaluate(dist);
  }
  fprintf(stderr, "Host sum = %18.12e\n", host_sum);
  fprintf(stderr, "CUDA sum = %18.12e\n", SumHost[0]);
#endif
}

template<class FT>
void OneBodyJastrowOrbitalBspline<FT>::NLratios(MCWalkerConfiguration& W,
                                                std::vector<NLjob>& jobList,
                                                std::vector<PosType>& quadPoints,
                                                std::vector<ValueType>& psi_ratios)
{
  auto& walkers         = W.WalkerList;
  float sim_cell_radius = W.getLattice().SimulationCellRadius;
  int njobs             = jobList.size();
  if (NL_JobListHost.size() < njobs)
  {
    NL_JobListHost.resize(njobs);
    NL_SplineCoefsListHost.resize(njobs);
    NL_NumCoefsHost.resize(njobs);
    NL_rMaxHost.resize(njobs);
  }
  if (NL_RatiosHost.size() < quadPoints.size())
  {
    int nalloc = quadPoints.size();
    NL_QuadPointsHost.resize(OHMMS_DIM * nalloc);
    NL_QuadPointsGPU.resize(OHMMS_DIM * nalloc, 1.25);
    NL_RatiosHost.resize(nalloc);
    NL_RatiosGPU.resize(nalloc, 1.25);
  }
  int iratio = 0;
  for (int ijob = 0; ijob < njobs; ijob++)
  {
    NLjob& job                      = jobList[ijob];
    NLjobGPU<CTS::RealType>& jobGPU = NL_JobListHost[ijob];
    jobGPU.R                        = (CTS::RealType*)walkers[job.walker]->R_GPU.data();
    jobGPU.Elec                     = job.elec;
    jobGPU.QuadPoints               = &(NL_QuadPointsGPU.data()[OHMMS_DIM * iratio]);
    jobGPU.NumQuadPoints            = job.numQuadPoints;
    jobGPU.Ratios                   = &(NL_RatiosGPU.data()[iratio]);
    iratio += job.numQuadPoints;
  }
  NL_JobListGPU = NL_JobListHost;
  // Copy quad points to GPU
  for (int iq = 0; iq < quadPoints.size(); iq++)
  {
    NL_RatiosHost[iq] = 1.0;
    for (int dim = 0; dim < OHMMS_DIM; dim++)
      NL_QuadPointsHost[OHMMS_DIM * iq + dim] = quadPoints[iq][dim];
  }
  NL_RatiosGPU     = NL_RatiosHost;
  NL_QuadPointsGPU = NL_QuadPointsHost;
  // Now, loop over electron groups
  for (int group = 0; group < NumCenterGroups; group++)
  {
    int first = CenterFirst[group];
    int last  = CenterLast[group];
    if (GPUSplines[group])
    {
      CudaSpline<CTS::RealType>& spline = *(GPUSplines[group]);
      if (UsePBC)
        one_body_NLratios_PBC(NL_JobListGPU.data(), C.data(), first, last, spline.coefs.data(), spline.coefs.size(),
                              spline.rMax, L.data(), Linv.data(), sim_cell_radius, njobs);
      else
        one_body_NLratios(NL_JobListGPU.data(), C.data(), first, last, spline.coefs.data(), spline.coefs.size(),
                          spline.rMax, njobs);
    }
  }
  NL_RatiosHost = NL_RatiosGPU;
  for (int i = 0; i < psi_ratios.size(); i++)
    psi_ratios[i] *= NL_RatiosHost[i];
}

template<class FT>
void OneBodyJastrowOrbitalBspline<FT>::calcGradient(MCWalkerConfiguration& W,
                                                    int iat,
                                                    int k,
                                                    std::vector<GradType>& grad)
{
  CTS::RealType sim_cell_radius = W.getLattice().SimulationCellRadius;
  auto& walkers                 = W.WalkerList;
  if (this->OneGradHost.size() < OHMMS_DIM * walkers.size())
  {
    this->OneGradHost.resize(walkers.size() * OHMMS_DIM);
    this->OneGradGPU.resize(walkers.size() * OHMMS_DIM, 1.25);
  }
  bool zero = true;
  for (int group = 0; group < NumCenterGroups; group++)
  {
    int first = CenterFirst[group];
    int last  = CenterLast[group];
    if (GPUSplines[group])
    {
      CudaSpline<CTS::RealType>& spline = *(GPUSplines[group]);
      if (UsePBC)
        one_body_gradient_PBC(W.RList_GPU.data(), iat, C.data(), first, last, spline.coefs.data(), spline.coefs.size(),
                              spline.rMax, L.data(), Linv.data(), zero, this->OneGradGPU.data(), walkers.size());
      else
        one_body_gradient(W.RList_GPU.data(), iat, C.data(), first, last, spline.coefs.data(), spline.coefs.size(),
                          spline.rMax, zero, this->OneGradGPU.data(), walkers.size());
      zero = false;
    }
  }
  // Copy data back to CPU memory
  gpu::streamsSynchronize();
  this->OneGradHost.asyncCopy(this->OneGradGPU);
  cudaEventRecord(gpu::gradientSyncOneBodyEvent, gpu::memoryStream);
}

template<class FT>
void OneBodyJastrowOrbitalBspline<FT>::addGradient(MCWalkerConfiguration& W, int iat, std::vector<GradType>& grad)
{
  auto& walkers = W.WalkerList;
  cudaEventSynchronize(gpu::gradientSyncOneBodyEvent);
  for (int iw = 0; iw < walkers.size(); iw++)
  {
#ifdef DEBUG_DELAYED
    fprintf(stderr, "1B Jastrow grad walker %i: (", iw);
#endif
    for (int dim = 0; dim < OHMMS_DIM; dim++)
    {
#ifdef DEBUG_DELAYED
      if (dim > 0)
        fprintf(stderr, ", ");
      fprintf(stderr, "%f (before: %f)", OneGradHost[OHMMS_DIM * iw + dim], grad[iw][dim]);
#endif
      grad[iw][dim] -= this->OneGradHost[OHMMS_DIM * iw + dim];
    }
#ifdef DEBUG_DELAYED
    fprintf(stderr, ")\n");
#endif
  }
}

template<class FT>
void OneBodyJastrowOrbitalBspline<FT>::gradLapl(MCWalkerConfiguration& W, GradMatrix& grad, ValueMatrix& lapl)
{
  auto& walkers = W.WalkerList;
  int numGL     = 4 * N * walkers.size();
  if (this->GradLaplGPU.size() < numGL)
  {
    this->GradLaplGPU.resize(numGL, 1.25);
    this->GradLaplHost.resize(numGL);
  }
  for (int i = 0; i < walkers.size() * 4 * this->N; i++)
    this->GradLaplHost[i] = 0.0;
  this->GradLaplGPU = this->GradLaplHost;
#ifdef CUDA_DEBUG
  std::vector<CTS::RealType> CPU_GradLapl(4 * N);
  DTD_BConds<double, 3, SUPERCELL_BULK> bconds;
  int iw = 0;
  for (int eptcl = 0; eptcl < N; eptcl++)
  {
    PosType grad(0.0, 0.0, 0.0);
    double lapl(0.0);
    for (int cptcl = 0; cptcl < CenterRef.getTotalNum(); cptcl++)
    {
      FT* func     = J1Unique[cptcl];
      PosType disp = walkers[iw]->R[eptcl] - CenterRef.R[cptcl];
      double dist  = std::sqrt(bconds.apply(ElecRef.Lattice, disp));
      double u, du, d2u;
      u = func->evaluate(dist, du, d2u);
      du /= dist;
      grad += disp * du;
      lapl += d2u + 2.0 * du;
    }
    CPU_GradLapl[eptcl * 4 + 0] = grad[0];
    CPU_GradLapl[eptcl * 4 + 1] = grad[1];
    CPU_GradLapl[eptcl * 4 + 2] = grad[2];
    CPU_GradLapl[eptcl * 4 + 3] = lapl;
  }
#endif
  for (int cgroup = 0; cgroup < NumCenterGroups; cgroup++)
  {
    int cfirst = CenterFirst[cgroup];
    int clast  = CenterLast[cgroup];
    int efirst = 0;
    int elast  = N - 1;
    if (GPUSplines[cgroup])
    {
      CudaSpline<CTS::RealType>& spline = *(GPUSplines[cgroup]);
      // spline.check_coefs();
      if (UsePBC)
        one_body_grad_lapl_PBC(C.data(), W.RList_GPU.data(), cfirst, clast, efirst, elast, spline.coefs.data(),
                               spline.coefs.size(), spline.rMax, L.data(), Linv.data(), this->GradLaplGPU.data(),
                               4 * this->N, walkers.size());
      else
        one_body_grad_lapl(C.data(), W.RList_GPU.data(), cfirst, clast, efirst, elast, spline.coefs.data(),
                           spline.coefs.size(), spline.rMax, this->GradLaplGPU.data(), 4 * this->N, walkers.size());
    }
  }
  // Copy data back to CPU memory
  this->GradLaplHost = this->GradLaplGPU;
#ifdef CUDA_DEBUG
  fprintf(stderr, "GPU  grad = %12.5e %12.5e %12.5e   Lapl = %12.5e\n", GradLaplHost[0], GradLaplHost[1],
          GradLaplHost[2], GradLaplHost[3]);
  fprintf(stderr, "CPU  grad = %12.5e %12.5e %12.5e   Lapl = %12.5e\n", CPU_GradLapl[0], CPU_GradLapl[1],
          CPU_GradLapl[2], CPU_GradLapl[3]);
#endif
  for (int iw = 0; iw < walkers.size(); iw++)
  {
    for (int ptcl = 0; ptcl < N; ptcl++)
    {
      for (int i = 0; i < OHMMS_DIM; i++)
        grad(iw, ptcl)[i] += this->GradLaplHost[4 * this->N * iw + 4 * ptcl + i];
      if (std::isnan(this->GradLaplHost[4 * this->N * iw + +4 * ptcl + 3]))
      {
        fprintf(stderr, "NAN in OneBodyJastrowOrbitalBspline<FT> laplacian.\n");
        abort();
      }
      lapl(iw, ptcl) += this->GradLaplHost[4 * this->N * iw + +4 * ptcl + 3];
    }
  }
}

template<class FT>
void OneBodyJastrowOrbitalBspline<FT>::resetParameters(const opt_variables_type& active)
{
  J1OrbitalSoA<BsplineFunctor<WaveFunctionComponent::RealType>>::resetParameters(active);
  for (int i = 0; i < NumCenterGroups; i++)
    if (JBase::J1UniqueFunctors[i] != nullptr)
      GPUSplines[i]->set(*this->J1UniqueFunctors[i]);
}

template<class FT>
void OneBodyJastrowOrbitalBspline<FT>::evaluateDerivatives(MCWalkerConfiguration& W,
                                                           const opt_variables_type& optvars,
                                                           RealMatrix_t& d_logpsi,
                                                           RealMatrix_t& dlapl_over_psi)
{
  CTS::RealType sim_cell_radius = W.getLattice().SimulationCellRadius;
  auto& walkers                 = W.WalkerList;
  int nw                        = walkers.size();
  if (DerivListGPU.size() < nw)
  {
    MaxCoefs = 0;
    for (int i = 0; i < UniqueSplines.size(); i++)
      MaxCoefs = std::max(MaxCoefs, (int)UniqueSplines[i]->coefs.size());
    // Round up to nearest 16 to allow coallesced GPU reads
    MaxCoefs = ((MaxCoefs + 15) / 16) * 16;
    SplineDerivsHost.resize(2 * MaxCoefs * nw);
    SplineDerivsGPU.resize(2 * MaxCoefs * nw);
    DerivListHost.resize(nw);
    DerivListGPU.resize(nw, 1.25);
    for (int iw = 0; iw < nw; iw++)
      DerivListHost[iw] = SplineDerivsGPU.data() + 2 * iw * MaxCoefs;
    DerivListGPU = DerivListHost;
  }
  int efirst = 0;
  int elast  = N - 1;
  for (int cgroup = 0; cgroup < NumCenterGroups; cgroup++)
  {
    int cfirst                        = CenterFirst[cgroup];
    int clast                         = CenterLast[cgroup];
    CudaSpline<CTS::RealType>& spline = *(GPUSplines[cgroup]);
    //       std::cerr << "cgroup = " << cgroup << std::endl;
    //       std::cerr << "cfirst = " << cfirst << "  clast = " << clast << std::endl;
    //       std::cerr << "spline.coefs.size() = " << spline.coefs.size() << std::endl;
    //       std::cerr << "spline.rMax = " << spline.rMax << std::endl;
    if (UsePBC)
      one_body_derivs_PBC(C.data(), W.RList_GPU.data(), W.GradList_GPU.data(), cfirst, clast, efirst, elast,
                          spline.coefs.size(), spline.rMax, L.data(), Linv.data(), sim_cell_radius, DerivListGPU.data(),
                          nw);
    else
      one_body_derivs(C.data(), W.RList_GPU.data(), W.GradList_GPU.data(), cfirst, clast, efirst, elast,
                      spline.coefs.size(), spline.rMax, DerivListGPU.data(), nw);
    // Copy data back to CPU memory
    SplineDerivsHost = SplineDerivsGPU;
    //       std::cerr << "SplineDerivsHost = " << std::endl;
    //       for (int i=0; i<SplineDerivsHost.size(); i++)
    //         std::cerr << SplineDerivsHost[i] << std::endl;
    opt_variables_type splineVars = this->J1UniqueFunctors[cgroup]->myVars;
    for (int iv = 0; iv < splineVars.size(); iv++)
    {
      int varIndex  = splineVars.Index[iv];
      int coefIndex = iv + 1;
      for (int iw = 0; iw < nw; iw++)
      {
        d_logpsi(iw, varIndex) += SplineDerivsHost[2 * (MaxCoefs * iw + coefIndex) + 0];
        dlapl_over_psi(iw, varIndex) += SplineDerivsHost[2 * (MaxCoefs * iw + coefIndex) + 1];
        //	  std::cerr << "dlapl_over_psi(iw,"<<varIndex<<") = " << d_logpsi(iw,varIndex) << std::endl;
      }
    }
    int varIndex  = splineVars.Index[1];
    int coefIndex = 0;
    for (int iw = 0; iw < nw; iw++)
    {
      d_logpsi(iw, varIndex) += SplineDerivsHost[2 * (MaxCoefs * iw + coefIndex) + 0];
      dlapl_over_psi(iw, varIndex) += SplineDerivsHost[2 * (MaxCoefs * iw + coefIndex) + 1];
    }
  }
}

// explicit instantiations of templates
template class OneBodyJastrowOrbitalBspline<BsplineFunctor<WaveFunctionComponent::RealType>>;

} // namespace qmcplusplus

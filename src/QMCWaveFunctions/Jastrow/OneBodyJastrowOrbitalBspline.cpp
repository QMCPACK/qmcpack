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

namespace qmcplusplus
{

void
OneBodyJastrowOrbitalBspline::recompute(MCWalkerConfiguration &W,
                                        bool firstTime)
{
}

void
OneBodyJastrowOrbitalBspline::reserve
(PointerPool<gpu::device_vector<CudaRealType> > &pool)
{
}

void
OneBodyJastrowOrbitalBspline::checkInVariables(opt_variables_type& active)
{
  OneBodyJastrowOrbital<BsplineFunctor<OrbitalBase::RealType> >::checkInVariables(active);
  for (int i=0; i<NumCenterGroups; i++)
    GPUSplines[i]->set(*Funique[i]);
}

void
OneBodyJastrowOrbitalBspline::addFunc(int ig, FT* j, int jg=-1)
{
  OneBodyJastrowOrbital<BsplineFunctor<OrbitalBase::RealType> >::addFunc(ig, j);
  CudaSpline<CudaReal> *newSpline = new CudaSpline<CudaReal>(*j);
  UniqueSplines.push_back(newSpline);
  // if(i==0) { //first time, assign everything
  //   for(int ig=0; ig<NumCenterGroups; ++ig)
  // 	if(GPUSplines[ig]==0) GPUSplines[ig]=newSpline;
  // }
  // else
  GPUSplines[ig]=newSpline;
}


void
OneBodyJastrowOrbitalBspline::addLog (MCWalkerConfiguration &W,
                                      std::vector<RealType> &logPsi)
{
  std::vector<Walker_t*> &walkers = W.WalkerList;
  if (SumHost.size() < 4*walkers.size())
  {
    SumGPU.resize(4*walkers.size());
    SumHost.resize(4*walkers.size());
    UpdateListHost.resize(walkers.size());
    UpdateListGPU.resize(walkers.size());
  }
  int numGL = 4*N*walkers.size();
  if (GradLaplGPU.size()  < numGL)
  {
    GradLaplGPU.resize(numGL);
    GradLaplHost.resize(numGL);
  }
  CudaReal RHost[OHMMS_DIM*N*walkers.size()];
  for (int iw=0; iw<walkers.size(); iw++)
  {
    Walker_t &walker = *(walkers[iw]);
    SumHost[iw] = 0.0;
  }
  SumGPU = SumHost;
  int efirst = 0;
  int elast = N-1;
  for (int cgroup=0; cgroup<NumCenterGroups; cgroup++)
  {
    int cfirst = CenterFirst[cgroup];
    int clast  = CenterLast[cgroup];
    CudaSpline<CudaReal> &spline = *(GPUSplines[cgroup]);
    if (GPUSplines[cgroup])
    {
      if (UsePBC)
        one_body_sum_PBC (C.data(), W.RList_GPU.data(),
                          cfirst, clast, efirst, elast,
                          spline.coefs.data(), spline.coefs.size(),
                          spline.rMax, L.data(), Linv.data(),
                          SumGPU.data(), walkers.size());
      else
        one_body_sum (C.data(), W.RList_GPU.data(),
                      cfirst, clast, efirst, elast,
                      spline.coefs.data(), spline.coefs.size(),
                      spline.rMax, SumGPU.data(), walkers.size());
    }
    // Copy data back to CPU memory
  }
  SumHost = SumGPU;
  for (int iw=0; iw<walkers.size(); iw++)
    logPsi[iw] -= SumHost[iw];
#ifdef CUDA_DEBUG
  DTD_BConds<double,3,SUPERCELL_BULK> bconds;
  double host_sum = 0.0;
  for (int cptcl=0; cptcl < CenterRef.getTotalNum(); cptcl++)
  {
    PosType c = CenterRef.R[cptcl];
    for (int eptcl=0; eptcl<N; eptcl++)
    {
      PosType r = walkers[0]->R[eptcl];
      PosType disp = r - c;
      double dist = std::sqrt(bconds.apply(ElecRef.Lattice, disp));
      host_sum -= Fs[cptcl]->evaluate(dist);
    }
  }
  fprintf (stderr, "host = %25.16f\n", host_sum);
  fprintf (stderr, "cuda = %25.16f\n", logPsi[0]);
#endif
}

void
OneBodyJastrowOrbitalBspline::update (std::vector<Walker_t*> &walkers, int iat)
{
  // for (int iw=0; iw<walkers.size(); iw++)
  //   UpdateListHost[iw] = (CudaReal*)walkers[iw]->R_GPU.data();
  // UpdateListGPU = UpdateListHost;
  // one_body_update(UpdateListGPU.data(), N, iat, walkers.size());
}

void
OneBodyJastrowOrbitalBspline::ratio
(MCWalkerConfiguration &W, int iat,
 std::vector<ValueType> &psi_ratios, std::vector<GradType>  &grad,
 std::vector<ValueType> &lapl)
{
  std::vector<Walker_t*> &walkers = W.WalkerList;
  bool zero = true;
  if (SumGPU.size() < 4*walkers.size())
    SumGPU.resize(4*walkers.size());
  for (int group=0; group<NumCenterGroups; group++)
  {
    int first = CenterFirst[group];
    int last  = CenterLast[group];
    if (GPUSplines[group])
    {
      CudaSpline<CudaReal> &spline = *(GPUSplines[group]);
// 	one_body_ratio_PBC (C.data(), W.RList_GPU.data(), first, last,
//  			    (CudaReal*)W.Rnew_GPU.data(), iat,
// 			    spline.coefs.data(), spline.coefs.size(),
// 			    spline.rMax, L.data(), Linv.data(),
// 			    SumGPU.data(), walkers.size());
      if (UsePBC)
      {
        bool use_fast_image = W.Lattice.SimulationCellRadius >= spline.rMax;
        one_body_ratio_grad_PBC (C.data(), W.RList_GPU.data(), first, last,
                                 (CudaReal*)W.Rnew_GPU.data(), iat,
                                 spline.coefs.data(), spline.coefs.size(),
                                 spline.rMax, L.data(), Linv.data(), zero,
                                 SumGPU.data(), walkers.size(), use_fast_image);
      }
      else
        one_body_ratio_grad (C.data(), W.RList_GPU.data(), first, last,
                             (CudaReal*)W.Rnew_GPU.data(), iat,
                             spline.coefs.data(), spline.coefs.size(),
                             spline.rMax, zero, SumGPU.data(), walkers.size());
      zero = false;
    }
  }
  // Copy data back to CPU memory
  SumHost = SumGPU;
  for (int iw=0; iw<walkers.size(); iw++)
  {
    psi_ratios[iw] *= std::exp(-SumHost[4*iw+0]);
    grad[iw][0] -= SumHost[4*iw+1];
    grad[iw][1] -= SumHost[4*iw+2];
    grad[iw][2] -= SumHost[4*iw+3];
  }
#ifdef CUDA_DEBUG
  DTD_BConds<double,3,SUPERCELL_BULK> bconds;
  int iw = 0;
  Walker_t &walker = *(walkers[iw]);
  double host_sum = 0.0;
  for (int cptcl=0; cptcl<CenterRef.getTotalNum(); cptcl++)
  {
    FT* func = Fs[cptcl];
    PosType disp = new_pos[iw] - CenterRef.R[cptcl];
    double dist = std::sqrt(bconds.apply(ElecRef.Lattice, disp));
    host_sum += func->evaluate(dist);
    disp = walkers[iw]->R[iat] - CenterRef.R[cptcl];
    dist = std::sqrt(bconds.apply(ElecRef.Lattice, disp));
    host_sum -= func->evaluate(dist);
  }
  fprintf (stderr, "Host sum = %18.12e\n", host_sum);
  fprintf (stderr, "CUDA sum = %18.12e\n", SumHost[0]);
#endif
}

void
OneBodyJastrowOrbitalBspline::calcRatio
(MCWalkerConfiguration &W, int iat,
 std::vector<ValueType> &psi_ratios, std::vector<GradType>  &grad,
 std::vector<ValueType> &lapl)
{
  std::vector<Walker_t*> &walkers = W.WalkerList;
  bool zero = true;
  if (SumGPU.size() < 4*walkers.size())
    SumGPU.resize(4*walkers.size());
  for (int group=0; group<NumCenterGroups; group++)
  {
    int first = CenterFirst[group];
    int last  = CenterLast[group];
    if (GPUSplines[group])
    {
      CudaSpline<CudaReal> &spline = *(GPUSplines[group]);
      if (UsePBC)
      {
        bool use_fast_image = W.Lattice.SimulationCellRadius >= spline.rMax;
        one_body_ratio_grad_PBC (C.data(), W.RList_GPU.data(), first, last,
                                 (CudaReal*)W.Rnew_GPU.data(), iat,
                                 spline.coefs.data(), spline.coefs.size(),
                                 spline.rMax, L.data(), Linv.data(), zero,
                                 SumGPU.data(), walkers.size(), use_fast_image);
      }
      else
        one_body_ratio_grad (C.data(), W.RList_GPU.data(), first, last,
                             (CudaReal*)W.Rnew_GPU.data(), iat,
                             spline.coefs.data(), spline.coefs.size(),
                             spline.rMax, zero, SumGPU.data(), walkers.size());
      zero = false;
    }
  }
  // Copy data back to CPU memory
  gpu::streamsSynchronize();
  SumHost.asyncCopy(SumGPU);
  cudaEventRecord(gpu::ratioSyncOneBodyEvent, gpu::memoryStream);
}

void
OneBodyJastrowOrbitalBspline::addRatio
(MCWalkerConfiguration &W, int iat,
 std::vector<ValueType> &psi_ratios, std::vector<GradType>  &grad,
 std::vector<ValueType> &lapl)
{
  std::vector<Walker_t*> &walkers = W.WalkerList;
  cudaEventSynchronize(gpu::ratioSyncOneBodyEvent);
  for (int iw=0; iw<walkers.size(); iw++)
  {
    psi_ratios[iw] *= std::exp(-SumHost[4*iw+0]);
    grad[iw][0] -= SumHost[4*iw+1];
    grad[iw][1] -= SumHost[4*iw+2];
    grad[iw][2] -= SumHost[4*iw+3];
  }
#ifdef CUDA_DEBUG
  DTD_BConds<double,3,SUPERCELL_BULK> bconds;
  int iw = 0;
  Walker_t &walker = *(walkers[iw]);
  double host_sum = 0.0;
  for (int cptcl=0; cptcl<CenterRef.getTotalNum(); cptcl++)
  {
    FT* func = Fs[cptcl];
    PosType disp = new_pos[iw] - CenterRef.R[cptcl];
    double dist = std::sqrt(bconds.apply(ElecRef.Lattice, disp));
    host_sum += func->evaluate(dist);
    disp = walkers[iw]->R[iat] - CenterRef.R[cptcl];
    dist = std::sqrt(bconds.apply(ElecRef.Lattice, disp));
    host_sum -= func->evaluate(dist);
  }
  fprintf (stderr, "Host sum = %18.12e\n", host_sum);
  fprintf (stderr, "CUDA sum = %18.12e\n", SumHost[0]);
#endif
}


void
OneBodyJastrowOrbitalBspline::NLratios
(MCWalkerConfiguration &W,  std::vector<NLjob> &jobList,
 std::vector<PosType> &quadPoints, std::vector<ValueType> &psi_ratios)
{
  std::vector<Walker_t*> &walkers = W.WalkerList;
  float sim_cell_radius = W.Lattice.SimulationCellRadius;
  int njobs = jobList.size();
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
    NL_QuadPointsHost.resize(OHMMS_DIM*nalloc);
    NL_QuadPointsGPU.resize(OHMMS_DIM*nalloc, 1.25);
    NL_RatiosHost.resize(nalloc);
    NL_RatiosGPU.resize(nalloc, 1.25);
  }
  int iratio = 0;
  for (int ijob=0; ijob < njobs; ijob++)
  {
    NLjob &job = jobList[ijob];
    NLjobGPU<CudaReal> &jobGPU = NL_JobListHost[ijob];
    jobGPU.R             = (CudaReal*)walkers[job.walker]->R_GPU.data();
    jobGPU.Elec          = job.elec;
    jobGPU.QuadPoints    = &(NL_QuadPointsGPU.data()[OHMMS_DIM*iratio]);
    jobGPU.NumQuadPoints = job.numQuadPoints;
    jobGPU.Ratios        = &(NL_RatiosGPU.data()[iratio]);
    iratio += job.numQuadPoints;
  }
  NL_JobListGPU         = NL_JobListHost;
  // Copy quad points to GPU
  for (int iq=0; iq<quadPoints.size(); iq++)
  {
    NL_RatiosHost[iq] = 1.0;
    for (int dim=0; dim<OHMMS_DIM; dim++)
      NL_QuadPointsHost[OHMMS_DIM*iq + dim] = quadPoints[iq][dim];
  }
  NL_RatiosGPU = NL_RatiosHost;
  NL_QuadPointsGPU = NL_QuadPointsHost;
  // Now, loop over electron groups
  for (int group=0; group<NumCenterGroups; group++)
  {
    int first = CenterFirst[group];
    int last  = CenterLast[group];
    if (GPUSplines[group])
    {
      CudaSpline<CudaReal> &spline = *(GPUSplines[group]);
      if (UsePBC)
        one_body_NLratios_PBC(NL_JobListGPU.data(), C.data(), first, last,
                              spline.coefs.data(), spline.coefs.size(),
                              spline.rMax, L.data(), Linv.data(), sim_cell_radius,
                              njobs);
      else
        one_body_NLratios(NL_JobListGPU.data(), C.data(), first, last,
                          spline.coefs.data(), spline.coefs.size(),
                          spline.rMax, njobs);
    }
  }
  NL_RatiosHost = NL_RatiosGPU;
  for (int i=0; i < psi_ratios.size(); i++)
    psi_ratios[i] *= NL_RatiosHost[i];
}


void OneBodyJastrowOrbitalBspline::calcGradient
(MCWalkerConfiguration &W, int iat, std::vector<GradType> &grad)
{
  CudaReal sim_cell_radius = W.Lattice.SimulationCellRadius;
  std::vector<Walker_t*> &walkers = W.WalkerList;
  if (OneGradHost.size() < OHMMS_DIM*walkers.size())
  {
    OneGradHost.resize (walkers.size()*OHMMS_DIM);
    OneGradGPU.resize (walkers.size()*OHMMS_DIM, 1.25);
  }
  bool zero = true;
  for (int group=0; group<NumCenterGroups; group++)
  {
    int first = CenterFirst[group];
    int last  = CenterLast[group];
    if (GPUSplines[group])
    {
      CudaSpline<CudaReal> &spline = *(GPUSplines[group]);
      if (UsePBC)
        one_body_gradient_PBC (W.RList_GPU.data(), iat, C.data(), first, last,
                               spline.coefs.data(), spline.coefs.size(),
                               spline.rMax, L.data(), Linv.data(),
                               zero, OneGradGPU.data(), walkers.size());
      else
        one_body_gradient (W.RList_GPU.data(), iat, C.data(), first, last,
                           spline.coefs.data(), spline.coefs.size(),
                           spline.rMax, zero, OneGradGPU.data(), walkers.size());
      zero = false;
    }
  }
  // Copy data back to CPU memory
  gpu::streamsSynchronize();
  OneGradHost.asyncCopy(OneGradGPU);
  cudaEventRecord(gpu::gradientSyncOneBodyEvent, gpu::memoryStream);
}

void OneBodyJastrowOrbitalBspline::addGradient
(MCWalkerConfiguration &W, int iat, std::vector<GradType> &grad)
{
  std::vector<Walker_t*> &walkers = W.WalkerList;
  cudaEventSynchronize(gpu::gradientSyncOneBodyEvent);
  for (int iw=0; iw<walkers.size(); iw++)
    for (int dim=0; dim<OHMMS_DIM; dim++)
      grad[iw][dim] -= OneGradHost[OHMMS_DIM*iw+dim];
}


void
OneBodyJastrowOrbitalBspline::gradLapl (MCWalkerConfiguration &W,
                                        GradMatrix_t &grad,
                                        ValueMatrix_t &lapl)
{
  std::vector<Walker_t*> &walkers = W.WalkerList;
  int numGL = 4*N*walkers.size();
  if (GradLaplGPU.size()  < numGL)
  {
    GradLaplGPU.resize(numGL, 1.25);
    GradLaplHost.resize(numGL);
  }
  for (int i=0; i<walkers.size()*4*N; i++)
    GradLaplHost[i] = 0.0;
  GradLaplGPU = GradLaplHost;
#ifdef CUDA_DEBUG
  std::vector<CudaReal> CPU_GradLapl(4*N);
  DTD_BConds<double,3,SUPERCELL_BULK> bconds;
  int iw = 0;
  for (int eptcl=0; eptcl<N; eptcl++)
  {
    PosType grad(0.0, 0.0, 0.0);
    double lapl(0.0);
    for (int cptcl=0; cptcl<CenterRef.getTotalNum(); cptcl++)
    {
      FT* func = Fs[cptcl];
      PosType disp = walkers[iw]->R[eptcl] - CenterRef.R[cptcl];
      double dist = std::sqrt(bconds.apply(ElecRef.Lattice, disp));
      double u, du, d2u;
      u = func->evaluate(dist, du, d2u);
      du /= dist;
      grad += disp * du;
      lapl += d2u + 2.0*du;
    }
    CPU_GradLapl[eptcl*4+0] = grad[0];
    CPU_GradLapl[eptcl*4+1] = grad[1];
    CPU_GradLapl[eptcl*4+2] = grad[2];
    CPU_GradLapl[eptcl*4+3] = lapl;
  }
#endif
  for (int cgroup=0; cgroup<NumCenterGroups; cgroup++)
  {
    int cfirst = CenterFirst[cgroup];
    int clast  = CenterLast[cgroup];
    int efirst = 0;
    int elast  = N-1;
    if (GPUSplines[cgroup])
    {
      CudaSpline<CudaReal> &spline = *(GPUSplines[cgroup]);
      // spline.check_coefs();
      if (UsePBC)
        one_body_grad_lapl_PBC (C.data(), W.RList_GPU.data(),
                                cfirst, clast, efirst, elast,
                                spline.coefs.data(), spline.coefs.size(),
                                spline.rMax, L.data(), Linv.data(),
                                GradLaplGPU.data(), 4*N, walkers.size());
      else
        one_body_grad_lapl (C.data(), W.RList_GPU.data(),
                            cfirst, clast, efirst, elast,
                            spline.coefs.data(), spline.coefs.size(),
                            spline.rMax, GradLaplGPU.data(),
                            4*N, walkers.size());
    }
  }
  // Copy data back to CPU memory
  GradLaplHost = GradLaplGPU;
#ifdef CUDA_DEBUG
  fprintf (stderr, "GPU  grad = %12.5e %12.5e %12.5e   Lapl = %12.5e\n",
           GradLaplHost[0],  GradLaplHost[1], GradLaplHost[2], GradLaplHost[3]);
  fprintf (stderr, "CPU  grad = %12.5e %12.5e %12.5e   Lapl = %12.5e\n",
           CPU_GradLapl[0],  CPU_GradLapl[1], CPU_GradLapl[2], CPU_GradLapl[3]);
#endif
  for (int iw=0; iw<walkers.size(); iw++)
  {
    for (int ptcl=0; ptcl<N; ptcl++)
    {
      for (int i=0; i<OHMMS_DIM; i++)
        grad(iw,ptcl)[i] += GradLaplHost[4*N*iw + 4*ptcl + i];
      if (std::isnan(GradLaplHost[4*N*iw+ + 4*ptcl +3]))
      {
        fprintf (stderr, "NAN in OneBodyJastrowOrbitalBspline laplacian.\n");
        abort();
      }
      lapl(iw,ptcl) += GradLaplHost[4*N*iw+ + 4*ptcl +3];
    }
  }
}

void
OneBodyJastrowOrbitalBspline::resetParameters(const opt_variables_type& active)
{
  OneBodyJastrowOrbital<BsplineFunctor<OrbitalBase::RealType> >::resetParameters(active);
  for (int i=0; i<NumCenterGroups; i++)
    GPUSplines[i]->set(*Funique[i]);
}

void
OneBodyJastrowOrbitalBspline::evaluateDerivatives
(MCWalkerConfiguration &W, const opt_variables_type& optvars,
 RealMatrix_t &d_logpsi, RealMatrix_t &dlapl_over_psi)
{
  CudaReal sim_cell_radius = W.Lattice.SimulationCellRadius;
  std::vector<Walker_t*> &walkers = W.WalkerList;
  int nw = walkers.size();
  if (DerivListGPU.size() < nw)
  {
    MaxCoefs = 0;
    for (int i=0; i<UniqueSplines.size(); i++)
      MaxCoefs = std::max(MaxCoefs, (int)UniqueSplines[i]->coefs.size());
    // Round up to nearest 16 to allow coallesced GPU reads
    MaxCoefs = ((MaxCoefs+15)/16)*16;
    SplineDerivsHost.resize(2*MaxCoefs*nw);
    SplineDerivsGPU.resize(2*MaxCoefs*nw);
    DerivListHost.resize(nw);
    DerivListGPU.resize(nw, 1.25);
    for (int iw=0; iw<nw; iw++)
      DerivListHost[iw] = SplineDerivsGPU.data() + 2*iw*MaxCoefs;
    DerivListGPU = DerivListHost;
  }
  int efirst = 0;
  int elast = N-1;
  for (int cgroup=0; cgroup < NumCenterGroups; cgroup++)
  {
    int cfirst = CenterFirst[cgroup];
    int clast  = CenterLast[cgroup];
    CudaSpline<CudaReal> &spline = *(GPUSplines[cgroup]);
//       std::cerr << "cgroup = " << cgroup << std::endl;
//       std::cerr << "cfirst = " << cfirst << "  clast = " << clast << std::endl;
//       std::cerr << "spline.coefs.size() = " << spline.coefs.size() << std::endl;
//       std::cerr << "spline.rMax = " << spline.rMax << std::endl;
    if (UsePBC)
      one_body_derivs_PBC (C.data(), W.RList_GPU.data(), W.GradList_GPU.data(),
                           cfirst, clast, efirst, elast,
                           spline.coefs.size(), spline.rMax, L.data(),
                           Linv.data(), sim_cell_radius, DerivListGPU.data(),nw);
    else
      one_body_derivs (C.data(), W.RList_GPU.data(), W.GradList_GPU.data(),
                       cfirst, clast, efirst, elast,
                       spline.coefs.size(), spline.rMax, DerivListGPU.data(),nw);
    // Copy data back to CPU memory
    SplineDerivsHost = SplineDerivsGPU;
//       std::cerr << "SplineDerivsHost = " << std::endl;
//       for (int i=0; i<SplineDerivsHost.size(); i++)
// 	cerr << SplineDerivsHost[i] << std::endl;
    opt_variables_type splineVars = Funique[cgroup]->myVars;
    for (int iv=0; iv<splineVars.size(); iv++)
    {
      int varIndex = splineVars.Index[iv];
      int coefIndex = iv+1;
      for (int iw=0; iw<nw; iw++)
      {
        d_logpsi(iw,varIndex) +=
          SplineDerivsHost[2*(MaxCoefs*iw+coefIndex)+0];
        dlapl_over_psi(iw,varIndex) +=
          SplineDerivsHost[2*(MaxCoefs*iw+coefIndex)+1];
        //	  std::cerr << "dlapl_over_psi(iw,"<<varIndex<<") = " << d_logpsi(iw,varIndex) << std::endl;
      }
    }
    int varIndex = splineVars.Index[1];
    int coefIndex = 0;
    for (int iw=0; iw<nw; iw++)
    {
      d_logpsi(iw,varIndex) +=
        SplineDerivsHost[2*(MaxCoefs*iw+coefIndex)+0];
      dlapl_over_psi(iw,varIndex) +=
        SplineDerivsHost[2*(MaxCoefs*iw+coefIndex)+1];
    }
  }
}
}

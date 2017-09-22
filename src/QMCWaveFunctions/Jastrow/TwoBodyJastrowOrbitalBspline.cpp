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
    
    
#include "TwoBodyJastrowOrbitalBspline.h"
#include "CudaSpline.h"
#include "Lattice/ParticleBConds.h"
#include "QMCWaveFunctions/Jastrow/BsplineJastrowCuda.h"
#include "QMCWaveFunctions/Jastrow/BsplineJastrowCudaPBC.h"
#include <unistd.h>


namespace qmcplusplus
{
void
TwoBodyJastrowOrbitalBspline::freeGPUmem()
{
  UpdateListGPU.clear();
  SumGPU.clear();
  GradLaplGPU.clear();
  OneGradGPU.clear();
  SplineDerivsGPU.clear();
  DerivListGPU.clear();
  NL_SplineCoefsListGPU.clear();
  NL_JobListGPU.clear();
  NL_NumCoefsGPU.clear();
  NL_NumQuadPointsGPU.clear();
  NL_rMaxGPU.clear();
  NL_QuadPointsGPU.clear();
  NL_RatiosGPU.clear();
}


void
TwoBodyJastrowOrbitalBspline::recompute(MCWalkerConfiguration &W,
                                        bool firstTime)
{
}

void
TwoBodyJastrowOrbitalBspline::reserve
(PointerPool<gpu::device_vector<CudaRealType> > &pool)
{
}

void
TwoBodyJastrowOrbitalBspline::checkInVariables(opt_variables_type& active)
{
  TwoBodyJastrowOrbital<BsplineFunctor<OrbitalBase::RealType> >::checkInVariables(active);
  for (int i=0; i<NumGroups*NumGroups; i++)
    GPUSplines[i]->set(*F[i]);
}

void
TwoBodyJastrowOrbitalBspline::addFunc(int ia, int ib, FT* j)
{
  TwoBodyJastrowOrbital<BsplineFunctor<OrbitalBase::RealType> >::addFunc(ia, ib, j);
  CudaSpline<CudaReal> *newSpline = new CudaSpline<CudaReal>(*j);
  UniqueSplines.push_back(newSpline);
  if(ia==ib)
  {
    if(ia==0)   //first time, assign everything
    {
      int ij=0;
      for(int ig=0; ig<NumGroups; ++ig)
        for(int jg=0; jg<NumGroups; ++jg, ++ij)
          if(GPUSplines[ij]==0) GPUSplines[ij]=newSpline;
    }
    else
      GPUSplines[ia*NumGroups+ib]=newSpline;
  }
  else
  {
    if(PtclRef.R.size()==2)
    {
      // a very special case, 1 up + 1 down
      // uu/dd was prevented by the builder
      for(int ig=0; ig<NumGroups; ++ig)
        for(int jg=0; jg<NumGroups; ++jg)
          GPUSplines[ig*NumGroups+jg]=newSpline;
    }
    else
    {
      // generic case
      GPUSplines[ia*NumGroups+ib]=newSpline;
      GPUSplines[ib*NumGroups+ia]=newSpline;
    }
  }
}


void
TwoBodyJastrowOrbitalBspline::addLog (MCWalkerConfiguration &W,
                                      std::vector<RealType> &logPsi)
{
  std::vector<Walker_t*> &walkers = W.WalkerList;
  if (SumGPU.size() < 4*walkers.size())
  {
    SumGPU.resize(4*walkers.size());
    SumHost.resize(4*walkers.size());
    OneGradGPU.resize (walkers.size()*OHMMS_DIM);
    OneGradHost.resize(walkers.size()*OHMMS_DIM);
    UpdateListHost.resize(walkers.size());
    UpdateListGPU.resize(walkers.size());
  }
  int numGL = 4*N*walkers.size();
  if (GradLaplGPU.size()  < numGL)
  {
    GradLaplGPU.resize(numGL);
    GradLaplHost.resize(numGL);
  }
  for (int iw=0; iw<walkers.size(); iw++)
  {
    Walker_t &walker = *(walkers[iw]);
    SumHost[iw] = 0.0;
    CudaReal *dest = (CudaReal*)walker.R_GPU.data();
  }
  SumGPU = SumHost;
//     DTD_BConds<double,3,SUPERCELL_BULK> bconds;
//     double host_sum = 0.0;
  for (int group1=0; group1<PtclRef.groups(); group1++)
  {
    int first1 = PtclRef.first(group1);
    int last1  = PtclRef.last(group1) -1;
    for (int group2=group1; group2<PtclRef.groups(); group2++)
    {
      int first2 = PtclRef.first(group2);
      int last2  = PtclRef.last(group2) -1;
// 	double factor = (group1 == group2) ? 0.5 : 1.0;
// 	for (int e1=first1; e1 <= last1; e1++)
// 	  for (int e2=first2; e2 <= last2; e2++) {
// 	    PosType disp = walkers[0]->R[e2] - walkers[0]->R[e1];
// 	    double dist = std::sqrt(bconds.apply(PtclRef.Lattice, disp));
// 	    if (e1 != e2)
// 	      host_sum -= factor * F[group2*NumGroups + group1]->evaluate(dist);
// 	  }
      CudaSpline<CudaReal> &spline = *(GPUSplines[group1*NumGroups+group2]);
      if (UsePBC)
        two_body_sum_PBC (W.RList_GPU.data(), first1, last1, first2, last2,
                          spline.coefs.data(), spline.coefs.size(),
                          spline.rMax, L.data(), Linv.data(),
                          SumGPU.data(), walkers.size());
      else
        two_body_sum (W.RList_GPU.data(), first1, last1, first2, last2,
                      spline.coefs.data(), spline.coefs.size(),
                      spline.rMax, SumGPU.data(), walkers.size());
    }
  }
  // Copy data back to CPU memory
  SumHost = SumGPU;
  for (int iw=0; iw<walkers.size(); iw++)
    logPsi[iw] -= SumHost[iw];
//     fprintf (stderr, "host = %25.16f\n", host_sum);
  // fprintf (stderr, "cuda = %25.16f\n", logPsi[10]);
}

void
TwoBodyJastrowOrbitalBspline::update (std::vector<Walker_t*> &walkers, int iat)
{
  // for (int iw=0; iw<walkers.size(); iw++)
  //   UpdateListHost[iw] = (CudaReal*)walkers[iw]->R_GPU.data();
  // UpdateListGPU = UpdateListHost;
  // two_body_update(UpdateListGPU.data(), N, iat, walkers.size());
}



// This currently does not actually compute the gradient or laplacian
void
TwoBodyJastrowOrbitalBspline::ratio
(MCWalkerConfiguration &W, int iat,
 std::vector<ValueType> &psi_ratios, std::vector<GradType>  &grad,
 std::vector<ValueType> &lapl)
{
  std::vector<Walker_t*> &walkers = W.WalkerList;
  if (SumGPU.size() < 4*walkers.size())
    SumGPU.resize(4*walkers.size());
#ifdef CPU_RATIO
  DTD_BConds<double,3,SUPERCELL_BULK> bconds;
  for (int iw=0; iw<walkers.size(); iw++)
  {
    Walker_t &walker = *(walkers[iw]);
    double sum = 0.0;
    int group2 = PtclRef.GroupID (iat);
    for (int group1=0; group1<PtclRef.groups(); group1++)
    {
      int first1 = PtclRef.first(group1);
      int last1  = PtclRef.last(group1);
      double factor = (group1 == group2) ? 0.5 : 1.0;
      int id = group1*NumGroups + group2;
      FT* func = F[id];
      for (int ptcl1=first1; ptcl1<last1; ptcl1++)
      {
        PosType disp = walkers[iw]->R[ptcl1] - walkers[iw]->R[iat];
        double dist = std::sqrt(bconds.apply(PtclRef.Lattice, disp));
        sum += factor*func->evaluate(dist);
        disp = walkers[iw]->R[ptcl1] - new_pos[iw];
        dist = std::sqrt(bconds.apply(PtclRef.Lattice, disp));
        sum -= factor*func->evaluate(dist);
      }
    }
    psi_ratios[iw] *= std::exp(-sum);
  }
#else
  int newGroup = PtclRef.GroupID[iat];
  bool zero = true;
  for (int group=0; group<PtclRef.groups(); group++)
  {
    int first = PtclRef.first(group);
    int last  = PtclRef.last(group) -1;
    CudaSpline<CudaReal> &spline = *(GPUSplines[group*NumGroups+newGroup]);
    // two_body_ratio (W.RList_GPU.data(), first, last,
    // 		      (CudaReal*)W.Rnew_GPU.data(), iat,
    // 		      spline.coefs.data(), spline.coefs.size(),
    // 		      spline.rMax, L.data(), Linv.data(),
    // 		      SumGPU.data(), walkers.size());
    if (UsePBC)
    {
      bool use_fast_image = W.Lattice.SimulationCellRadius >= spline.rMax;
      two_body_ratio_grad_PBC (W.RList_GPU.data(), first, last,
                               (CudaReal*)W.Rnew_GPU.data(), iat,
                               spline.coefs.data(), spline.coefs.size(),
                               spline.rMax, L.data(), Linv.data(), zero,
                               SumGPU.data(), walkers.size(), use_fast_image);
    }
    else
      two_body_ratio_grad (W.RList_GPU.data(), first, last,
                           (CudaReal*)W.Rnew_GPU.data(), iat,
                           spline.coefs.data(), spline.coefs.size(),
                           spline.rMax, zero, SumGPU.data(),
                           walkers.size());
    zero = false;
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
#endif
}

// This currently does not actually compute the gradient or laplacian
void
TwoBodyJastrowOrbitalBspline::calcRatio
(MCWalkerConfiguration &W, int iat,
 std::vector<ValueType> &psi_ratios, std::vector<GradType>  &grad,
 std::vector<ValueType> &lapl)
{
  std::vector<Walker_t*> &walkers = W.WalkerList;
  if (SumGPU.size() < 4*walkers.size())
    SumGPU.resize(4*walkers.size());
#ifdef CPU_RATIO
  DTD_BConds<double,3,SUPERCELL_BULK> bconds;
  for (int iw=0; iw<walkers.size(); iw++)
  {
    Walker_t &walker = *(walkers[iw]);
    double sum = 0.0;
    int group2 = PtclRef.GroupID (iat);
    for (int group1=0; group1<PtclRef.groups(); group1++)
    {
      int first1 = PtclRef.first(group1);
      int last1  = PtclRef.last(group1);
      double factor = (group1 == group2) ? 0.5 : 1.0;
      int id = group1*NumGroups + group2;
      FT* func = F[id];
      for (int ptcl1=first1; ptcl1<last1; ptcl1++)
      {
        PosType disp = walkers[iw]->R[ptcl1] - walkers[iw]->R[iat];
        double dist = std::sqrt(bconds.apply(PtclRef.Lattice, disp));
        sum += factor*func->evaluate(dist);
        disp = walkers[iw]->R[ptcl1] - new_pos[iw];
        dist = std::sqrt(bconds.apply(PtclRef.Lattice, disp));
        sum -= factor*func->evaluate(dist);
      }
    }
    psi_ratios[iw] *= std::exp(-sum);
  }
#else
  int newGroup = PtclRef.GroupID[iat];
  bool zero = true;
  for (int group=0; group<PtclRef.groups(); group++)
  {
    int first = PtclRef.first(group);
    int last  = PtclRef.last(group) -1;
    CudaSpline<CudaReal> &spline = *(GPUSplines[group*NumGroups+newGroup]);
    if (UsePBC)
    {
      bool use_fast_image = W.Lattice.SimulationCellRadius >= spline.rMax;
      two_body_ratio_grad_PBC (W.RList_GPU.data(), first, last,
                               (CudaReal*)W.Rnew_GPU.data(), iat,
                               spline.coefs.data(), spline.coefs.size(),
                               spline.rMax, L.data(), Linv.data(), zero,
                               SumGPU.data(), walkers.size(), use_fast_image);
    }
    else
      two_body_ratio_grad (W.RList_GPU.data(), first, last,
                           (CudaReal*)W.Rnew_GPU.data(), iat,
                           spline.coefs.data(), spline.coefs.size(),
                           spline.rMax, zero, SumGPU.data(),
                           walkers.size());
    zero = false;
  }
  // Copy data back to CPU memory
  gpu::streamsSynchronize();
  SumHost.asyncCopy(SumGPU);
  cudaEventRecord(gpu::ratioSyncTwoBodyEvent, gpu::memoryStream);
#endif
}
void
TwoBodyJastrowOrbitalBspline::addRatio
(MCWalkerConfiguration &W, int iat,
 std::vector<ValueType> &psi_ratios, std::vector<GradType>  &grad,
 std::vector<ValueType> &lapl)
{
#ifndef CPU_RATIO
  std::vector<Walker_t*> &walkers = W.WalkerList;
  cudaEventSynchronize(gpu::ratioSyncTwoBodyEvent);
  for (int iw=0; iw<walkers.size(); iw++)
  {
    psi_ratios[iw] *= std::exp(-SumHost[4*iw+0]);
    grad[iw][0] -= SumHost[4*iw+1];
    grad[iw][1] -= SumHost[4*iw+2];
    grad[iw][2] -= SumHost[4*iw+3];
  }
#endif
}




void
TwoBodyJastrowOrbitalBspline::NLratios
(MCWalkerConfiguration &W,  std::vector<NLjob> &jobList,
 std::vector<PosType> &quadPoints, std::vector<ValueType> &psi_ratios)
{
  CudaReal sim_cell_radius = W.Lattice.SimulationCellRadius;
  std::vector<Walker_t*> &walkers = W.WalkerList;
  int njobs = jobList.size();
  if (NL_JobListHost.size() < njobs)
  {
    NL_JobListHost.resize        (njobs);
    NL_SplineCoefsListHost.resize(njobs);
    NL_NumCoefsHost.resize       (njobs);
    NL_rMaxHost.resize           (njobs);
  }
  if (NL_RatiosHost.size() < quadPoints.size())
  {
    int nq = quadPoints.size();
    NL_QuadPointsHost.resize(OHMMS_DIM*nq);
    NL_QuadPointsGPU.resize (OHMMS_DIM*nq, 1.25);
    NL_RatiosHost.resize    (nq);
    NL_RatiosGPU.resize     (nq, 1.25);
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
  for (int group=0; group<PtclRef.groups(); group++)
  {
    int first = PtclRef.first(group);
    int last  = PtclRef.last(group) -1;
    for (int ijob=0; ijob<jobList.size(); ijob++)
    {
      int newGroup = PtclRef.GroupID[jobList[ijob].elec];
      CudaSpline<CudaReal> &spline = *(GPUSplines[group*NumGroups+newGroup]);
      NL_SplineCoefsListHost[ijob] = spline.coefs.data();
      NL_NumCoefsHost[ijob] = spline.coefs.size();
      NL_rMaxHost[ijob]     = spline.rMax;
    }
    NL_SplineCoefsListGPU = NL_SplineCoefsListHost;
    NL_NumCoefsGPU        = NL_NumCoefsHost;
    NL_rMaxGPU            = NL_rMaxHost;
    if (UsePBC)
      two_body_NLratios_PBC(NL_JobListGPU.data(), first, last,
                            NL_SplineCoefsListGPU.data(), NL_NumCoefsGPU.data(),
                            NL_rMaxGPU.data(), L.data(), Linv.data(),
                            sim_cell_radius, njobs);
    else
      two_body_NLratios(NL_JobListGPU.data(), first, last,
                        NL_SplineCoefsListGPU.data(), NL_NumCoefsGPU.data(),
                        NL_rMaxGPU.data(), njobs);
  }
  NL_RatiosHost = NL_RatiosGPU;
  for (int i=0; i < psi_ratios.size(); i++)
    psi_ratios[i] *= NL_RatiosHost[i];
}


void TwoBodyJastrowOrbitalBspline::calcGradient(MCWalkerConfiguration &W, int iat,
    std::vector<GradType> &grad)
{
  CudaReal sim_cell_radius = W.Lattice.SimulationCellRadius;
  std::vector<Walker_t*> &walkers = W.WalkerList;
  int newGroup = PtclRef.GroupID[iat];
  if (OneGradHost.size() < OHMMS_DIM*walkers.size())
  {
    OneGradHost.resize (walkers.size()*OHMMS_DIM);
    OneGradGPU.resize (walkers.size()*OHMMS_DIM);
  }
  for (int group=0; group<PtclRef.groups(); group++)
  {
    int first = PtclRef.first(group);
    int last  = PtclRef.last(group) -1;
    CudaSpline<CudaReal> &spline = *(GPUSplines[group*NumGroups+newGroup]);
    if (UsePBC)
      two_body_gradient_PBC (W.RList_GPU.data(), first, last, iat,
                             spline.coefs.data(), spline.coefs.size(),
                             spline.rMax, L.data(), Linv.data(), sim_cell_radius,
                             group==0, OneGradGPU.data(), walkers.size());
    else
      two_body_gradient (W.RList_GPU.data(), first, last, iat,
                         spline.coefs.data(), spline.coefs.size(),
                         spline.rMax, group==0, OneGradGPU.data(), walkers.size());
  }
  // Copy data back to CPU memory
  gpu::streamsSynchronize();
  OneGradHost.asyncCopy(OneGradGPU);
  cudaEventRecord(gpu::gradientSyncTwoBodyEvent, gpu::memoryStream);
}

void TwoBodyJastrowOrbitalBspline::addGradient(MCWalkerConfiguration &W, int iat,
    std::vector<GradType> &grad)
{
  std::vector<Walker_t*> &walkers = W.WalkerList;
  cudaEventSynchronize(gpu::gradientSyncTwoBodyEvent);
  for (int iw=0; iw<walkers.size(); iw++)
    for (int dim=0; dim<OHMMS_DIM; dim++)
      grad[iw][dim] -= OneGradHost[OHMMS_DIM*iw+dim];
}

void
TwoBodyJastrowOrbitalBspline::gradLapl (MCWalkerConfiguration &W,
                                        GradMatrix_t &grad,
                                        ValueMatrix_t &lapl)
{
  CudaReal sim_cell_radius = W.Lattice.SimulationCellRadius;
  std::vector<Walker_t*> &walkers = W.WalkerList;
  int numGL = 4*N*walkers.size();
  if (GradLaplGPU.size()  < numGL)
  {
    GradLaplGPU.resize(numGL);
    GradLaplHost.resize(numGL);
  }
  for (int iw=0; iw<walkers.size(); iw++)
  {
    Walker_t &walker = *(walkers[iw]);
    SumHost[iw] = 0.0;
  }
  SumGPU = SumHost;
  for (int i=0; i<walkers.size()*4*N; i++)
    GradLaplHost[i] = 0.0;
  GradLaplGPU = GradLaplHost;
#ifdef CUDA_DEBUG
  std::vector<CudaReal> CPU_GradLapl(4*N);
  DTD_BConds<double,3,SUPERCELL_BULK> bconds;
  int iw = 0;
  for (int group1=0; group1<PtclRef.groups(); group1++)
  {
    int first1 = PtclRef.first(group1);
    int last1  = PtclRef.last(group1) -1;
    for (int ptcl1=first1; ptcl1<=last1; ptcl1++)
    {
      PosType grad(0.0, 0.0, 0.0);
      double lapl(0.0);
      for (int group2=0; group2<PtclRef.groups(); group2++)
      {
        int first2 = PtclRef.first(group2);
        int last2  = PtclRef.last(group2) -1;
        int id = group2*NumGroups + group1;
        FT* func = F[id];
        for (int ptcl2=first2; ptcl2<=last2; ptcl2++)
        {
          if (ptcl1 != ptcl2 )
          {
            PosType disp = walkers[iw]->R[ptcl2] - walkers[iw]->R[ptcl1];
            double dist = std::sqrt(bconds.apply(PtclRef.Lattice, disp));
            double u, du, d2u;
            u = func->evaluate(dist, du, d2u);
            du /= dist;
            grad += disp * du;
            lapl += d2u + 2.0*du;
          }
        }
      }
      CPU_GradLapl[ptcl1*4+0] = grad[0];
      CPU_GradLapl[ptcl1*4+1] = grad[1];
      CPU_GradLapl[ptcl1*4+2] = grad[2];
      CPU_GradLapl[ptcl1*4+3] = lapl;
    }
  }
#endif
  for (int group1=0; group1<PtclRef.groups(); group1++)
  {
    int first1 = PtclRef.first(group1);
    int last1  = PtclRef.last(group1) -1;
    for (int group2=0; group2<PtclRef.groups(); group2++)
    {
      int first2 = PtclRef.first(group2);
      int last2  = PtclRef.last(group2) -1;
      CudaSpline<CudaReal> &spline = *(GPUSplines[group1*NumGroups+group2]);
      if (UsePBC)
        two_body_grad_lapl_PBC (W.RList_GPU.data(), first1, last1, first2, last2,
                                spline.coefs.data(), spline.coefs.size(),
                                spline.rMax, L.data(), Linv.data(), sim_cell_radius,
                                GradLaplGPU.data(), 4*N, walkers.size());
      else
        two_body_grad_lapl (W.RList_GPU.data(), first1, last1, first2, last2,
                            spline.coefs.data(), spline.coefs.size(),
                            spline.rMax, GradLaplGPU.data(), 4*N, walkers.size());
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
        char buff[500];
        gethostname(buff, 500);
        fprintf (stderr, "NAN in TwoBodyJastrowOrbitalBspline laplacian.  Host=%s\n", buff);
        abort();
      }
      lapl(iw,ptcl) += GradLaplHost[4*N*iw+ + 4*ptcl +3];
    }
  }
}


void
TwoBodyJastrowOrbitalBspline::resetParameters(const opt_variables_type& active)
{
  TwoBodyJastrowOrbital<BsplineFunctor<OrbitalBase::RealType> >::resetParameters(active);
  for (int i=0; i<NumGroups*NumGroups; i++)
    GPUSplines[i]->set(*F[i]);
}

void
TwoBodyJastrowOrbitalBspline::evaluateDerivatives
(MCWalkerConfiguration &W, const opt_variables_type& optvars,
 RealMatrix_t &d_logpsi, RealMatrix_t &dlapl_over_psi)
{
  CudaReal sim_cell_radius = W.Lattice.SimulationCellRadius;
  std::vector<Walker_t*> &walkers = W.WalkerList;
  int nw = walkers.size();
  for (int i=0; i<UniqueSplines.size(); i++)
    if (DerivListGPU.size() < nw)
    {
      MaxCoefs = 0;
      for (int i=0; i<UniqueSplines.size(); i++)
        MaxCoefs = std::max(MaxCoefs, (int)UniqueSplines[i]->coefs.size());
      // Round up to nearest 16 to allow coallesced GPU reads
      MaxCoefs = ((MaxCoefs+7)/8)*8;
      SplineDerivsHost.resize(2*MaxCoefs*nw);
      SplineDerivsGPU.resize(2*MaxCoefs*nw);
      DerivListHost.resize(nw);
      DerivListGPU.resize(nw);
      for (int iw=0; iw<nw; iw++)
        DerivListHost[iw] = SplineDerivsGPU.data() + 2*iw*MaxCoefs;
      DerivListGPU = DerivListHost;
    }
  std::vector<TinyVector<RealType, 2> > derivs;
  for (int group1=0; group1<PtclRef.groups(); group1++)
  {
    int first1 = PtclRef.first(group1);
    int last1  = PtclRef.last(group1) -1;
    for (int group2=0; group2<PtclRef.groups(); group2++)
    {
      int first2 = PtclRef.first(group2);
      int last2  = PtclRef.last(group2) -1;
      int ptype = group1*NumGroups+group2;
      CudaSpline<CudaReal> &spline = *(GPUSplines[group1*NumGroups+group2]);
      if (UsePBC)
        two_body_derivs_PBC (W.RList_GPU.data(), W.GradList_GPU.data(),
                             first1, last1, first2, last2,
                             spline.coefs.size(), spline.rMax, L.data(),
                             Linv.data(), sim_cell_radius, DerivListGPU.data(),nw);
      else
        two_body_derivs (W.RList_GPU.data(), W.GradList_GPU.data(),
                         first1, last1, first2, last2,
                         spline.coefs.size(), spline.rMax, DerivListGPU.data(),nw);
      // Copy data back to CPU memory
      SplineDerivsHost = SplineDerivsGPU;
      opt_variables_type splineVars = F[ptype]->myVars;
      for (int iv=0; iv<splineVars.size(); iv++)
      {
        // std::cerr << "groups = (" << group1 << "," << group2
        //      << ") Index=" << splineVars.Index[iv] << std::endl;
        int varIndex = splineVars.Index[iv];
        int coefIndex = iv+1;
        for (int iw=0; iw<nw; iw++)
        {
          d_logpsi(iw,varIndex) +=
            SplineDerivsHost[2*(MaxCoefs*iw+coefIndex)+0];
          dlapl_over_psi(iw,varIndex) +=
            SplineDerivsHost[2*(MaxCoefs*iw+coefIndex)+1];
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



}

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
#include "QMCWaveFunctions/detail/CUDA_legacy/BsplineJastrowCuda.h"
#include "QMCWaveFunctions/detail/CUDA_legacy/BsplineJastrowCudaPBC.h"
#include <unistd.h>


namespace qmcplusplus
{
template<class FT>
void TwoBodyJastrowOrbitalBspline<FT>::freeGPUmem()
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
};

template<class FT>
void TwoBodyJastrowOrbitalBspline<FT>::recompute(MCWalkerConfiguration& W, bool firstTime)
{}

template<class FT>
void TwoBodyJastrowOrbitalBspline<FT>::reserve(PointerPool<gpu::device_vector<CTS::RealType>>& pool)
{}

template<class FT>
void TwoBodyJastrowOrbitalBspline<FT>::checkInVariables(opt_variables_type& active)
{
  J2OrbitalSoA<BsplineFunctor<WaveFunctionComponent::RealType>>::checkInVariables(active);
  for (int i = 0; i < this->NumGroups * this->NumGroups; i++)
    GPUSplines[i]->set(*this->F[i]);
}

template<class FT>
void TwoBodyJastrowOrbitalBspline<FT>::addFunc(int ia, int ib, std::unique_ptr<FT> j)
{
  CudaSpline<CTS::RealType>* newSpline = new CudaSpline<CTS::RealType>(*j);
  J2OrbitalSoA<BsplineFunctor<WaveFunctionComponent::RealType>>::addFunc(ia, ib, std::move(j));
  UniqueSplines.push_back(newSpline);
  if (ia == ib)
  {
    if (ia == 0) //first time, assign everything
    {
      int ij = 0;
      for (int ig = 0; ig < this->NumGroups; ++ig)
        for (int jg = 0; jg < this->NumGroups; ++jg, ++ij)
          if (GPUSplines[ij] == 0)
            GPUSplines[ij] = newSpline;
    }
    else
      GPUSplines[ia * this->NumGroups + ib] = newSpline;
  }
  else
  {
    // a very special case, 1 particle of each type (e.g. 1 up + 1 down)
    // uu/dd/etc. was prevented by the builder
    if (PtclRef.R.size() == this->NumGroups)
      for (int ig = 0; ig < this->NumGroups; ++ig)
        GPUSplines[ig * this->NumGroups + ig] = newSpline;
    // generic case
    GPUSplines[ia * this->NumGroups + ib] = newSpline;
    GPUSplines[ib * this->NumGroups + ia] = newSpline;
  }
}

template<class FT>
void TwoBodyJastrowOrbitalBspline<FT>::addLog(MCWalkerConfiguration& W, std::vector<RealType>& logPsi)
{
  auto& walkers = W.WalkerList;
  if (SumGPU.size() < 4 * walkers.size())
  {
    SumGPU.resize(4 * walkers.size());
    SumHost.resize(4 * walkers.size());
    OneGradGPU.resize(walkers.size() * OHMMS_DIM);
    OneGradHost.resize(walkers.size() * OHMMS_DIM);
    UpdateListHost.resize(walkers.size());
    UpdateListGPU.resize(walkers.size());
  }
  int numGL = 4 * this->N * walkers.size();
  if (GradLaplGPU.size() < numGL)
  {
    GradLaplGPU.resize(numGL);
    GradLaplHost.resize(numGL);
  }
  for (int iw = 0; iw < walkers.size(); iw++)
  {
    Walker_t& walker    = *(walkers[iw]);
    SumHost[iw]         = 0.0;
    CTS::RealType* dest = (CTS::RealType*)walker.R_GPU.data();
  }
  SumGPU = SumHost;
  //     DTD_BConds<double,3,SUPERCELL_BULK> bconds;
  //     double host_sum = 0.0;
  for (int group1 = 0; group1 < PtclRef.groups(); group1++)
  {
    int first1 = PtclRef.first(group1);
    int last1  = PtclRef.last(group1) - 1;
    for (int group2 = group1; group2 < PtclRef.groups(); group2++)
    {
      int first2 = PtclRef.first(group2);
      int last2  = PtclRef.last(group2) - 1;
      // 	double factor = (group1 == group2) ? 0.5 : 1.0;
      // 	for (int e1=first1; e1 <= last1; e1++)
      // 	  for (int e2=first2; e2 <= last2; e2++) {
      // 	    PosType disp = walkers[0]->R[e2] - walkers[0]->R[e1];
      // 	    double dist = std::sqrt(bconds.apply(PtclRef.Lattice, disp));
      // 	    if (e1 != e2)
      // 	      host_sum -= factor * F[group2*NumGroups + group1]->evaluate(dist);
      // 	  }
      CudaSpline<CTS::RealType>& spline = *(GPUSplines[group1 * this->NumGroups + group2]);
      if (UsePBC)
        two_body_sum_PBC(W.RList_GPU.data(), first1, last1, first2, last2, spline.coefs.data(), spline.coefs.size(),
                         spline.rMax, L.data(), Linv.data(), SumGPU.data(), walkers.size());
      else
        two_body_sum(W.RList_GPU.data(), first1, last1, first2, last2, spline.coefs.data(), spline.coefs.size(),
                     spline.rMax, SumGPU.data(), walkers.size());
    }
  }
  // Copy data back to CPU memory
  SumHost = SumGPU;
  for (int iw = 0; iw < walkers.size(); iw++)
    logPsi[iw] -= SumHost[iw];
  //     fprintf (stderr, "host = %25.16f\n", host_sum);
  // fprintf (stderr, "cuda = %25.16f\n", logPsi[10]);
}

template<class FT>
void TwoBodyJastrowOrbitalBspline<FT>::update(MCWalkerConfiguration* W,
                                              std::vector<Walker_t*>& walkers,
                                              int iat,
                                              std::vector<bool>* acc,
                                              int k)
{
  // for (int iw=0; iw<walkers.size(); iw++)
  //   UpdateListHost[iw] = (CTS::RealType*)walkers[iw]->R_GPU.data();
  // UpdateListGPU = UpdateListHost;
  // two_body_update(UpdateListGPU.data(), N, iat, walkers.size());
}

// #define DEBUG_DELAYED

// This currently does not actually compute the gradient or laplacian
template<class FT>
void TwoBodyJastrowOrbitalBspline<FT>::ratio(MCWalkerConfiguration& W,
                                             int iat,
                                             std::vector<ValueType>& psi_ratios,
                                             std::vector<GradType>& grad,
                                             std::vector<ValueType>& lapl)
{
  auto& walkers = W.WalkerList;
  int N         = W.Rnew_GPU.size();
  int nw        = walkers.size();
  if (SumGPU.size() < 4 * nw)
    SumGPU.resize(4 * nw);
#ifdef CPU_RATIO
  DTD_BConds<double, 3, SUPERCELL_BULK> bconds;
  for (int iw = 0; iw < walkers.size(); iw++)
  {
    Walker_t& walker = *(walkers[iw]);
    double sum       = 0.0;
    int group2       = PtclRef.GroupID(iat);
    for (int group1 = 0; group1 < PtclRef.groups(); group1++)
    {
      int first1    = PtclRef.first(group1);
      int last1     = PtclRef.last(group1);
      double factor = (group1 == group2) ? 0.5 : 1.0;
      int id        = group1 * this->NumGroups + group2;
      FT* func      = F[id];
      for (int ptcl1 = first1; ptcl1 < last1; ptcl1++)
      {
        PosType disp = walkers[iw]->R[ptcl1] - walkers[iw]->R[iat];
        double dist  = std::sqrt(bconds.apply(PtclRef.Lattice, disp));
        sum += factor * func->evaluate(dist);
        disp = walkers[iw]->R[ptcl1] - new_pos[iw];
        dist = std::sqrt(bconds.apply(PtclRef.Lattice, disp));
        sum -= factor * func->evaluate(dist);
      }
    }
    psi_ratios[iw] *= std::exp(-sum);
  }
#else
  int newGroup = PtclRef.GroupID[iat];
  bool zero    = true;
  for (int group = 0; group < PtclRef.groups(); group++)
  {
    int first                         = PtclRef.first(group);
    int last                          = PtclRef.last(group) - 1;
    CudaSpline<CTS::RealType>& spline = *(GPUSplines[group * this->NumGroups + newGroup]);
    // two_body_ratio (W.RList_GPU.data(), first, last,
    // 		      (CTS::RealType*)W.Rnew_GPU.data(), iat,
    // 		      spline.coefs.data(), spline.coefs.size(),
    // 		      spline.rMax, L.data(), Linv.data(),
    // 		      SumGPU.data(), walkers.size());
    if (UsePBC)
    {
      bool use_fast_image = W.getLattice().SimulationCellRadius >= spline.rMax;
      two_body_ratio_grad_PBC(W.RList_GPU.data(), first, last, (CTS::RealType*)W.Rnew_GPU.data(), iat, kcurr * nw,
                              spline.coefs.data(), spline.coefs.size(), spline.rMax, L.data(), Linv.data(), zero,
                              SumGPU.data(), nw, use_fast_image);
    }
    else
      two_body_ratio_grad(W.RList_GPU.data(), first, last, (CTS::RealType*)W.Rnew_GPU.data(), iat, kcurr * nw,
                          spline.coefs.data(), spline.coefs.size(), spline.rMax, zero, SumGPU.data(), nw);
    zero = false;
  }
  // Copy data back to CPU memory
  SumHost = SumGPU;
  for (int iw = 0; iw < nw; iw++)
  {
#ifdef DEBUG_DELAYED
    if (iw % nw == 0)
      fprintf(stderr, "k = %i:\n", kcurr);
    fprintf(stderr, "walker %i ratio 2B Jastrow: %f|(%f,%f,%f) -> ", iw % nw, psi_ratios[iw], grad[iw][0], grad[iw][1],
            grad[iw][2]);
#endif
    psi_ratios[nw * kcurr + iw] *= std::exp(-SumHost[4 * iw + 0]);
    grad[nw * kcurr + iw][0] -= SumHost[4 * iw + 1];
    grad[nw * kcurr + iw][1] -= SumHost[4 * iw + 2];
    grad[nw * kcurr + iw][2] -= SumHost[4 * iw + 3];
#ifdef DEBUG_DELAYED
    fprintf(stderr, "%f|(%f,%f,%f)\n", psi_ratios[nw * kcurr + iw], grad[iw][0], grad[iw][1], grad[iw][2]);
#endif
  }
#endif
}

// This currently does not actually compute the gradient or laplacian
template<class FT>
void TwoBodyJastrowOrbitalBspline<FT>::calcRatio(MCWalkerConfiguration& W,
                                                 int iat,
                                                 std::vector<ValueType>& psi_ratios,
                                                 std::vector<GradType>& grad,
                                                 std::vector<ValueType>& lapl)
{
  auto& walkers = W.WalkerList;
  int N         = W.Rnew_GPU.size();
  int nw        = walkers.size();
  int kd        = W.getkDelay();
  int k         = W.getkcurr() - (kd > 1);
  if (k < 0)
    k += W.getkupdate();
  int offset = 0;
  if (W.getklinear())
    offset = k * nw;
  if (SumGPU.size() < 4 * nw)
    SumGPU.resize(4 * nw);
#ifdef CPU_RATIO
  DTD_BConds<double, 3, SUPERCELL_BULK> bconds;
  for (int iw = 0; iw < walkers.size(); iw++)
  {
    Walker_t& walker = *(walkers[iw]);
    double sum       = 0.0;
    int group2       = PtclRef.GroupID(iat);
    for (int group1 = 0; group1 < PtclRef.groups(); group1++)
    {
      int first1    = PtclRef.first(group1);
      int last1     = PtclRef.last(group1);
      double factor = (group1 == group2) ? 0.5 : 1.0;
      int id        = group1 * this->NumGroups + group2;
      FT* func      = F[id];
      for (int ptcl1 = first1; ptcl1 < last1; ptcl1++)
      {
        PosType disp = walkers[iw]->R[ptcl1] - walkers[iw]->R[iat];
        double dist  = std::sqrt(bconds.apply(PtclRef.Lattice, disp));
        sum += factor * func->evaluate(dist);
        disp = walkers[iw]->R[ptcl1] - new_pos[iw];
        dist = std::sqrt(bconds.apply(PtclRef.Lattice, disp));
        sum -= factor * func->evaluate(dist);
      }
    }
    psi_ratios[iw] *= std::exp(-sum);
  }
#else
  int newGroup = PtclRef.GroupID[iat];
  bool zero    = true;
  for (int group = 0; group < PtclRef.groups(); group++)
  {
    int first                         = PtclRef.first(group);
    int last                          = PtclRef.last(group) - 1;
    CudaSpline<CTS::RealType>& spline = *(GPUSplines[group * this->NumGroups + newGroup]);
    if (UsePBC)
    {
      bool use_fast_image = W.getLattice().SimulationCellRadius >= spline.rMax;
      two_body_ratio_grad_PBC(W.RList_GPU.data(), first, last, &(((CTS::RealType*)W.Rnew_GPU.data())[3 * offset]), iat,
                              kcurr * nw, spline.coefs.data(), spline.coefs.size(), spline.rMax, L.data(), Linv.data(),
                              zero, SumGPU.data(), nw, use_fast_image);
    }
    else
      two_body_ratio_grad(W.RList_GPU.data(), first, last, &(((CTS::RealType*)W.Rnew_GPU.data())[3 * offset]), iat,
                          kcurr * nw, spline.coefs.data(), spline.coefs.size(), spline.rMax, zero, SumGPU.data(), nw);
    zero = false;
  }
  // Copy data back to CPU memory
  gpu::streamsSynchronize();
  SumHost.asyncCopy(SumGPU);
  cudaEventRecord(gpu::ratioSyncTwoBodyEvent, gpu::memoryStream);
#endif
}

template<class FT>
void TwoBodyJastrowOrbitalBspline<FT>::addRatio(MCWalkerConfiguration& W,
                                                int iat,
                                                int k,
                                                std::vector<ValueType>& psi_ratios,
                                                std::vector<GradType>& grad,
                                                std::vector<ValueType>& lapl)
{
#ifndef CPU_RATIO
  auto& walkers = W.WalkerList;
  cudaEventSynchronize(gpu::ratioSyncTwoBodyEvent);
  for (int iw = 0; iw < walkers.size(); iw++)
  {
#ifdef DEBUG_DELAYED
    fprintf(stderr, "-> 2B Jastrow walker %i: ratio = %f ; grad = (%f,%f,%f) -> ", iw, psi_ratios[iw], grad[iw][0],
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
#endif
}

template<class FT>
void TwoBodyJastrowOrbitalBspline<FT>::NLratios(MCWalkerConfiguration& W,
                                                std::vector<NLjob>& jobList,
                                                std::vector<PosType>& quadPoints,
                                                std::vector<ValueType>& psi_ratios)
{
  CTS::RealType sim_cell_radius = W.getLattice().SimulationCellRadius;
  auto& walkers                 = W.WalkerList;
  int njobs                     = jobList.size();
  if (NL_JobListHost.size() < njobs)
  {
    NL_JobListHost.resize(njobs);
    NL_SplineCoefsListHost.resize(njobs);
    NL_NumCoefsHost.resize(njobs);
    NL_rMaxHost.resize(njobs);
  }
  if (NL_RatiosHost.size() < quadPoints.size())
  {
    int nq = quadPoints.size();
    NL_QuadPointsHost.resize(OHMMS_DIM * nq);
    NL_QuadPointsGPU.resize(OHMMS_DIM * nq, 1.25);
    NL_RatiosHost.resize(nq);
    NL_RatiosGPU.resize(nq, 1.25);
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
  for (int group = 0; group < PtclRef.groups(); group++)
  {
    int first = PtclRef.first(group);
    int last  = PtclRef.last(group) - 1;
    for (int ijob = 0; ijob < jobList.size(); ijob++)
    {
      int newGroup                      = PtclRef.GroupID[jobList[ijob].elec];
      CudaSpline<CTS::RealType>& spline = *(GPUSplines[group * this->NumGroups + newGroup]);
      NL_SplineCoefsListHost[ijob]      = spline.coefs.data();
      NL_NumCoefsHost[ijob]             = spline.coefs.size();
      NL_rMaxHost[ijob]                 = spline.rMax;
    }
    NL_SplineCoefsListGPU = NL_SplineCoefsListHost;
    NL_NumCoefsGPU        = NL_NumCoefsHost;
    NL_rMaxGPU            = NL_rMaxHost;
    if (UsePBC)
      two_body_NLratios_PBC(NL_JobListGPU.data(), first, last, NL_SplineCoefsListGPU.data(), NL_NumCoefsGPU.data(),
                            NL_rMaxGPU.data(), L.data(), Linv.data(), sim_cell_radius, njobs);
    else
      two_body_NLratios(NL_JobListGPU.data(), first, last, NL_SplineCoefsListGPU.data(), NL_NumCoefsGPU.data(),
                        NL_rMaxGPU.data(), njobs);
  }
  NL_RatiosHost = NL_RatiosGPU;
  for (int i = 0; i < psi_ratios.size(); i++)
    psi_ratios[i] *= NL_RatiosHost[i];
}

template<class FT>
void TwoBodyJastrowOrbitalBspline<FT>::calcGradient(MCWalkerConfiguration& W,
                                                    int iat,
                                                    int k,
                                                    std::vector<GradType>& grad)
{
  CTS::RealType sim_cell_radius = W.getLattice().SimulationCellRadius;
  auto& walkers                 = W.WalkerList;
  int newGroup                  = PtclRef.GroupID[iat];
  if (OneGradHost.size() < OHMMS_DIM * walkers.size())
  {
    OneGradHost.resize(walkers.size() * OHMMS_DIM);
    OneGradGPU.resize(walkers.size() * OHMMS_DIM);
  }
  for (int group = 0; group < PtclRef.groups(); group++)
  {
    int first                         = PtclRef.first(group);
    int last                          = PtclRef.last(group) - 1;
    CudaSpline<CTS::RealType>& spline = *(GPUSplines[group * this->NumGroups + newGroup]);
    if (UsePBC)
      two_body_gradient_PBC(W.RList_GPU.data(), first, last, iat, spline.coefs.data(), spline.coefs.size(), spline.rMax,
                            L.data(), Linv.data(), sim_cell_radius, group == 0, OneGradGPU.data(), walkers.size());
    else
      two_body_gradient(W.RList_GPU.data(), first, last, iat, spline.coefs.data(), spline.coefs.size(), spline.rMax,
                        group == 0, OneGradGPU.data(), walkers.size());
  }
  // Copy data back to CPU memory
  gpu::streamsSynchronize();
  OneGradHost.asyncCopy(OneGradGPU);
  cudaEventRecord(gpu::gradientSyncTwoBodyEvent, gpu::memoryStream);
}

template<class FT>
void TwoBodyJastrowOrbitalBspline<FT>::addGradient(MCWalkerConfiguration& W, int iat, std::vector<GradType>& grad)
{
  auto& walkers = W.WalkerList;
  cudaEventSynchronize(gpu::gradientSyncTwoBodyEvent);
  //  for (int iw=0; iw<walkers.size(); iw++)
  //    for (int dim=0; dim<OHMMS_DIM; dim++)
  //      grad[iw][dim] -= OneGradHost[OHMMS_DIM*iw+dim];
  for (int iw = 0; iw < walkers.size(); iw++)
  {
#ifdef DEBUG_DELAYED
    fprintf(stderr, "2B Jastrow grad walker %i: (", iw);
#endif
    for (int dim = 0; dim < OHMMS_DIM; dim++)
    {
#ifdef DEBUG_DELAYED
      if (dim > 0)
        fprintf(stderr, ", ");
      fprintf(stderr, "%f (before: %f)", OneGradHost[OHMMS_DIM * iw + dim], grad[iw][dim]);
#endif
      grad[iw][dim] -= OneGradHost[OHMMS_DIM * iw + dim];
    }
#ifdef DEBUG_DELAYED
    fprintf(stderr, ")\n");
#endif
  }
}

template<class FT>
void TwoBodyJastrowOrbitalBspline<FT>::gradLapl(MCWalkerConfiguration& W, GradMatrix& grad, ValueMatrix& lapl)
{
  CTS::RealType sim_cell_radius = W.getLattice().SimulationCellRadius;
  auto& walkers                 = W.WalkerList;
  int numGL                     = 4 * this->N * walkers.size();
  if (GradLaplGPU.size() < numGL)
  {
    GradLaplGPU.resize(numGL);
    GradLaplHost.resize(numGL);
  }
  for (int iw = 0; iw < walkers.size(); iw++)
  {
    Walker_t& walker = *(walkers[iw]);
    SumHost[iw]      = 0.0;
  }
  SumGPU = SumHost;
  for (int i = 0; i < walkers.size() * 4 * this->N; i++)
    GradLaplHost[i] = 0.0;
  GradLaplGPU = GradLaplHost;
#ifdef CUDA_DEBUG
  std::vector<CTS::RealType> CPU_GradLapl(4 * this->N);
  DTD_BConds<double, 3, SUPERCELL_BULK> bconds;
  int iw = 0;
  for (int group1 = 0; group1 < PtclRef.groups(); group1++)
  {
    int first1 = PtclRef.first(group1);
    int last1  = PtclRef.last(group1) - 1;
    for (int ptcl1 = first1; ptcl1 <= last1; ptcl1++)
    {
      PosType grad(0.0, 0.0, 0.0);
      double lapl(0.0);
      for (int group2 = 0; group2 < PtclRef.groups(); group2++)
      {
        int first2 = PtclRef.first(group2);
        int last2  = PtclRef.last(group2) - 1;
        int id     = group2 * this->NumGroups + group1;
        FT* func   = F[id];
        for (int ptcl2 = first2; ptcl2 <= last2; ptcl2++)
        {
          if (ptcl1 != ptcl2)
          {
            PosType disp = walkers[iw]->R[ptcl2] - walkers[iw]->R[ptcl1];
            double dist  = std::sqrt(bconds.apply(PtclRef.Lattice, disp));
            double u, du, d2u;
            u = func->evaluate(dist, du, d2u);
            du /= dist;
            grad += disp * du;
            lapl += d2u + 2.0 * du;
          }
        }
      }
      CPU_GradLapl[ptcl1 * 4 + 0] = grad[0];
      CPU_GradLapl[ptcl1 * 4 + 1] = grad[1];
      CPU_GradLapl[ptcl1 * 4 + 2] = grad[2];
      CPU_GradLapl[ptcl1 * 4 + 3] = lapl;
    }
  }
#endif
  for (int group1 = 0; group1 < PtclRef.groups(); group1++)
  {
    int first1 = PtclRef.first(group1);
    int last1  = PtclRef.last(group1) - 1;
    for (int group2 = 0; group2 < PtclRef.groups(); group2++)
    {
      int first2                        = PtclRef.first(group2);
      int last2                         = PtclRef.last(group2) - 1;
      CudaSpline<CTS::RealType>& spline = *(GPUSplines[group1 * this->NumGroups + group2]);
      if (UsePBC)
        two_body_grad_lapl_PBC(W.RList_GPU.data(), first1, last1, first2, last2, spline.coefs.data(),
                               spline.coefs.size(), spline.rMax, L.data(), Linv.data(), sim_cell_radius,
                               GradLaplGPU.data(), 4 * this->N, walkers.size());
      else
        two_body_grad_lapl(W.RList_GPU.data(), first1, last1, first2, last2, spline.coefs.data(), spline.coefs.size(),
                           spline.rMax, GradLaplGPU.data(), 4 * this->N, walkers.size());
    }
  }
  // Copy data back to CPU memory
  GradLaplHost = GradLaplGPU;
#ifdef CUDA_DEBUG
  fprintf(stderr, "GPU  grad = %12.5e %12.5e %12.5e   Lapl = %12.5e\n", GradLaplHost[0], GradLaplHost[1],
          GradLaplHost[2], GradLaplHost[3]);
  fprintf(stderr, "CPU  grad = %12.5e %12.5e %12.5e   Lapl = %12.5e\n", CPU_GradLapl[0], CPU_GradLapl[1],
          CPU_GradLapl[2], CPU_GradLapl[3]);
#endif
  for (int iw = 0; iw < walkers.size(); iw++)
  {
    for (int ptcl = 0; ptcl < this->N; ptcl++)
    {
      for (int i = 0; i < OHMMS_DIM; i++)
        grad(iw, ptcl)[i] += GradLaplHost[4 * this->N * iw + 4 * ptcl + i];
      if (std::isnan(GradLaplHost[4 * this->N * iw + +4 * ptcl + 3]))
      {
        char buff[500];
        gethostname(buff, 500);
        fprintf(stderr, "NAN in TwoBodyJastrowOrbitalBspline laplacian.  Host=%s\n", buff);
        abort();
      }
      lapl(iw, ptcl) += GradLaplHost[4 * this->N * iw + +4 * ptcl + 3];
    }
  }
}

template<class FT>
void TwoBodyJastrowOrbitalBspline<FT>::resetParameters(const opt_variables_type& active)
{
  J2OrbitalSoA<BsplineFunctor<WaveFunctionComponent::RealType>>::resetParameters(active);
  for (int i = 0; i < this->NumGroups * this->NumGroups; i++)
    GPUSplines[i]->set(*this->F[i]);
}

template<class FT>
void TwoBodyJastrowOrbitalBspline<FT>::evaluateDerivatives(
    MCWalkerConfiguration& W,
    const opt_variables_type& optvars,
    TwoBodyJastrowOrbitalBspline<FT>::RealMatrix_t& d_logpsi,
    TwoBodyJastrowOrbitalBspline<FT>::RealMatrix_t& dlapl_over_psi)
{
  CTS::RealType sim_cell_radius = W.getLattice().SimulationCellRadius;
  auto& walkers                 = W.WalkerList;
  int nw                        = walkers.size();
  for (int i = 0; i < UniqueSplines.size(); i++)
    if (DerivListGPU.size() < nw)
    {
      MaxCoefs = 0;
      for (int i = 0; i < UniqueSplines.size(); i++)
        MaxCoefs = std::max(MaxCoefs, (int)UniqueSplines[i]->coefs.size());
      // Round up to nearest 16 to allow coallesced GPU reads
      MaxCoefs = ((MaxCoefs + 7) / 8) * 8;
      SplineDerivsHost.resize(2 * MaxCoefs * nw);
      SplineDerivsGPU.resize(2 * MaxCoefs * nw);
      DerivListHost.resize(nw);
      DerivListGPU.resize(nw);
      for (int iw = 0; iw < nw; iw++)
        DerivListHost[iw] = SplineDerivsGPU.data() + 2 * iw * MaxCoefs;
      DerivListGPU = DerivListHost;
    }
  std::vector<TinyVector<RealType, 2>> derivs;
  for (int group1 = 0; group1 < PtclRef.groups(); group1++)
  {
    int first1 = PtclRef.first(group1);
    int last1  = PtclRef.last(group1) - 1;
    for (int group2 = 0; group2 < PtclRef.groups(); group2++)
    {
      int first2                        = PtclRef.first(group2);
      int last2                         = PtclRef.last(group2) - 1;
      int ptype                         = group1 * this->NumGroups + group2;
      CudaSpline<CTS::RealType>& spline = *(GPUSplines[group1 * this->NumGroups + group2]);
      if (UsePBC)
        two_body_derivs_PBC(W.RList_GPU.data(), W.GradList_GPU.data(), first1, last1, first2, last2,
                            spline.coefs.size(), spline.rMax, L.data(), Linv.data(), sim_cell_radius,
                            DerivListGPU.data(), nw);
      else
        two_body_derivs(W.RList_GPU.data(), W.GradList_GPU.data(), first1, last1, first2, last2, spline.coefs.size(),
                        spline.rMax, DerivListGPU.data(), nw);
      // Copy data back to CPU memory
      SplineDerivsHost              = SplineDerivsGPU;
      opt_variables_type splineVars = this->F[ptype]->myVars;
      for (int iv = 0; iv < splineVars.size(); iv++)
      {
        // std::cerr << "groups = (" << group1 << "," << group2
        //      << ") Index=" << splineVars.Index[iv] << std::endl;
        int varIndex  = splineVars.Index[iv];
        int coefIndex = iv + 1;
        for (int iw = 0; iw < nw; iw++)
        {
          d_logpsi(iw, varIndex) += SplineDerivsHost[2 * (MaxCoefs * iw + coefIndex) + 0];
          dlapl_over_psi(iw, varIndex) += SplineDerivsHost[2 * (MaxCoefs * iw + coefIndex) + 1];
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
}

// explicit instantiations of templates
template class TwoBodyJastrowOrbitalBspline<BsplineFunctor<WaveFunctionComponent::RealType>>;


} // namespace qmcplusplus

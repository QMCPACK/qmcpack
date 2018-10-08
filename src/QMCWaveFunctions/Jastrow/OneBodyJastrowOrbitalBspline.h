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
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#ifndef ONE_BODY_JASTROW_ORBITAL_BSPLINE_H
#define ONE_BODY_JASTROW_ORBITAL_BSPLINE_H

#include "Particle/DistanceTableData.h"
#include "Particle/DistanceTable.h"
#include "QMCWaveFunctions/Jastrow/J1OrbitalSoA.h"
#include "QMCWaveFunctions/Jastrow/BsplineFunctor.h"
#include "QMCWaveFunctions/Jastrow/CudaSpline.h"
#include "QMCWaveFunctions/Jastrow/BsplineJastrowCuda.h"
#include "QMCWaveFunctions/Jastrow/BsplineJastrowCudaPBC.h"
#include "Configuration.h"

namespace qmcplusplus
{

template<class FT>
class OneBodyJastrowOrbitalBspline :
  public J1OrbitalSoA<FT>
{
private:
  bool UsePBC;
  // The following is so we can refer to type aliases(defs) below the
  // templated base class in the object hierarchy
  // Mostly QMCTraits here
  using JBase = J1OrbitalSoA<FT>;
  using CudaRealType = typename JBase::CudaRealType;
  // Duplication that should be removed
  using CudaReal = CudaRealType;
  using RealType = typename JBase::RealType;
  using ValueType = typename JBase::ValueType;
  using GradType = typename JBase::GradType;
  using PosType = typename JBase::PosType;
  using GradMatrix_t = typename JBase::GradMatrix_t;
  using ValueMatrix_t = typename JBase::ValueMatrix_t;
  using RealMatrix_t = typename JBase::RealMatrix_t;

  std::vector<CudaSpline<CudaReal>*> GPUSplines, UniqueSplines;
  int MaxCoefs;
  ParticleSet &ElecRef;
  gpu::device_vector<CudaReal> L, Linv;

  // Holds center positions
  gpu::device_vector<CudaReal> C;

  gpu::device_vector<CudaReal*> UpdateListGPU;
  gpu::device_vector<CudaReal> SumGPU, GradLaplGPU, OneGradGPU;

  gpu::host_vector<CudaReal*> UpdateListHost;
  gpu::host_vector<CudaReal> SumHost, GradLaplHost, OneGradHost;
  int NumCenterGroups, NumElecGroups;
  std::vector<int> CenterFirst, CenterLast;
  gpu::host_vector<CudaReal> SplineDerivsHost;
  gpu::device_vector<CudaReal> SplineDerivsGPU;
  gpu::host_vector<CudaReal*> DerivListHost;
  gpu::device_vector<CudaReal*> DerivListGPU;

  gpu::host_vector<CudaReal*> NL_SplineCoefsListHost;
  gpu::device_vector<CudaReal*> NL_SplineCoefsListGPU;
  gpu::host_vector<NLjobGPU<CudaReal> > NL_JobListHost;
  gpu::device_vector<NLjobGPU<CudaReal> > NL_JobListGPU;
  gpu::host_vector<int> NL_NumCoefsHost, NL_NumQuadPointsHost;
  gpu::device_vector<int> NL_NumCoefsGPU,  NL_NumQuadPointsGPU;
  gpu::host_vector<CudaReal> NL_rMaxHost, NL_QuadPointsHost, NL_RatiosHost;
  gpu::device_vector<CudaReal> NL_rMaxGPU,  NL_QuadPointsGPU,  NL_RatiosGPU;

  int N;
public:
  typedef ParticleSet::Walker_t     Walker_t;

  GPU_XRAY_TRACE void resetParameters(const opt_variables_type& active);
  GPU_XRAY_TRACE void checkInVariables(opt_variables_type& active);
  GPU_XRAY_TRACE void addFunc(int ig, FT* j, int jg=-1);
  GPU_XRAY_TRACE void recompute(MCWalkerConfiguration &W, bool firstTime);
  GPU_XRAY_TRACE void reserve (PointerPool<gpu::device_vector<CudaRealType> > &pool);
  GPU_XRAY_TRACE void addLog (MCWalkerConfiguration &W, std::vector<RealType> &logPsi);
  GPU_XRAY_TRACE void update (std::vector<Walker_t*> &walkers, int iat);
  void update (const std::vector<Walker_t*> &walkers, const std::vector<int> &iatList)
  {
    /* This function doesn't really need to return the ratio */
  }
  GPU_XRAY_TRACE void ratio (MCWalkerConfiguration &W, int iat,
              std::vector<ValueType> &psi_ratios,	std::vector<GradType>  &grad,
              std::vector<ValueType> &lapl);
  GPU_XRAY_TRACE void calcRatio (MCWalkerConfiguration &W, int iat,
                  std::vector<ValueType> &psi_ratios,	std::vector<GradType>  &grad,
                  std::vector<ValueType> &lapl);
  GPU_XRAY_TRACE void addRatio (MCWalkerConfiguration &W, int iat,
                 std::vector<ValueType> &psi_ratios,	std::vector<GradType>  &grad,
                 std::vector<ValueType> &lapl);
  GPU_XRAY_TRACE void ratio (std::vector<Walker_t*> &walkers,    std::vector<int> &iatList,
              std::vector<PosType> &rNew, std::vector<ValueType> &psi_ratios,
              std::vector<GradType>  &grad, std::vector<ValueType> &lapl)
  {
    /* This function doesn't really need to return the ratio */
  }


  GPU_XRAY_TRACE void calcGradient(MCWalkerConfiguration &W, int iat,
                    std::vector<GradType> &grad);
  GPU_XRAY_TRACE void addGradient(MCWalkerConfiguration &W, int iat,
                   std::vector<GradType> &grad);
  GPU_XRAY_TRACE void gradLapl (MCWalkerConfiguration &W, GradMatrix_t &grads,
                 ValueMatrix_t &lapl);
  GPU_XRAY_TRACE void NLratios (MCWalkerConfiguration &W,  std::vector<NLjob> &jobList,
                 std::vector<PosType> &quadPoints, std::vector<ValueType> &psi_ratios);
  GPU_XRAY_TRACE void evaluateDerivatives (MCWalkerConfiguration &W,
                            const opt_variables_type& optvars,
                            RealMatrix_t &dlogpsi,
                            RealMatrix_t &dlapl_over_psi);
  OneBodyJastrowOrbitalBspline(ParticleSet &centers, ParticleSet& elecs) :
    J1OrbitalSoA<FT>(centers,elecs),
    ElecRef(elecs),
    L("OneBodyJastrowOrbitalBspline::L"),
    Linv("OneBodyJastrowOrbitalBspline::Linv"),
    C("OneBodyJastrowOrbitalBspline::C"),
    UpdateListGPU("OneBodyJastrowOrbitalBspline::UpdateListGPU"),
    SumGPU("OneBodyJastrowOrbitalBspline::SumGPU"),
    GradLaplGPU("OneBodyJastrowOrbitalBspline::GradLaplGPU"),
    OneGradGPU("OneBodyJastrowOrbitalBspline::OneGradGPU"),
    SplineDerivsGPU("OneBodyJastrowOrbitalBspline::SplineDerivsGPU"),
    DerivListGPU("OneBodyJastrowOrbitalBspline::DerivListGPU"),
    NL_SplineCoefsListGPU("OneBodyJastrowOrbitalBspline::NL_SplineCoefsListGPU"),
    NL_JobListGPU("OneBodyJastrowOrbitalBspline::NL_JobListGPU"),
    NL_NumCoefsGPU("OneBodyJastrowOrbitalBspline::NL_NumCoefsGPU"),
    NL_NumQuadPointsGPU("OneBodyJastrowOrbitalBspline::NL_NumQuadPointsGPU"),
    NL_rMaxGPU("OneBodyJastrowOrbitalBspline::NL_rMaxGPU"),
    NL_QuadPointsGPU("OneBodyJastrowOrbitalBspline::NL_QuadPointsGPU"),
    NL_RatiosGPU("OneBodyJastrowOrbitalBspline::NL_RatiosGPU")
  {
    UsePBC = elecs.Lattice.SuperCellEnum;
    NumElecGroups = elecs.groups();
    SpeciesSet &sSet = centers.getSpeciesSet();
    NumCenterGroups = sSet.getTotalNum();
    //      NumCenterGroups = centers.groups();
    // std::cerr << "NumCenterGroups = " << NumCenterGroups << std::endl;
    GPUSplines.resize(NumCenterGroups,0);
    if (UsePBC)
    {
      gpu::host_vector<CudaReal> LHost(OHMMS_DIM*OHMMS_DIM),
          LinvHost(OHMMS_DIM*OHMMS_DIM);
      for (int i=0; i<OHMMS_DIM; i++)
        for (int j=0; j<OHMMS_DIM; j++)
        {
          LHost[OHMMS_DIM*i+j]    = (CudaReal)elecs.Lattice.a(i)[j];
          LinvHost[OHMMS_DIM*i+j] = (CudaReal)elecs.Lattice.b(j)[i];
        }
      L = LHost;
      Linv = LinvHost;
    }
    N = elecs.getTotalNum();
    // Copy center positions to GPU, sorting by GroupID
    gpu::host_vector<CudaReal> C_host(OHMMS_DIM*centers.getTotalNum());
    int index=0;
    for (int cgroup=0; cgroup<NumCenterGroups; cgroup++)
    {
      CenterFirst.push_back(index);
      for (int i=0; i<centers.getTotalNum(); i++)
      {
        if (centers.GroupID[i] == cgroup)
        {
          for (int dim=0; dim<OHMMS_DIM; dim++)
            C_host[OHMMS_DIM*index+dim] = centers.R[i][dim];
          index++;
        }
      }
      CenterLast.push_back(index-1);
    }
    // gpu::host_vector<CudaReal> C_host(OHMMS_DIM*centers.getTotalNum());
    // for (int i=0; i<centers.getTotalNum(); i++)
    // 	for (int dim=0; dim<OHMMS_DIM; dim++)
    // 	  C_host[OHMMS_DIM*i+dim] = centers.R[i][dim];
    C = C_host;
  }
};
}


#endif

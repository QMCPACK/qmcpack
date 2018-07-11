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
    
    
#ifndef TWO_BODY_JASTROW_ORBITAL_BSPLINE_AOS_H
#define TWO_BODY_JASTROW_ORBITAL_BSPLINE_AOS_H

#include "Particle/DistanceTableData.h"
#include "Particle/DistanceTable.h"
#include "QMCWaveFunctions/Jastrow/TwoBodyJastrowOrbital.h"
#include "QMCWaveFunctions/Jastrow/BsplineFunctor.h"
#include "Configuration.h"
#include "QMCWaveFunctions/Jastrow/CudaSpline.h"
#include "NLjobGPU.h"

namespace qmcplusplus
{

class TwoBodyJastrowOrbitalBsplineAoS :
  public TwoBodyJastrowOrbital<BsplineFunctor<OrbitalBase::RealType> >
{
private:
  bool UsePBC;
  typedef CUDA_PRECISION CudaReal;
  //typedef double CudaReal;

  std::vector<CudaSpline<CudaReal>*> GPUSplines, UniqueSplines;
  int MaxCoefs;
  ParticleSet &PtclRef;
  gpu::device_vector<CudaReal> L, Linv;

  gpu::device_vector<CudaReal*> UpdateListGPU;
  gpu::device_vector<CudaReal> SumGPU, GradLaplGPU, OneGradGPU;

  gpu::host_vector<CudaReal*> UpdateListHost;
  gpu::host_vector<CudaReal> SumHost, GradLaplHost, OneGradHost;
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
public:
  typedef BsplineFunctor<OrbitalBase::RealType> FT;
  typedef ParticleSet::Walker_t     Walker_t;

  GPU_XRAY_TRACE void freeGPUmem();
  GPU_XRAY_TRACE void checkInVariables(opt_variables_type& active);
  //void addFunc(const std::string& aname, int ia, int ib, FT* j);
  GPU_XRAY_TRACE void addFunc(int ia, int ib, FT* j);
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
  void ratio (std::vector<Walker_t*> &walkers,    std::vector<int> &iatList,
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

  GPU_XRAY_TRACE void resetParameters(const opt_variables_type& active);

  // Evaluates the derivatives of log psi and laplacian log psi w.r.t.
  // the parameters for optimization.  First index of the ValueMatrix is
  // the parameter.  The second is the walker.
  GPU_XRAY_TRACE void
  evaluateDerivatives (MCWalkerConfiguration &W,
                       const opt_variables_type& optvars,
                       RealMatrix_t &dlogpsi,
                       RealMatrix_t &dlapl_over_psi);

  //TwoBodyJastrowOrbitalBsplineAoS(ParticleSet& pset, bool is_master) :
  //  TwoBodyJastrowOrbital<BsplineFunctor<OrbitalBase::RealType> > (pset, is_master),
  TwoBodyJastrowOrbitalBsplineAoS(ParticleSet& pset, int tid) :
    TwoBodyJastrowOrbital<BsplineFunctor<OrbitalBase::RealType> > (pset, tid),
    PtclRef(pset),
    UpdateListGPU        ("TwoBodyJastrowOrbitalBsplineAoS::UpdateListGPU"),
    L                    ("TwoBodyJastrowOrbitalBsplineAoS::L"),
    Linv                 ("TwoBodyJastrowOrbitalBsplineAoS::Linv"),
    SumGPU               ("TwoBodyJastrowOrbitalBsplineAoS::SumGPU"),
    GradLaplGPU          ("TwoBodyJastrowOrbitalBsplineAoS::GradLaplGPU"),
    OneGradGPU           ("TwoBodyJastrowOrbitalBsplineAoS::OneGradGPU"),
    SplineDerivsGPU      ("TwoBodyJastrowOrbitalBsplineAoS::SplineDerivsGPU"),
    DerivListGPU         ("TwoBodyJastrowOrbitalBsplineAoS::DerivListGPU"),
    NL_SplineCoefsListGPU("TwoBodyJastrowOrbitalBsplineAoS::NL_SplineCoefsListGPU"),
    NL_JobListGPU        ("TwoBodyJastrowOrbitalBsplineAoS::NL_JobListGPU"),
    NL_NumCoefsGPU       ("TwoBodyJastrowOrbitalBsplineAoS::NL_NumCoefsGPU"),
    NL_NumQuadPointsGPU  ("TwoBodyJastrowOrbitalBsplineAoS::NL_NumQuadPointsGPU"),
    NL_rMaxGPU           ("TwoBodyJastrowOrbitalBsplineAoS::NL_rMaxGPU"),
    NL_QuadPointsGPU     ("TwoBodyJastrowOrbitalBsplineAoS::NL_QuadPointsGPU"),
    NL_RatiosGPU         ("TwoBodyJastrowOrbitalBsplineAoS::NL_RatiosGPU")
  {
    UsePBC = pset.Lattice.SuperCellEnum;
    app_log() << "UsePBC = " << UsePBC << std::endl;
    int nsp = NumGroups = pset.groups();
    GPUSplines.resize(nsp*nsp,0);
    if (UsePBC)
    {
      gpu::host_vector<CudaReal> LHost(OHMMS_DIM*OHMMS_DIM),
          LinvHost(OHMMS_DIM*OHMMS_DIM);
      for (int i=0; i<OHMMS_DIM; i++)
        for (int j=0; j<OHMMS_DIM; j++)
        {
          LHost[OHMMS_DIM*i+j]    = (CudaReal)pset.Lattice.a(i)[j];
          LinvHost[OHMMS_DIM*i+j] = (CudaReal)pset.Lattice.b(j)[i];
        }
      // for (int i=0; i<OHMMS_DIM; i++)
      // 	for (int j=0; j<OHMMS_DIM; j++) {
      // 	  double sum = 0.0;
      // 	  for (int k=0; k<OHMMS_DIM; k++)
      // 	    sum += LHost[OHMMS_DIM*i+k]*LinvHost[OHMMS_DIM*k+j];
      // 	  if (i == j) sum -= 1.0;
      // 	  if (std::abs(sum) > 1.0e-5) {
      // 	    app_error() << "sum = " << sum << std::endl;
      // 	    app_error() << "Linv * L != identity.\n";
      // 	    abort();
      // 	  }
      // 	}
      //       fprintf (stderr, "Identity should follow:\n");
      //       for (int i=0; i<3; i++){
      // 	for (int j=0; j<3; j++) {
      // 	  CudaReal val = 0.0f;
      // 	  for (int k=0; k<3; k++)
      // 	    val += LinvHost[3*i+k]*LHost[3*k+j];
      // 	  fprintf (stderr, "  %8.3f", val);
      // 	}
      // 	fprintf (stderr, "\n");
      //       }
      L = LHost;
      Linv = LinvHost;
    }
  }
};
}


#endif

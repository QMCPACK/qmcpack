#ifndef ONE_BODY_JASTROW_ORBITAL_BSPLINE_H
#define ONE_BODY_JASTROW_ORBITAL_BSPLINE_H

#include "Particle/DistanceTableData.h"
#include "Particle/DistanceTable.h"
#include "QMCWaveFunctions/Jastrow/OneBodyJastrowOrbital.h"
#include "QMCWaveFunctions/Jastrow/BsplineFunctor.h"
#include "QMCWaveFunctions/Jastrow/CudaSpline.h"
#include "QMCWaveFunctions/Jastrow/BsplineJastrowCuda.h"
#include "QMCWaveFunctions/Jastrow/BsplineJastrowCudaPBC.h"
#include "Configuration.h"

namespace qmcplusplus
{

class OneBodyJastrowOrbitalBspline :
  public OneBodyJastrowOrbital<BsplineFunctor<OrbitalBase::RealType> >
{
private:
  bool UsePBC;
  typedef CUDA_PRECISION CudaReal;
  //typedef double CudaReal;

  vector<CudaSpline<CudaReal>*> GPUSplines, UniqueSplines;
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
  vector<int> CenterFirst, CenterLast;
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
  typedef BsplineFunctor<OrbitalBase::RealType> FT;
  typedef ParticleSet::Walker_t     Walker_t;

  void resetParameters(const opt_variables_type& active);
  void checkInVariables(opt_variables_type& active);
  void addFunc(int ig, FT* j, int jg);
  void recompute(MCWalkerConfiguration &W, bool firstTime);
  void reserve (PointerPool<gpu::device_vector<CudaRealType> > &pool);
  void addLog (MCWalkerConfiguration &W, vector<RealType> &logPsi);
  void update (vector<Walker_t*> &walkers, int iat);
  void update (const vector<Walker_t*> &walkers, const vector<int> &iatList)
  {
    /* This function doesn't really need to return the ratio */
  }
  void ratio (MCWalkerConfiguration &W, int iat,
              vector<ValueType> &psi_ratios,	vector<GradType>  &grad,
              vector<ValueType> &lapl);
  void calcRatio (MCWalkerConfiguration &W, int iat,
                  vector<ValueType> &psi_ratios,	vector<GradType>  &grad,
                  vector<ValueType> &lapl);
  void addRatio (MCWalkerConfiguration &W, int iat,
                 vector<ValueType> &psi_ratios,	vector<GradType>  &grad,
                 vector<ValueType> &lapl);
  void ratio (vector<Walker_t*> &walkers,    vector<int> &iatList,
              vector<PosType> &rNew, vector<ValueType> &psi_ratios,
              vector<GradType>  &grad, vector<ValueType> &lapl)
  {
    /* This function doesn't really need to return the ratio */
  }


  void calcGradient(MCWalkerConfiguration &W, int iat,
                    vector<GradType> &grad);
  void addGradient(MCWalkerConfiguration &W, int iat,
                   vector<GradType> &grad);
  void gradLapl (MCWalkerConfiguration &W, GradMatrix_t &grads,
                 ValueMatrix_t &lapl);
  void NLratios (MCWalkerConfiguration &W,  vector<NLjob> &jobList,
                 vector<PosType> &quadPoints, vector<ValueType> &psi_ratios);
  void evaluateDerivatives (MCWalkerConfiguration &W,
                            const opt_variables_type& optvars,
                            ValueMatrix_t &dlogpsi,
                            ValueMatrix_t &dlapl_over_psi);
  OneBodyJastrowOrbitalBspline(ParticleSet &centers, ParticleSet& elecs) :
    OneBodyJastrowOrbital<BsplineFunctor<OrbitalBase::RealType> > (centers,elecs),
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
    // cerr << "NumCenterGroups = " << NumCenterGroups << endl;
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

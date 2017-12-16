//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//
// File created by: Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
    
    
#include "QMCHamiltonians/StressKinetic.h"
#include "Particle/DistanceTable.h"
#include "Particle/DistanceTableData.h"
#include "Utilities/ProgressReportEngine.h"
#include "QMCWaveFunctions/OrbitalBase.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include <numeric>

namespace qmcplusplus
{

StressKinetic::StressKinetic(ParticleSet& ref, TrialWaveFunction& Psi0) :
   ForceBase(ref,ref), Ps(ref), Psi(Psi0)
{
  ReportEngine PRE("StressKinetic","StressKinetic");
  //save source tag
  SourceID=ref.tag();
  //create a distance table: just to get the table name
 // DistanceTableData *d_aa = DistanceTable::add(ref);
 // PtclRefName=d_aa->Name;
  prefix="S_kin";
}

StressKinetic:: ~StressKinetic() { }


SymTensor<StressKinetic::RealType,OHMMS_DIM> StressKinetic::evaluateKineticSymTensor(ParticleSet& P)
{
  OrbitalBase::HessVector_t grad_grad_psi;
  Psi.evaluateHessian(P,grad_grad_psi);
  SymTensor<RealType,OHMMS_DIM> kinetic_tensor;
  SymTensor<ComplexType, OHMMS_DIM> complex_ktensor;

  
  for (int iat=0; iat<P.getTotalNum(); iat++)
  {
      complex_ktensor+= (grad_grad_psi[iat] + outerProduct(P.G[iat],P.G[iat]))*(-1.0/(2*P.Mass[iat]));
  }
  
  for (int i=0; i<OHMMS_DIM; i++)
	for(int j=i; j<OHMMS_DIM; j++)
	{
	  kinetic_tensor(i,j)=complex_ktensor(i,j).real();
	}
  return kinetic_tensor;
}


StressKinetic::Return_t
StressKinetic::evaluate(ParticleSet& P)
{
  stress=evaluateKineticSymTensor(P);
 // return Value;
 return 0.0;
}

QMCHamiltonianBase* StressKinetic::makeClone(ParticleSet& qp, TrialWaveFunction& psi)
{
  return new StressKinetic(*this);
}
}



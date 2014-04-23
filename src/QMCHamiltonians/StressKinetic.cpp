//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim and Kris Delaney
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
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

StressKinetic::Return_t
StressKinetic::evaluatePbyP(ParticleSet& P, int active)
{

    return 0.0;
}



QMCHamiltonianBase* StressKinetic::makeClone(ParticleSet& qp, TrialWaveFunction& psi)
{
  return new StressKinetic(*this);
}
}

/***************************************************************************
 * $RCSfile$   $Author: jtkrogel $
 * $Revision: 5976 $   $Date: 2013-09-13 13:39:44 -0500 (Fri, 13 Sep 2013) $
 * $Id: StressKinetic.cpp 5976 2013-09-13 18:39:44Z jtkrogel $
 ***************************************************************************/


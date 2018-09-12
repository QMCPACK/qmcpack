//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2018 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//   initially refactored from DiracDeterminantBase.cpp
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "QMCWaveFunctions/Fermion/DiracDeterminantSingle.h"
#include "Numerics/MatrixOperators.h"

namespace qmcplusplus
{
DiracDeterminantSingle::DiracDeterminantSingle(SPOSetPtr const &spos, int first):
  DiracDeterminantBase(first), Phi(spos)
{
  Optimizable=false;
  if(Phi->Optimizable)
    Optimizable=true;
  OrbitalName="DiracDeterminantSingle";
  registerTimers();
}

DiracDeterminantSingle::ValueType
DiracDeterminantSingle::ratioGrad(ParticleSet& P, int iat, GradType& grad_iat)
{
  SPOVGLTimer.start();
  Phi->evaluate(P, iat, psiV, dpsiV, d2psiV);
  SPOVGLTimer.stop();
  RatioTimer.start();
  WorkingIndex = iat-FirstIndex;
  UpdateMode=ORB_PBYP_PARTIAL;
  curRatio=simd::dot(psiM[WorkingIndex],psiV.data(),NumOrbitals);
  GradType rv=simd::dot(psiM[WorkingIndex],dpsiV.data(),NumOrbitals);
  grad_iat += ((RealType)1.0/curRatio) * rv;
  RatioTimer.stop();
  return curRatio;
}

void DiracDeterminantSingle::updateAfterSweep(ParticleSet& P,
					      ParticleSet::ParticleGradient_t& G,
					      ParticleSet::ParticleLaplacian_t& L)
{
  if(UpdateMode == ORB_PBYP_RATIO)
  { //need to compute dpsiM and d2psim. Use Phi->t_logpsi. Do not touch psiM!
    SPOVGLTimer.start();
    Phi->evaluate_notranspose(P,FirstIndex,LastIndex,psiM_temp,dpsiM,d2psiM);
    SPOVGLTimer.stop();
  }

  if(NumPtcls==1)
  {
    ValueType y = psiM(0,0);
    GradType rv = y*dpsiM(0,0);
    G[FirstIndex]+=rv;
    L[FirstIndex]+=y*d2psiM(0,0)-dot(rv,rv);
  }
  else
  {
    for(size_t i=0,iat=FirstIndex; i<NumPtcls; ++i,++iat)
    {
      mValueType dot_temp=simd::dot(psiM[i],d2psiM[i],NumOrbitals);
      mGradType rv=simd::dot(psiM[i],dpsiM[i],NumOrbitals);
      G[iat]+=rv;
      L[iat]+=dot_temp-dot(rv,rv);
    }
  }
}

DiracDeterminantSingle::ValueType
DiracDeterminantSingle::ratio(ParticleSet& P, int iat)
{
  UpdateMode=ORB_PBYP_RATIO;
  WorkingIndex = iat-FirstIndex;
  SPOVTimer.start();
  Phi->evaluate(P, iat, psiV);
  SPOVTimer.stop();
  RatioTimer.start();
  curRatio=simd::dot(psiM[WorkingIndex],psiV.data(),NumOrbitals);
  //curRatio = DetRatioByRow(psiM, psiV,WorkingIndex);
  RatioTimer.stop();
  return curRatio;
}

void DiracDeterminantSingle::evaluateRatiosAlltoOne(ParticleSet& P, std::vector<ValueType>& ratios)
{
  SPOVTimer.start();
  Phi->evaluate(P, -1, psiV);
  SPOVTimer.stop();
  MatrixOperators::product(psiM,psiV.data(),&ratios[FirstIndex]);
}

DiracDeterminantSingle::GradType DiracDeterminantSingle::evalGradSource(ParticleSet& P,
						ParticleSet& source,
						int iat)
{
  Phi->evaluateGradSource (P, FirstIndex, LastIndex, source, iat, grad_source_psiM);
  return simd::dot(psiM.data(),grad_source_psiM.data(),psiM.size());
}

DiracDeterminantSingle::GradType
DiracDeterminantSingle::evalGradSourcep
(ParticleSet& P, ParticleSet& source,int iat,
 TinyVector<ParticleSet::ParticleGradient_t, OHMMS_DIM> &grad_grad,
 TinyVector<ParticleSet::ParticleLaplacian_t,OHMMS_DIM> &lapl_grad)
{
  Phi->evaluateGradSource (P, FirstIndex, LastIndex, source, iat,
                           grad_source_psiM, grad_grad_source_psiM,
                           grad_lapl_source_psiM);
  Phi->evaluate_notranspose(P, FirstIndex, LastIndex, psiM_temp, dpsiM, d2psiM);

  invertPsiM(psiM_temp,psiM);

  GradMatrix_t &Phi_alpha(grad_source_psiM);
  GradMatrix_t &Grad_phi(dpsiM);
  ValueMatrix_t &Grad2_phi(d2psiM);
  HessMatrix_t &Grad_phi_alpha(grad_grad_source_psiM);
  GradMatrix_t &Grad2_phi_alpha(grad_lapl_source_psiM);
  GradType Psi_alpha_over_psi;
  Psi_alpha_over_psi = evalGradSource(P, source, iat);
  std::ofstream outfile;
  outfile.open("grad_psi_alpha_over_psi",std::ios::app);
  ValueMatrix_t toDet;
  ValueMatrix_t toDet_l;
  toDet.resize(2,2);
  toDet_l.resize(2,2);
  for (int ptcl=0; ptcl<NumPtcls; ptcl++)
  {
    ValueType Grad2_psi_over_psi(0.0);
    GradType Grad_psi_over_psi(0.0);
    HessType Grad_psi_alpha_over_psi(0.0);
    HessType one_row_change(0.0);
    HessType two_row_change(0.0);
    GradType one_row_change_l(0.0);
    GradType two_row_change_l(0.0);
    for (int el_dim=0; el_dim<OHMMS_DIM; el_dim++)
    {
      for (int orbital=0; orbital<NumOrbitals; orbital++)
      {
        Grad_psi_over_psi[el_dim]+=Grad_phi(ptcl,orbital)[el_dim]*psiM(ptcl,orbital);
        if (el_dim==0)
          Grad2_psi_over_psi+=Grad2_phi(ptcl,orbital)*psiM(ptcl,orbital);
      }
      for (int dim=0; dim<OHMMS_DIM; dim++)
      {
        one_row_change(dim,el_dim)=0.0;
        for (int orbital=0; orbital<NumOrbitals; orbital++)
        {
          one_row_change(dim,el_dim)+=Grad_phi_alpha(ptcl,orbital)(dim,el_dim)*psiM(ptcl,orbital);
          if (el_dim==0)
            one_row_change_l[dim]+=Grad2_phi_alpha(ptcl,orbital)[dim]*psiM(ptcl,orbital);
        }
        for (int ptcl2=0; ptcl2<NumPtcls; ptcl2++)
        {
          if (ptcl!=ptcl2)
          {
            toDet=0.0;
            toDet_l=0.0;
            for (int orbital=0; orbital<NumOrbitals; orbital++)
            {
              toDet(0,0)+=Grad_phi(ptcl,orbital)[el_dim]*psiM(ptcl,orbital);
              toDet_l(0,0)+=Grad2_phi(ptcl,orbital)*psiM(ptcl,orbital);
              toDet(0,1)+=Grad_phi(ptcl,orbital)[el_dim]*psiM(ptcl2,orbital);
              toDet_l(0,1)+=Grad2_phi(ptcl,orbital)*psiM(ptcl2,orbital);
              toDet(1,0)+=Phi_alpha(ptcl2,orbital)[dim]*psiM(ptcl,orbital);
              toDet_l(1,0)+=Phi_alpha(ptcl2,orbital)[dim]*psiM(ptcl,orbital);
              toDet(1,1)+=Phi_alpha(ptcl2,orbital)[dim]*psiM(ptcl2,orbital);
              toDet_l(1,1)+=Phi_alpha(ptcl2,orbital)[dim]*psiM(ptcl2,orbital);
            }
            two_row_change(dim,el_dim)+=toDet(0,0)*toDet(1,1)-toDet(1,0)*toDet(0,1);
            if (el_dim==0)
              two_row_change_l[dim]+=toDet_l(0,0)*toDet_l(1,1)-toDet_l(1,0)*toDet_l(0,1);
          }
        }
        Grad_psi_alpha_over_psi(dim,el_dim)=one_row_change(dim,el_dim)+two_row_change(dim,el_dim);
        outfile<<Grad_psi_alpha_over_psi(dim,el_dim)<< std::endl;
        grad_grad[dim][ptcl][el_dim]=one_row_change(dim,el_dim)+two_row_change(dim,el_dim)-
                                     Grad_psi_over_psi[el_dim]*Psi_alpha_over_psi[dim];
      }
    }
    for (int dim=0; dim<OHMMS_DIM; dim++)
    {
      lapl_grad[dim][ptcl]=0.0;
      lapl_grad[dim][ptcl]+=one_row_change_l[dim]+two_row_change_l[dim]- Psi_alpha_over_psi[dim]*Grad2_psi_over_psi;
      for (int el_dim=0; el_dim<OHMMS_DIM; el_dim++)
      {
        lapl_grad[dim][ptcl]-= (RealType)2.0*Grad_psi_alpha_over_psi(dim,el_dim)*Grad_psi_over_psi[el_dim];
        lapl_grad[dim][ptcl]+= (RealType)2.0*Psi_alpha_over_psi[dim]*(Grad_psi_over_psi[el_dim]*Grad_psi_over_psi[el_dim]);
      }
    }
  }
  outfile.close();
  return Psi_alpha_over_psi;
}


void DiracDeterminantSingle::evaluateRatios(VirtualParticleSet& VP, std::vector<ValueType>& ratios)
{
  const int nVP = VP.getTotalNum();
  const size_t memory_needed = nVP*NumOrbitals+Phi->estimateMemory(nVP);
  //std::cout << "debug " << memory_needed << " pool " << memoryPool.size() << std::endl;
  if(memoryPool.size()<memory_needed)
  {
    // usually in small systems
    for(int iat=0; iat<nVP; iat++)
    {
      SPOVTimer.start();
      Phi->evaluate(VP, iat, psiV);
      SPOVTimer.stop();
      RatioTimer.start();
      ratios[iat]=simd::dot(psiM[VP.refPtcl-FirstIndex],psiV.data(),NumOrbitals);
      RatioTimer.stop();
    }
  }
  else
  {
    const size_t offset = memory_needed-nVP*NumOrbitals;
    // SPO value result matrix. Always use existing memory
    Matrix<ValueType> psiT(memoryPool.data()+offset, nVP, NumOrbitals);
    // SPO scratch memory. Always use existing memory
    SPOSet::ValueAlignedVector_t SPOMem;
    SPOMem.attachReference((ValueType*)memoryPool.data(),offset);
    SPOVTimer.start();
    Phi->evaluateValues(VP, psiT, SPOMem, Phi->OrbitalSetSize);
    SPOVTimer.stop();
    RatioTimer.start();
    MatrixOperators::product(psiT, psiM[VP.refPtcl-FirstIndex], ratios.data());
    RatioTimer.stop();
  }
}

void DiracDeterminantSingle::evaluateHessian(ParticleSet& P, HessVector_t& grad_grad_psi)
{
  //IM A HACK.  Assumes evaluateLog has already been executed.
  Phi->evaluate_notranspose(P, FirstIndex, LastIndex, psiM_temp, dpsiM, grad_grad_source_psiM);
  invertPsiM(psiM_temp,psiM);

  phi_alpha_Minv = 0.0;
  grad_phi_Minv = 0.0;
  lapl_phi_Minv = 0.0;
  grad_phi_alpha_Minv = 0.0;
  //grad_grad_psi.resize(NumPtcls);

  for(int i=0, iat=FirstIndex; i<NumPtcls; i++, iat++)
  {
    GradType rv=simd::dot(psiM[i],dpsiM[i],NumOrbitals);
    //  HessType hess_tmp=simd::dot(psiM[i],grad_grad_source_psiM[i],NumOrbitals);
    HessType hess_tmp;
    hess_tmp=0.0;
    hess_tmp=simd::dot(psiM[i],grad_grad_source_psiM[i],NumOrbitals);
    grad_grad_psi[iat]=hess_tmp-outerProduct(rv,rv);
  }
}

DiracDeterminantSingle::GradType
DiracDeterminantSingle::evalGradSource
(ParticleSet& P, ParticleSet& source,int iat,
 TinyVector<ParticleSet::ParticleGradient_t, OHMMS_DIM> &grad_grad,
 TinyVector<ParticleSet::ParticleLaplacian_t,OHMMS_DIM> &lapl_grad)
{
  Phi->evaluateGradSource (P, FirstIndex, LastIndex, source, iat,
                           grad_source_psiM, grad_grad_source_psiM,
                           grad_lapl_source_psiM);
  // HACK HACK HACK
  // Phi->evaluate(P, FirstIndex, LastIndex, psiM, dpsiM, d2psiM);
  // psiM_temp = psiM;
  // LogValue=InvertWithLog(psiM.data(),NumPtcls,NumOrbitals,
  // 			   WorkSpace.data(),Pivot.data(),PhaseValue);
  // for (int i=0; i<NumPtcls; i++)
  //   for (int j=0; j<NumPtcls; j++) {
  // 	double val = 0.0;
  // 	for (int k=0; k<NumPtcls; k++)
  // 	  val += psiM(i,k) * psiM_temp(k,j);
  // 	val -= (i == j) ? 1.0 : 0.0;
  // 	if (std::abs(val) > 1.0e-12)
  // 	  std::cerr << "Error in inverse.\n";
  //   }
  // for (int i=0; i<NumPtcls; i++) {
  //   P.G[FirstIndex+i] = GradType();
  //   for (int j=0; j<NumOrbitals; j++)
  // 	P.G[FirstIndex+i] += psiM(i,j)*dpsiM(i,j);
  // }
  // Compute matrices
  phi_alpha_Minv = 0.0;
  grad_phi_Minv = 0.0;
  lapl_phi_Minv = 0.0;
  grad_phi_alpha_Minv = 0.0;
  for (int i=0; i<NumPtcls; i++)
    for (int j=0; j<NumOrbitals; j++)
    {
      lapl_phi_Minv(i,j) = 0.0;
      for (int k=0; k<NumOrbitals; k++)
        lapl_phi_Minv(i,j) += d2psiM(i,k)*psiM(j,k);
    }
  for (int dim=0; dim<OHMMS_DIM; dim++)
  {
    for (int i=0; i<NumPtcls; i++)
      for (int j=0; j<NumOrbitals; j++)
      {
        for (int k=0; k<NumOrbitals; k++)
        {
          phi_alpha_Minv(i,j)[dim] += grad_source_psiM(i,k)[dim] * psiM(j,k);
          grad_phi_Minv(i,j)[dim] += dpsiM(i,k)[dim] * psiM(j,k);
          for (int dim_el=0; dim_el<OHMMS_DIM; dim_el++)
            grad_phi_alpha_Minv(i,j)(dim, dim_el) +=
              grad_grad_source_psiM(i,k)(dim,dim_el)*psiM(j,k);
        }
      }
  }
  GradType gradPsi;
  for(int i=0, iel=FirstIndex; i<NumPtcls; i++, iel++)
  {
    HessType dval (0.0);
    GradType d2val(0.0);
    for (int dim=0; dim<OHMMS_DIM; dim++)
      for (int dim_el=0; dim_el<OHMMS_DIM; dim_el++)
        dval(dim,dim_el) = grad_phi_alpha_Minv(i,i)(dim,dim_el);
    for(int j=0; j<NumOrbitals; j++)
    {
      gradPsi += grad_source_psiM(i,j) * psiM(i,j);
      for (int dim=0; dim<OHMMS_DIM; dim++)
        for (int k=0; k<OHMMS_DIM; k++)
          dval(dim,k) -= phi_alpha_Minv(j,i)[dim]*grad_phi_Minv(i,j)[k];
    }
    for (int dim=0; dim<OHMMS_DIM; dim++)
    {
      for (int k=0; k<OHMMS_DIM; k++)
        grad_grad[dim][iel][k] += dval(dim,k);
      for (int j=0; j<NumOrbitals; j++)
      {
        // First term, eq 9
        lapl_grad[dim][iel] += grad_lapl_source_psiM(i,j)[dim] *
                               psiM(i,j);
        // Second term, eq 9
        if (j == i)
          for (int dim_el=0; dim_el<OHMMS_DIM; dim_el++)
            lapl_grad[dim][iel] -= (RealType)2.0 * grad_phi_alpha_Minv(j,i)(dim,dim_el)
                                   * grad_phi_Minv(i,j)[dim_el];
        // Third term, eq 9
        // First term, eq 10
        lapl_grad[dim][iel] -= phi_alpha_Minv(j,i)[dim]*lapl_phi_Minv(i,j);
        // Second term, eq 11
        for (int dim_el=0; dim_el<OHMMS_DIM; dim_el++)
          lapl_grad[dim][iel] += (RealType)2.0*phi_alpha_Minv(j,i)[dim] *
                                 grad_phi_Minv(i,i)[dim_el]*grad_phi_Minv(i,j)[dim_el];
      }
    }
  }
  return gradPsi;
}

void DiracDeterminantSingle::recompute(ParticleSet& P)
{
  SPOVGLTimer.start();
  Phi->evaluate_notranspose(P, FirstIndex, LastIndex, psiM_temp, dpsiM, d2psiM);
  SPOVGLTimer.stop();
  if(NumPtcls==1)
  {
    //CurrentDet=psiM(0,0);
    ValueType det=psiM_temp(0,0);
    psiM(0,0)=RealType(1)/det;
    LogValue = evaluateLogAndPhase(det,PhaseValue);
  }
  else
  {
    invertPsiM(psiM_temp,psiM);
  }
}

DiracDeterminantSingle* DiracDeterminantSingle::makeCopy(SPOSetPtr spo) const
{
  DiracDeterminantSingle* dclone= new DiracDeterminantSingle(spo);
  dclone->set(FirstIndex,LastIndex-FirstIndex);
  return dclone;
}


}

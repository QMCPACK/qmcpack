// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Bryan Clark, bclark@Princeton.edu, Princeton University
//                    Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    

#include "QMCWaveFunctions/Fermion/DiracDeterminantBase.h"
#include "Numerics/DeterminantOperators.h"
#include "Numerics/OhmmsBlas.h"
#include "Numerics/MatrixOperators.h"
#include "simd/simd.hpp"

namespace qmcplusplus
{

/** constructor
 *@param spos the single-particle orbital set
 *@param first index of the first particle
 */
DiracDeterminantBase::DiracDeterminantBase(SPOSetBasePtr const &spos, int first):
  NP(0), Phi(spos), FirstIndex(first)
  ,UpdateTimer("DiracDeterminantBase::update",timer_level_fine)
  ,RatioTimer("DiracDeterminantBase::ratio",timer_level_fine)
  ,InverseTimer("DiracDeterminantBase::inverse",timer_level_fine)
  ,BufferTimer("DiracDeterminantBase::buffer",timer_level_fine)
  ,SPOVTimer("DiracDeterminantBase::spoval",timer_level_fine)
  ,SPOVGLTimer("DiracDeterminantBase::spovgl",timer_level_fine)
{
  Optimizable=false;
  if(Phi->Optimizable)
    Optimizable=true;
  OrbitalName="DiracDeterminantBase";
  registerTimers();
}

///default destructor
DiracDeterminantBase::~DiracDeterminantBase() {}

#if 0
DiracDeterminantBase& DiracDeterminantBase::operator=(const DiracDeterminantBase& s)
{
  Bytes_in_WFBuffer=s.Bytes_in_WFBuffer;
  NP=0;
  resize(s.NumPtcls, s.NumOrbitals);
  return *this;
}
#endif

/** set the index of the first particle in the determinant and reset the size of the determinant
 *@param first index of first particle
 *@param nel number of particles in the determinant
 */
void DiracDeterminantBase::set(int first, int nel)
{
  FirstIndex = first;
  resize(nel,nel);
}

void DiracDeterminantBase::invertPsiM(const ValueMatrix_t& logdetT, ValueMatrix_t& invMat)
{
  InverseTimer.start();
#ifdef MIXED_PRECISION
  simd::transpose(logdetT.data(), NumOrbitals, logdetT.cols(), 
      psiM_hp.data(), NumOrbitals, psiM_hp.cols());
  ParticleSet::Scalar_t PhaseValue_hp;
  detEng_hp.invert(psiM_hp,true);
  LogValue = static_cast<RealType>(detEng_hp.LogDet);
  PhaseValue = static_cast<RealType>(detEng_hp.Phase);
  invMat = psiM_hp;
#else
  simd::transpose(logdetT.data(), NumOrbitals, logdetT.cols(), 
      invMat.data(), NumOrbitals, invMat.cols());
  detEng.invert(invMat,true);
  LogValue = detEng.LogDet;
  PhaseValue = detEng.Phase;
#endif
  InverseTimer.stop();
}



///reset the size: with the number of particles and number of orbtials
void DiracDeterminantBase::resize(int nel, int morb)
{
  int norb=morb;
  if(norb <= 0)
    norb = nel; // for morb == -1 (default)
  psiM.resize(nel,norb);
  psiM_temp.resize(nel,norb);
  dpsiM.resize(nel,norb);
  d2psiM.resize(nel,norb);
  psiV.resize(norb);
#ifdef MIXED_PRECISION
  psiM_hp.resize(nel,norb);
#endif
  LastIndex = FirstIndex + nel;
  NumPtcls=nel;
  NumOrbitals=norb;

  dpsiV.resize(NumOrbitals);
  d2psiV.resize(NumOrbitals);
  FirstAddressOfdV = &(dpsiM(0,0)[0]); //(*dpsiM.begin())[0]);
  LastAddressOfdV = FirstAddressOfdV + NumPtcls*NumOrbitals*DIM;
  // For forces
  /* Ye Luo, Apr 18th 2015
   * To save the memory used by every walker, the resizing the following giant matrices are commented.
   * When ZVZB forces and stresses are ready for deployment, R. Clay will take care of those matrices.
   */
  /*
  grad_source_psiM.resize(nel,norb);
  grad_lapl_source_psiM.resize(nel,norb);
  grad_grad_source_psiM.resize(nel,norb);
  phi_alpha_Minv.resize(nel,norb);
  grad_phi_Minv.resize(nel,norb);
  lapl_phi_Minv.resize(nel,norb);
  grad_phi_alpha_Minv.resize(nel,norb);
  */
}

DiracDeterminantBase::GradType
DiracDeterminantBase::evalGrad(ParticleSet& P, int iat)
{
  WorkingIndex = iat-FirstIndex;
  RatioTimer.start();
  DiracDeterminantBase::GradType g = simd::dot(psiM[WorkingIndex],dpsiM[WorkingIndex],NumOrbitals);
  RatioTimer.stop();
  return g;
}

DiracDeterminantBase::ValueType
DiracDeterminantBase::ratioGrad(ParticleSet& P, int iat, GradType& grad_iat)
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

/** move was accepted, update the real container
*/
void DiracDeterminantBase::acceptMove(ParticleSet& P, int iat)
{
  PhaseValue += evaluatePhase(curRatio);
  LogValue +=std::log(std::abs(curRatio));
  UpdateTimer.start();
  detEng.updateRow(psiM,psiV.data(),WorkingIndex,curRatio);
  if(UpdateMode == ORB_PBYP_PARTIAL)
  {
    simd::copy(dpsiM[WorkingIndex],  dpsiV.data(),  NumOrbitals);
    simd::copy(d2psiM[WorkingIndex], d2psiV.data(), NumOrbitals);
  }
  UpdateTimer.stop();
  curRatio=1.0;
}

/** move was rejected. copy the real container to the temporary to move on
*/
void DiracDeterminantBase::restore(int iat)
{
  curRatio=1.0;
}

void DiracDeterminantBase::updateAfterSweep(ParticleSet& P,
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

void
DiracDeterminantBase::registerData(ParticleSet& P, WFBufferType& buf)
{
  // Ye: no idea about NP.
  if(NP == 0) //first time, allocate once
  {
    NP=P.getTotalNum();
  }

  if ( Bytes_in_WFBuffer == 0 )
  {
    //add the data: inverse, gradient and laplacian
    Bytes_in_WFBuffer = buf.current();
    buf.add(psiM.first_address(),psiM.last_address());
    buf.add(FirstAddressOfdV,LastAddressOfdV);
    buf.add(d2psiM.first_address(),d2psiM.last_address());
    Bytes_in_WFBuffer = buf.current()-Bytes_in_WFBuffer;
    // free local space
    psiM.free();
    dpsiM.free();
    d2psiM.free();
  }
  else
  {
    buf.forward(Bytes_in_WFBuffer);
  }
  buf.add(LogValue);
  buf.add(PhaseValue);
}

DiracDeterminantBase::RealType DiracDeterminantBase::updateBuffer(ParticleSet& P,
    WFBufferType& buf, bool fromscratch)
{
  if(fromscratch)
  {
    LogValue=evaluateLog(P,P.G,P.L);
  }
  else
  {
    updateAfterSweep(P,P.G,P.L);
  }
  BufferTimer.start();
  buf.forward(Bytes_in_WFBuffer);
  buf.put(LogValue);
  buf.put(PhaseValue);
  BufferTimer.stop();
  return LogValue;
}

void DiracDeterminantBase::copyFromBuffer(ParticleSet& P, WFBufferType& buf)
{
  BufferTimer.start();
  psiM.attachReference(buf.lendReference<ValueType>(psiM.size()));
  dpsiM.attachReference(buf.lendReference<GradType>(dpsiM.size()));
  d2psiM.attachReference(buf.lendReference<ValueType>(d2psiM.size()));
  buf.get(LogValue);
  buf.get(PhaseValue);
  BufferTimer.stop();
}

/** return the ratio only for the  iat-th partcle move
 * @param P current configuration
 * @param iat the particle thas is being moved
 */
DiracDeterminantBase::ValueType DiracDeterminantBase::ratio(ParticleSet& P, int iat)
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

void DiracDeterminantBase::evaluateRatios(VirtualParticleSet& VP, std::vector<ValueType>& ratios)
{
  Matrix<ValueType> psiT(ratios.size(),NumOrbitals);
  Phi->evaluateValues(VP,psiT);
  MatrixOperators::product(psiT,psiM[VP.activePtcl-FirstIndex],&ratios[0]);
}

void DiracDeterminantBase::evaluateRatiosAlltoOne(ParticleSet& P, std::vector<ValueType>& ratios)
{
  SPOVTimer.start();
  Phi->evaluate(P, -1, psiV);
  SPOVTimer.stop();
  MatrixOperators::product(psiM,psiV.data(),&ratios[FirstIndex]);
}

DiracDeterminantBase::GradType
DiracDeterminantBase::evalGradSource(ParticleSet& P, ParticleSet& source,
                                     int iat)
{
  Phi->evaluateGradSource (P, FirstIndex, LastIndex, source, iat, grad_source_psiM);
  return simd::dot(psiM.data(),grad_source_psiM.data(),psiM.size());
}

DiracDeterminantBase::GradType
DiracDeterminantBase::evalGradSourcep
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

void DiracDeterminantBase::evaluateHessian(ParticleSet& P, HessVector_t& grad_grad_psi)
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

DiracDeterminantBase::GradType
DiracDeterminantBase::evalGradSource
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


/** Calculate the value of the Dirac determinant for particles
 *@param P input configuration containing N particles
 *@param G a vector containing N gradients
 *@param L a vector containing N laplacians
 *@return the value of the determinant
 *
 *\f$ (first,first+nel). \f$  Add the gradient and laplacian
 *contribution of the determinant to G(radient) and L(aplacian)
 *for local energy calculations.
 */
DiracDeterminantBase::ValueType
DiracDeterminantBase::evaluate(ParticleSet& P,
                               ParticleSet::ParticleGradient_t& G,
                               ParticleSet::ParticleLaplacian_t& L)
{
  RealType logval = evaluateLog(P, G, L);
#if defined(QMC_COMPLEX)
  RealType ratioMag = std::exp(logval);
  return std::complex<OHMMS_PRECISION>(std::cos(PhaseValue)*ratioMag,std::sin(PhaseValue)*ratioMag);
#else
  return std::cos(PhaseValue)*std::exp(logval);
#endif
}


DiracDeterminantBase::RealType
DiracDeterminantBase::evaluateLog(ParticleSet& P,
                                  ParticleSet::ParticleGradient_t& G,
                                  ParticleSet::ParticleLaplacian_t& L)
{
  recompute(P);

  if(NumPtcls==1)
  {
    ValueType y=psiM(0,0);
    GradType rv = y*dpsiM(0,0);
    G[FirstIndex] += rv;
    L[FirstIndex] += y*d2psiM(0,0) - dot(rv,rv);
  }
  else
  {
    for(int i=0, iat=FirstIndex; i<NumPtcls; i++, iat++)
    {
      mGradType rv=simd::dot(psiM[i],dpsiM[i],NumOrbitals);
      mValueType lap=simd::dot(psiM[i],d2psiM[i],NumOrbitals);
      G[iat] += rv;
      L[iat] += lap - dot(rv,rv);
    }
  }
  return LogValue;
}

void
DiracDeterminantBase::recompute(ParticleSet& P)
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

void
DiracDeterminantBase::evaluateDerivatives(ParticleSet& P,
    const opt_variables_type& active,
    std::vector<RealType>& dlogpsi,
    std::vector<RealType>& dhpsioverpsi)
{
}

OrbitalBasePtr DiracDeterminantBase::makeClone(ParticleSet& tqp) const
{
  APP_ABORT(" Illegal action. Cannot use DiracDeterminantBase::makeClone");
  return 0;
}

DiracDeterminantBase* DiracDeterminantBase::makeCopy(SPOSetBasePtr spo) const
{
  DiracDeterminantBase* dclone= new DiracDeterminantBase(spo);
  dclone->set(FirstIndex,LastIndex-FirstIndex);
  return dclone;
}

DiracDeterminantBase::DiracDeterminantBase(const DiracDeterminantBase& s)
  : OrbitalBase(s), NP(0), Phi(s.Phi), FirstIndex(s.FirstIndex)
  ,UpdateTimer(s.UpdateTimer)
  ,RatioTimer(s.RatioTimer)
  ,InverseTimer(s.InverseTimer)
  ,BufferTimer(s.BufferTimer)
  ,SPOVTimer(s.SPOVTimer)
  ,SPOVGLTimer(s.SPOVGLTimer)
{
  registerTimers();
  this->resize(s.NumPtcls,s.NumOrbitals);
}

//SPOSetBasePtr  DiracDeterminantBase::clonePhi() const
//{
//  return Phi->makeClone();
//}

void DiracDeterminantBase::registerTimers()
{
  UpdateTimer.reset();
  RatioTimer.reset();
  TimerManager.addTimer (&UpdateTimer);
  TimerManager.addTimer (&RatioTimer);
  TimerManager.addTimer (&InverseTimer);
  TimerManager.addTimer (&BufferTimer);
  TimerManager.addTimer (&SPOVTimer);
  TimerManager.addTimer (&SPOVGLTimer);
}

}

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
    
    

#include "QMCWaveFunctions/Fermion/DiracDeterminant.h"
#include "Numerics/DeterminantOperators.h"
#include "Numerics/OhmmsBlas.h"
#include "Numerics/BlasThreadingEnv.h"
#include "Numerics/MatrixOperators.h"
#include "simd/simd.hpp"

namespace qmcplusplus
{

/** constructor
 *@param spos the single-particle orbital set
 *@param first index of the first particle
 */
DiracDeterminant::DiracDeterminant(SPOSetPtr const spos, int first):
  DiracDeterminantBase(spos,first), ndelay(1)
{
  ClassName = "DiracDeterminant";
}

///default destructor
DiracDeterminant::~DiracDeterminant() {}

/** set the index of the first particle in the determinant and reset the size of the determinant
 *@param first index of first particle
 *@param nel number of particles in the determinant
 */
void DiracDeterminant::set(int first, int nel, int delay)
{
  FirstIndex = first;
  ndelay = delay;
  resize(nel,nel);
}

void DiracDeterminant::invertPsiM(const ValueMatrix_t& logdetT, ValueMatrix_t& invMat)
{
  InverseTimer.start();
  {
    BlasThreadingEnv knob(getNumThreadsNested());
#ifdef MIXED_PRECISION
    simd::transpose(logdetT.data(), NumOrbitals, logdetT.cols(),
                    psiM_hp.data(), NumOrbitals, psiM_hp.cols());
    detEng.invert(psiM_hp,true);
    LogValue = static_cast<RealType>(detEng.LogDet);
    PhaseValue = static_cast<RealType>(detEng.Phase);
    invMat = psiM_hp;
#else
    simd::transpose(logdetT.data(), NumOrbitals, logdetT.cols(),
                    invMat.data(), NumOrbitals, invMat.cols());
    detEng.invert(invMat,true);
    LogValue = detEng.LogDet;
    PhaseValue = detEng.Phase;
#endif
  } // end of BlasThreadingEnv
  InverseTimer.stop();
}



///reset the size: with the number of particles and number of orbtials
void DiracDeterminant::resize(int nel, int morb)
{
  int norb=morb;
  if(norb <= 0)
    norb = nel; // for morb == -1 (default)
  updateEng.resize(norb,ndelay);
  psiM.resize(nel,norb);
  dpsiM.resize(nel,norb);
  d2psiM.resize(nel,norb);
  psiV.resize(norb);
  psiM_temp.resize(nel,norb);
  if( typeid(ValueType) != typeid(mValueType) )
    psiM_hp.resize(nel,norb);
  LastIndex = FirstIndex + nel;
  NumPtcls=nel;
  NumOrbitals=norb;

  dpsiV.resize(NumOrbitals);
  d2psiV.resize(NumOrbitals);
  FirstAddressOfdV = &(dpsiM(0,0)[0]); //(*dpsiM.begin())[0]);
  LastAddressOfdV = FirstAddressOfdV + NumPtcls*NumOrbitals*DIM;
  
  if(ionDerivs)
  {
    grad_source_psiM.resize(nel,norb);
    grad_lapl_source_psiM.resize(nel,norb);
    grad_grad_source_psiM.resize(nel,norb);
    phi_alpha_Minv.resize(nel,norb);
    grad_phi_Minv.resize(nel,norb);
    lapl_phi_Minv.resize(nel,norb);
    grad_phi_alpha_Minv.resize(nel,norb);
  }
}

DiracDeterminant::GradType
DiracDeterminant::evalGrad(ParticleSet& P, int iat)
{
  const int WorkingIndex = iat-FirstIndex;
  RatioTimer.start();
  GradType g = updateEng.evalGrad(psiM, WorkingIndex, dpsiM[WorkingIndex]);
  RatioTimer.stop();
  return g;
}

DiracDeterminant::ValueType
DiracDeterminant::ratioGrad(ParticleSet& P, int iat, GradType& grad_iat)
{
  SPOVGLTimer.start();
  Phi->evaluate(P, iat, psiV, dpsiV, d2psiV);
  SPOVGLTimer.stop();
  RatioTimer.start();
  const int WorkingIndex = iat-FirstIndex;
  UpdateMode=ORB_PBYP_PARTIAL;
  GradType rv;
  curRatio = updateEng.ratioGrad(psiM, WorkingIndex, psiV, dpsiV, rv);
  grad_iat += ((RealType)1.0/curRatio) * rv;
  RatioTimer.stop();
  return curRatio;
}

/** move was accepted, update the real container
*/
void DiracDeterminant::acceptMove(ParticleSet& P, int iat)
{
  const int WorkingIndex = iat-FirstIndex;
  PhaseValue += evaluatePhase(curRatio);
  LogValue +=std::log(std::abs(curRatio));
  UpdateTimer.start();
  updateEng.acceptRow(psiM,WorkingIndex,psiV);
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
void DiracDeterminant::restore(int iat)
{
  curRatio=1.0;
}

void DiracDeterminant::completeUpdates()
{
  UpdateTimer.start();
  updateEng.updateInvMat(psiM);
  UpdateTimer.stop();
}

void DiracDeterminant::updateAfterSweep(ParticleSet& P,
      ParticleSet::ParticleGradient_t& G,
      ParticleSet::ParticleLaplacian_t& L)
{
  if(UpdateMode == ORB_PBYP_RATIO)
  { //need to compute dpsiM and d2psiM. Do not touch psiM!
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
DiracDeterminant::registerData(ParticleSet& P, WFBufferType& buf)
{
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

DiracDeterminant::RealType DiracDeterminant::updateBuffer(ParticleSet& P,
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

void DiracDeterminant::copyFromBuffer(ParticleSet& P, WFBufferType& buf)
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
DiracDeterminant::ValueType DiracDeterminant::ratio(ParticleSet& P, int iat)
{
  UpdateMode=ORB_PBYP_RATIO;
  const int WorkingIndex = iat-FirstIndex;
  SPOVTimer.start();
  Phi->evaluate(P, iat, psiV);
  SPOVTimer.stop();
  RatioTimer.start();
  curRatio = updateEng.ratio(psiM, WorkingIndex, psiV);
  RatioTimer.stop();
  return curRatio;
}

void DiracDeterminant::evaluateRatios(VirtualParticleSet& VP, std::vector<ValueType>& ratios)
{
  ValueVector_t psiM_row(psiM[VP.refPtcl-FirstIndex], psiM.cols());
  SPOVTimer.start();
  Phi->evaluateDetRatios(VP, psiV, psiM_row, ratios);
  SPOVTimer.stop();
}

void DiracDeterminant::evaluateRatiosAlltoOne(ParticleSet& P, std::vector<ValueType>& ratios)
{
  SPOVTimer.start();
  Phi->evaluate(P, -1, psiV);
  SPOVTimer.stop();
  MatrixOperators::product(psiM,psiV.data(),&ratios[FirstIndex]);
}

DiracDeterminant::GradType
DiracDeterminant::evalGradSource(ParticleSet& P, ParticleSet& source,
                                     int iat)
{
  if(!ionDerivs) APP_ABORT("DiracDeterminant::evalGradSource.  Determinant not initialized for force computations.");
  Phi->evaluateGradSource (P, FirstIndex, LastIndex, source, iat, grad_source_psiM);
  return simd::dot(psiM.data(),grad_source_psiM.data(),psiM.size());
}

DiracDeterminant::GradType
DiracDeterminant::evalGradSourcep
(ParticleSet& P, ParticleSet& source,int iat,
 TinyVector<ParticleSet::ParticleGradient_t, OHMMS_DIM> &grad_grad,
 TinyVector<ParticleSet::ParticleLaplacian_t,OHMMS_DIM> &lapl_grad)
{
  if(!ionDerivs) APP_ABORT("DiracDeterminant::evalGradSourcep.  Determinant not initialized for force computations.");
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

void DiracDeterminant::evaluateHessian(ParticleSet& P, HessVector_t& grad_grad_psi)
{
  // Hessian is not often used, so only resize/allocate if used
  grad_grad_source_psiM.resize(psiM.rows(),psiM.cols());
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

DiracDeterminant::GradType
DiracDeterminant::evalGradSource
(ParticleSet& P, ParticleSet& source,int iat,
 TinyVector<ParticleSet::ParticleGradient_t, OHMMS_DIM> &grad_grad,
 TinyVector<ParticleSet::ParticleLaplacian_t,OHMMS_DIM> &lapl_grad)
{
  if(!ionDerivs) APP_ABORT("DiracDeterminant::evalGradSource.  Determinant not initialized for force computations.");
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


/** Calculate the log value of the Dirac determinant for particles
 *@param P input configuration containing N particles
 *@param G a vector containing N gradients
 *@param L a vector containing N laplacians
 *@return the value of the determinant
 *
 *\f$ (first,first+nel). \f$  Add the gradient and laplacian
 *contribution of the determinant to G(radient) and L(aplacian)
 *for local energy calculations.
 */
DiracDeterminant::RealType
DiracDeterminant::evaluateLog(ParticleSet& P,
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
DiracDeterminant::recompute(ParticleSet& P)
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
DiracDeterminant::evaluateDerivatives(ParticleSet& P,
    const opt_variables_type& active,
    std::vector<RealType>& dlogpsi,
    std::vector<RealType>& dhpsioverpsi)
{
}

DiracDeterminant* DiracDeterminant::makeCopy(SPOSetPtr spo) const
{
  DiracDeterminant* dclone= new DiracDeterminant(spo);
  dclone->set(FirstIndex,LastIndex-FirstIndex,ndelay);
  return dclone;
}

}

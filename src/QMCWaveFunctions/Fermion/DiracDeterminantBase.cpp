//////////////////////////////////////////////////////////////////
// (c) Copyright 1998-2002,2003- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
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
  ,UpdateTimer("DiracDeterminantBase::update")
  ,RatioTimer("DiracDeterminantBase::ratio")
  ,InverseTimer("DiracDeterminantBase::inverse")
  ,BufferTimer("DiracDeterminantBase::buffer")
  ,SPOVTimer("DiracDeterminantBase::spoval")
  ,SPOVGLTimer("DiracDeterminantBase::spovgl")
{
  Optimizable=false;
  if(Phi->Optimizable)
    Optimizable=true;
  OrbitalName="DiracDeterminantBase";
  registerTimers();
}

///default destructor
DiracDeterminantBase::~DiracDeterminantBase() {}

DiracDeterminantBase& DiracDeterminantBase::operator=(const DiracDeterminantBase& s)
{
  NP=0;
  resize(s.NumPtcls, s.NumOrbitals);
  return *this;
}

/** set the index of the first particle in the determinant and reset the size of the determinant
 *@param first index of first particle
 *@param nel number of particles in the determinant
 */
void DiracDeterminantBase::set(int first, int nel)
{
  FirstIndex = first;
  resize(nel,nel);
}


///reset the size: with the number of particles and number of orbtials
void DiracDeterminantBase::resize(int nel, int morb)
{
  int norb=morb;
  if(norb <= 0)
    norb = nel; // for morb == -1 (default)
  psiM.resize(nel,norb);
  dpsiM.resize(nel,norb);
  d2psiM.resize(nel,norb);
  psiM_temp.resize(nel,norb);
  dpsiM_temp.resize(nel,norb);
  d2psiM_temp.resize(nel,norb);
  psiMinv.resize(nel,norb);
  psiV.resize(norb);
  WorkSpace.resize(nel);
  Pivot.resize(nel);
  LastIndex = FirstIndex + nel;
  NumPtcls=nel;
  NumOrbitals=norb;
  // For forces
  grad_source_psiM.resize(nel,norb);
  grad_lapl_source_psiM.resize(nel,norb);
  grad_grad_source_psiM.resize(nel,norb);
  phi_alpha_Minv.resize(nel,norb);
  grad_phi_Minv.resize(nel,norb);
  lapl_phi_Minv.resize(nel,norb);
  grad_phi_alpha_Minv.resize(nel,norb);
}

DiracDeterminantBase::RealType
DiracDeterminantBase::registerData(ParticleSet& P, PooledData<RealType>& buf)
{
  if(NP == 0)
    //first time, allocate once
  {
    //int norb = cols();
    dpsiV.resize(NumOrbitals);
    d2psiV.resize(NumOrbitals);
    workV1.resize(NumOrbitals);
    workV2.resize(NumOrbitals);
    NP=P.getTotalNum();
    myG.resize(NP);
    myL.resize(NP);
    myG_temp.resize(NP);
    myL_temp.resize(NP);
    FirstAddressOfG = &myG[0][0];
    LastAddressOfG = FirstAddressOfG + NP*DIM;
    FirstAddressOfdV = &(dpsiM(0,0)[0]); //(*dpsiM.begin())[0]);
    LastAddressOfdV = FirstAddressOfdV + NumPtcls*NumOrbitals*DIM;
  }
  myG=0.0;
  myL=0.0;
  //ValueType x=evaluate(P,myG,myL);
  LogValue=evaluateLog(P,myG,myL);
  P.G += myG;
  P.L += myL;
  //add the data: determinant, inverse, gradient and laplacians
  buf.add(psiM.first_address(),psiM.last_address());
  buf.add(FirstAddressOfdV,LastAddressOfdV);
  buf.add(d2psiM.first_address(),d2psiM.last_address());
  buf.add(myL.first_address(), myL.last_address());
  buf.add(FirstAddressOfG,LastAddressOfG);
  buf.add(LogValue);
  buf.add(PhaseValue);
  return LogValue;
}

DiracDeterminantBase::RealType DiracDeterminantBase::updateBuffer(ParticleSet& P,
    PooledData<RealType>& buf, bool fromscratch)
{
  //myG=0.0;
  //myL=0.0;
  if(fromscratch)
  {
    LogValue=evaluateLog(P,myG,myL);
    UpdateTimer.start();
  }
  else
  {
    if(UpdateMode == ORB_PBYP_RATIO)
    {
      SPOVGLTimer.start();
      Phi->evaluate(P, FirstIndex, LastIndex, psiM_temp,dpsiM, d2psiM);
      SPOVGLTimer.stop();
    }
    UpdateTimer.start();
    if(NumPtcls==1)
    {
      // ValueType y=1.0/psiM_temp(0,0);
      // psiM(0,0)=y;
      // GradType rv = y*dpsiM(0,0);
      // myG(FirstIndex) += rv;
      // myL(FirstIndex) += y*d2psiM(0,0) - dot(rv,rv);
      ValueType y = psiM(0,0);
      GradType rv = y*dpsiM(0,0);
      P.G[FirstIndex]+=(myG[FirstIndex]=rv);
      P.L[FirstIndex]+=(myL[FirstIndex]=y*d2psiM(0,0)-dot(rv,rv));
      //myG(FirstIndex) += rv;
      //myL(FirstIndex) += y*d2psiM(0,0) - dot(rv,rv);
      //P.G += myG;
      //P.L += myL;
    }
    else
    {
      for(int i=0,iat=FirstIndex; i<NumPtcls; ++i,++iat)
      {
        myG[iat]=simd::dot(psiM[i],dpsiM[i],NumOrbitals);
        myL[iat]=simd::dot(psiM[i],d2psiM[i],NumOrbitals)-dot(myG[iat],myG[iat]);
      }
      for(int iat=FirstIndex; iat<LastIndex; ++iat)
        P.G[iat] += myG[iat];
      for(int iat=FirstIndex; iat<LastIndex; ++iat)
        P.L[iat] += myL[iat];
    }
  }
  UpdateTimer.stop();
  BufferTimer.start();
  //copy psiM to psiM_temp
  //psiM_temp=psiM;
  simd::copy(psiM_temp.data(),psiM.data(),psiM.size());
  buf.put(psiM.first_address(),psiM.last_address());
  buf.put(FirstAddressOfdV,LastAddressOfdV);
  buf.put(d2psiM.first_address(),d2psiM.last_address());
  buf.put(myL.first_address(), myL.last_address());
  buf.put(FirstAddressOfG,LastAddressOfG);
  buf.put(LogValue);
  buf.put(PhaseValue);
  BufferTimer.stop();
  return LogValue;
}

void DiracDeterminantBase::copyFromBuffer(ParticleSet& P, PooledData<RealType>& buf)
{
  BufferTimer.start();
  buf.get(psiM.first_address(),psiM.last_address());
  buf.get(FirstAddressOfdV,LastAddressOfdV);
  buf.get(d2psiM.first_address(),d2psiM.last_address());
  buf.get(myL.first_address(), myL.last_address());
  buf.get(FirstAddressOfG,LastAddressOfG);
  buf.get(LogValue);
  buf.get(PhaseValue);
  //re-evaluate it for testing
  //Phi.evaluate(P, FirstIndex, LastIndex, psiM, dpsiM, d2psiM);
  //CurrentDet = Invert(psiM.data(),NumPtcls,NumOrbitals);
  //need extra copy for gradient/laplacian calculations without updating it
  //psiM_temp = psiM;
  //dpsiM_temp = dpsiM;
  //d2psiM_temp = d2psiM;
  simd::copy(psiM_temp.data(),  psiM.data(),  psiM.size());
  simd::copy(dpsiM_temp.data(), dpsiM.data(), dpsiM.size());
  simd::copy(d2psiM_temp.data(),d2psiM.data(),d2psiM.size());
  BufferTimer.stop();
}

/** dump the inverse to the buffer
*/
void DiracDeterminantBase::dumpToBuffer(ParticleSet& P, PooledData<RealType>& buf)
{
  APP_ABORT("DiracDeterminantBase::dumpToBuffer");
  buf.add(psiM.first_address(),psiM.last_address());
}

/** copy the inverse from the buffer
*/
void DiracDeterminantBase::dumpFromBuffer(ParticleSet& P, PooledData<RealType>& buf)
{
  APP_ABORT("DiracDeterminantBase::dumpFromBuffer");
  buf.get(psiM.first_address(),psiM.last_address());
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
  curRatio = DetRatioByRow(psiM, psiV,WorkingIndex);
  RatioTimer.stop();
  return curRatio;
}

void DiracDeterminantBase::evaluateRatios(VirtualParticleSet& VP, vector<ValueType>& ratios)
{
  Matrix<ValueType> psiT(ratios.size(),NumOrbitals);
  Phi->evaluateValues(VP,psiT);
  MatrixOperators::product(psiT,psiM[VP.activePtcl-FirstIndex],&ratios[0]);
}

void DiracDeterminantBase::get_ratios(ParticleSet& P, vector<ValueType>& ratios)
{
  SPOVTimer.start();
  Phi->evaluate(P, 0, psiV);
  SPOVTimer.stop();
  MatrixOperators::product(psiM,psiV.data(),&ratios[FirstIndex]);
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

DiracDeterminantBase::GradType
DiracDeterminantBase::evalGradSource(ParticleSet& P, ParticleSet& source,
                                     int iat)
{
  Phi->evaluateGradSource (P, FirstIndex, LastIndex, source, iat, grad_source_psiM);
//     Phi->evaluate(P, FirstIndex, LastIndex, psiM, dpsiM, d2psiM);
//     LogValue=InvertWithLog(psiM.data(),NumPtcls,NumOrbitals,WorkSpace.data(),Pivot.data(),PhaseValue);
  return simd::dot(psiM.data(),grad_source_psiM.data(),psiM.size());
  //const ValueType* restrict yptr=psiM[0];
  //const GradType* restrict dyptr=grad_source_psiM[0];
  //for (int i=0; i<NumPtcls; i++)
  //  for(int j=0; j<NumOrbitals; j++)
  //    //rv += (*yptr++) *(*dyptr++);
  //    rv += grad_source_psiM(i,j) * psiM(i,j);
  //// HACK HACK
  ////return (grad_source_psiM(1,3));
  //return rv;
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
  Phi->evaluate(P, FirstIndex, LastIndex, psiM, dpsiM, d2psiM);
  InverseTimer.start();
  LogValue=InvertWithLog(psiM.data(),NumPtcls,NumOrbitals,
                         WorkSpace.data(),Pivot.data(),PhaseValue);
  InverseTimer.stop();
  GradMatrix_t &Phi_alpha(grad_source_psiM);
  GradMatrix_t &Grad_phi(dpsiM);
  ValueMatrix_t &Grad2_phi(d2psiM);
  HessMatrix_t &Grad_phi_alpha(grad_grad_source_psiM);
  GradMatrix_t &Grad2_phi_alpha(grad_lapl_source_psiM);
  //    vector<ValueType> grad_psi_over_psi_vector;
  //    vector<ValueType> grad_psi_alpha_over_psi_vector;
  GradType Psi_alpha_over_psi;
  Psi_alpha_over_psi = evalGradSource(P, source, iat);
  ofstream outfile;
  outfile.open("grad_psi_alpha_over_psi",ios::app);
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
        Grad_psi_over_psi(el_dim)+=Grad_phi(ptcl,orbital)(el_dim)*psiM(ptcl,orbital);
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
            one_row_change_l(dim)+=Grad2_phi_alpha(ptcl,orbital)(dim)*psiM(ptcl,orbital);
        }
        for (int ptcl2=0; ptcl2<NumPtcls; ptcl2++)
        {
          if (ptcl!=ptcl2)
          {
            toDet=0.0;
            toDet_l=0.0;
            for (int orbital=0; orbital<NumOrbitals; orbital++)
            {
              toDet(0,0)+=Grad_phi(ptcl,orbital)(el_dim)*psiM(ptcl,orbital);
              toDet_l(0,0)+=Grad2_phi(ptcl,orbital)*psiM(ptcl,orbital);
              toDet(0,1)+=Grad_phi(ptcl,orbital)(el_dim)*psiM(ptcl2,orbital);
              toDet_l(0,1)+=Grad2_phi(ptcl,orbital)*psiM(ptcl2,orbital);
              toDet(1,0)+=Phi_alpha(ptcl2,orbital)(dim)*psiM(ptcl,orbital);
              toDet_l(1,0)+=Phi_alpha(ptcl2,orbital)(dim)*psiM(ptcl,orbital);
              toDet(1,1)+=Phi_alpha(ptcl2,orbital)(dim)*psiM(ptcl2,orbital);
              toDet_l(1,1)+=Phi_alpha(ptcl2,orbital)(dim)*psiM(ptcl2,orbital);
            }
            two_row_change(dim,el_dim)+=toDet(0,0)*toDet(1,1)-toDet(1,0)*toDet(0,1);
            if (el_dim==0)
              two_row_change_l(dim)+=toDet_l(0,0)*toDet_l(1,1)-toDet_l(1,0)*toDet_l(0,1);
          }
        }
        Grad_psi_alpha_over_psi(dim,el_dim)=one_row_change(dim,el_dim)+two_row_change(dim,el_dim);
        outfile<<Grad_psi_alpha_over_psi(dim,el_dim)<<endl;
        grad_grad(dim)(ptcl)(el_dim)=one_row_change(dim,el_dim)+two_row_change(dim,el_dim)-
                                     Grad_psi_over_psi(el_dim)*Psi_alpha_over_psi(dim);
      }
    }
    for (int dim=0; dim<OHMMS_DIM; dim++)
    {
      lapl_grad(dim)(ptcl)=0.0;
      lapl_grad(dim)(ptcl)+=one_row_change_l(dim)+two_row_change_l(dim)- Psi_alpha_over_psi(dim)*Grad2_psi_over_psi;
      for (int el_dim=0; el_dim<OHMMS_DIM; el_dim++)
      {
        lapl_grad(dim)(ptcl)-= 2.0*Grad_psi_alpha_over_psi(dim,el_dim)*Grad_psi_over_psi(el_dim);
        lapl_grad(dim)(ptcl)+= 2.0*Psi_alpha_over_psi(dim)*(Grad_psi_over_psi(el_dim)*Grad_psi_over_psi(el_dim));
      }
    }
  }
  outfile.close();
  return Psi_alpha_over_psi;
}

void DiracDeterminantBase::evaluateHessian(ParticleSet& P, HessVector_t& grad_grad_psi)
{
	//IM A HACK.  Assumes evaluateLog has already been executed.
	Phi->evaluate(P, FirstIndex, LastIndex, psiM, dpsiM, grad_grad_source_psiM);
	phi_alpha_Minv = 0.0;
    grad_phi_Minv = 0.0;
    lapl_phi_Minv = 0.0;
    grad_phi_alpha_Minv = 0.0;
    //grad_grad_psi.resize(NumPtcls);

	InvertWithLog(psiM.data(),NumPtcls,NumOrbitals,WorkSpace.data(),Pivot.data(),PhaseValue);
   // InverseTimer.stop();
 //   RatioTimer.start();
 //   for(int i=0, iat=FirstIndex; i<NumPtcls; i++, iat++)
 //   {
 //     GradType rv=simd::dot(psiM[i],dpsiM[i],NumOrbitals);
 //     ValueType lap=simd::dot(psiM[i],d2psiM[i],NumOrbitals);
  //    G(iat) += rv;
  //    L(iat) += lap - dot(rv,rv);
 //   }
 //   RatioTimer.stop();

 

    for(int i=0, iat=FirstIndex; i<NumPtcls; i++, iat++)
    {
      GradType rv=simd::dot(psiM[i],dpsiM[i],NumOrbitals);
      
    //  HessType hess_tmp=simd::dot(psiM[i],grad_grad_source_psiM[i],NumOrbitals);
      HessType hess_tmp;
      hess_tmp=0.0;
      hess_tmp=simd::dot(psiM[i],grad_grad_source_psiM[i],NumOrbitals);
      
   //   for (int orb=0; orb<NumOrbitals; orb++)
  //    {
//		hess_tmp+=psiM(i,orb)*grad_grad_source_psiM(i,orb);
 //     }
     // G(iat) += rv;
     // L(iat) += lap - dot(rv,rv);
  //    app_log()<<"rv = "<<rv<<endl;
   //   app_log()<<"outer = "<<outerProduct(rv,rv)<<endl;
   //   app_log()<<"hess_tmp = "<<hess_tmp<<endl;
      grad_grad_psi[iat]=hess_tmp-outerProduct(rv,rv);
    }
	 psiM_temp = psiM;
	
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
  // 	if (std::fabs(val) > 1.0e-12)
  // 	  cerr << "Error in inverse.\n";
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
            lapl_grad[dim][iel] -= 2.0 * grad_phi_alpha_Minv(j,i)(dim,dim_el)
                                   * grad_phi_Minv(i,j)[dim_el];
        // Third term, eq 9
        // First term, eq 10
        lapl_grad[dim][iel] -= phi_alpha_Minv(j,i)[dim]*lapl_phi_Minv(i,j);
        // Second term, eq 11
        for (int dim_el=0; dim_el<OHMMS_DIM; dim_el++)
          lapl_grad[dim][iel] += 2.0*phi_alpha_Minv(j,i)[dim] *
                                 grad_phi_Minv(i,i)[dim_el]*grad_phi_Minv(i,j)[dim_el];
      }
    }
  }
  return gradPsi;
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
  grad_iat += (1.0/curRatio) * rv;
  RatioTimer.stop();
  return curRatio;
  ////////////////////////////////////////
  ////THIS WILL BE REMOVED. ONLY FOR DEBUG DUE TO WAVEFUNCTIONTEST
  //{
  //  int kat=FirstIndex;
  //  for(int j=0; j<NumOrbitals; j++) {
  //    dpsiM_temp(WorkingIndex,j)=dpsiV[j];
  //    d2psiM_temp(WorkingIndex,j)=d2psiV[j];
  //  }
  //  const ValueType* restrict yptr=psiM_temp.data();
  //  const ValueType* restrict d2yptr=d2psiM_temp.data();
  //  const GradType* restrict dyptr=dpsiM_temp.data();
  //  for(int i=0; i<NumPtcls; i++,kat++) {
  //    //This mimics gemm with loop optimization
  //    GradType rv;
  //    ValueType lap=0.0;
  //    for(int j=0; j<NumOrbitals; j++,yptr++) {
  //      rv += *yptr * *dyptr++;
  //      lap += *yptr * *d2yptr++;
  //    }
  //    lap -= dot(rv,rv);
  //    myG_temp[kat]=rv;
  //    myL_temp[kat]=lap;
  //  }
  //}
  ////////////////////////////////////////
}

/** return the ratio
 * @param P current configuration
 * @param iat particle whose position is moved
 * @param dG differential Gradients
 * @param dL differential Laplacians
 *
 * Data member *_temp contain the data assuming that the move is accepted
 * and are used to evaluate differential Gradients and Laplacians.
 */
DiracDeterminantBase::ValueType DiracDeterminantBase::ratio(ParticleSet& P, int iat,
    ParticleSet::ParticleGradient_t& dG,
    ParticleSet::ParticleLaplacian_t& dL)
{
  UpdateMode=ORB_PBYP_ALL;
  SPOVGLTimer.start();
  Phi->evaluate(P, iat, psiV, dpsiV, d2psiV);
  SPOVGLTimer.stop();
  RatioTimer.start();
  WorkingIndex = iat-FirstIndex;
  //psiM_temp = psiM;
  curRatio= DetRatioByRow(psiM_temp, psiV, WorkingIndex);
  RatioTimer.stop();
  if(abs(curRatio)<numeric_limits<RealType>::epsilon())
  {
    UpdateMode=ORB_PBYP_RATIO; //singularity! do not update inverse
    return 0.0;
  }
  UpdateTimer.start();
  //update psiM_temp with the row substituted
  InverseUpdateByRow(psiM_temp,psiV,workV1,workV2,WorkingIndex,curRatio);
  //update dpsiM_temp and d2psiM_temp
  //std::copy(dpsiV.begin(),dpsiV.end(),dpsiM_temp[WorkingIndex]);
  //std::copy(d2psiV.begin(),d2psiV.end(),d2psiM_temp[WorkingIndex]);
  simd::copy(dpsiM_temp[WorkingIndex],dpsiV.data(),NumOrbitals);
  simd::copy(d2psiM_temp[WorkingIndex],d2psiV.data(),NumOrbitals);
  UpdateTimer.stop();
  RatioTimer.start();
  for(int i=0,kat=FirstIndex; i<NumPtcls; i++,kat++)
  {
    //using inline dot functions
    GradType rv=simd::dot(psiM_temp[i],dpsiM_temp[i],NumOrbitals);
    ValueType lap=simd::dot(psiM_temp[i],d2psiM_temp[i],NumOrbitals);
    lap -= dot(rv,rv);
    dG[kat] += rv - myG[kat];
    myG_temp[kat]=rv;
    dL[kat] += lap -myL[kat];
    myL_temp[kat]=lap;
  }
  RatioTimer.stop();
  return curRatio;
}

DiracDeterminantBase::ValueType DiracDeterminantBase::logRatio(ParticleSet& P, int iat,
    ParticleSet::ParticleGradient_t& dG,
    ParticleSet::ParticleLaplacian_t& dL)
{
  APP_ABORT("  logRatio is not allowed");
  //THIS SHOULD NOT BE CALLED
  ValueType r=ratio(P,iat,dG,dL);
  return LogValue = evaluateLogAndPhase(r,PhaseValue);
}


/** move was accepted, update the real container
*/
void DiracDeterminantBase::acceptMove(ParticleSet& P, int iat)
{
  PhaseValue += evaluatePhase(curRatio);
  LogValue +=std::log(std::abs(curRatio));
  UpdateTimer.start();
  switch(UpdateMode)
  {
  case ORB_PBYP_RATIO:
    InverseUpdateByRow(psiM,psiV,workV1,workV2,WorkingIndex,curRatio);
    break;
  case ORB_PBYP_PARTIAL:
    InverseUpdateByRow(psiM,psiV,workV1,workV2,WorkingIndex,curRatio);
    //std::copy(dpsiV.begin(),dpsiV.end(),dpsiM[WorkingIndex]);
    //std::copy(d2psiV.begin(),d2psiV.end(),d2psiM[WorkingIndex]);
    simd::copy(dpsiM[WorkingIndex],  dpsiV.data(),  NumOrbitals);
    simd::copy(d2psiM[WorkingIndex], d2psiV.data(), NumOrbitals);
    //////////////////////////////////////
    ////THIS WILL BE REMOVED. ONLY FOR DEBUG DUE TO WAVEFUNCTIONTEST
    //myG = myG_temp;
    //myL = myL_temp;
    ///////////////////////
    break;
  default:
    myG = myG_temp;
    myL = myL_temp;
    //psiM = psiM_temp;
    //std::copy(dpsiV.begin(),dpsiV.end(),dpsiM[WorkingIndex]);
    //std::copy(d2psiV.begin(),d2psiV.end(),d2psiM[WorkingIndex]);
    simd::copy(psiM.data(),     psiM_temp.data(),   psiM.size());
    simd::copy(dpsiM[WorkingIndex],  dpsiV.data(),  NumOrbitals);
    simd::copy(d2psiM[WorkingIndex], d2psiV.data(), NumOrbitals);
    break;
  }
  UpdateTimer.stop();
  curRatio=1.0;
}

/** move was rejected. copy the real container to the temporary to move on
*/
void DiracDeterminantBase::restore(int iat)
{
  if(UpdateMode == ORB_PBYP_ALL)
  {
    //psiM_temp = psiM;
    //std::copy(dpsiM[WorkingIndex],dpsiM[WorkingIndex+1],dpsiM_temp[WorkingIndex]);
    //std::copy(d2psiM[WorkingIndex],d2psiM[WorkingIndex+1],d2psiM_temp[WorkingIndex]);
    simd::copy(psiM_temp.data(),     psiM.data(),   psiM.size());
    simd::copy(dpsiM_temp[WorkingIndex],  dpsiM[WorkingIndex],  NumOrbitals);
    simd::copy(d2psiM_temp[WorkingIndex], d2psiM[WorkingIndex], NumOrbitals);
  }
  curRatio=1.0;
}

void DiracDeterminantBase::update(ParticleSet& P,
                                  ParticleSet::ParticleGradient_t& dG,
                                  ParticleSet::ParticleLaplacian_t& dL,
                                  int iat)
{
  UpdateTimer.start();
  InverseUpdateByRow(psiM,psiV,workV1,workV2,WorkingIndex,curRatio);
  //for(int j=0; j<NumOrbitals; j++) {
  //  dpsiM(WorkingIndex,j)=dpsiV[j];
  //  d2psiM(WorkingIndex,j)=d2psiV[j];
  //}
  simd::copy(dpsiM[WorkingIndex],  dpsiV.data(),  NumOrbitals);
  simd::copy(d2psiM[WorkingIndex], d2psiV.data(), NumOrbitals);
  UpdateTimer.stop();
  RatioTimer.start();
  int kat=FirstIndex;
  for(int i=0; i<NumPtcls; i++,kat++)
  {
    GradType rv=simd::dot(psiM[i],dpsiM[i],NumOrbitals);
    ValueType lap=simd::dot(psiM[i],d2psiM[i],NumOrbitals);
    lap -= dot(rv,rv);
    dG[kat] += rv - myG[kat];
    myG[kat]=rv;
    dL[kat] += lap -myL[kat];
    myL[kat]=lap;
  }
  RatioTimer.stop();
  PhaseValue += evaluatePhase(curRatio);
  LogValue +=std::log(std::abs(curRatio));
  curRatio=1.0;
}

DiracDeterminantBase::RealType
DiracDeterminantBase::evaluateLog(ParticleSet& P, PooledData<RealType>& buf)
{
  buf.put(psiM.first_address(),psiM.last_address());
  buf.put(FirstAddressOfdV,LastAddressOfdV);
  buf.put(d2psiM.first_address(),d2psiM.last_address());
  buf.put(myL.first_address(), myL.last_address());
  buf.put(FirstAddressOfG,LastAddressOfG);
  buf.put(LogValue);
  buf.put(PhaseValue);
  return LogValue;
}

void DiracDeterminantBase::registerDataForDerivatives(ParticleSet& P, PooledData<RealType>& buf, int storageType)
{
  DerivStorageType=storageType;
  if(storageType == 0)
    // store everything
  {
    buf.add(psiM.first_address(),psiM.last_address());
//       buf.add(FirstAddressOfdV,LastAddressOfdV);
//       buf.add(d2psiM.first_address(),d2psiM.last_address());
    buf.add(myL.first_address(), myL.last_address());
    buf.add(FirstAddressOfG,LastAddressOfG);
    buf.add(LogValue);
    buf.add(PhaseValue);
  }
  else
    // only store the inverse
  {
    buf.add(psiM.first_address(),psiM.last_address());
  }
}

void DiracDeterminantBase::copyToDerivativeBuffer(ParticleSet& P, PooledData<RealType>& buf)
{
  if(DerivStorageType==0)
  {
    buf.put(psiM.first_address(),psiM.last_address());
//     buf.put(FirstAddressOfdV,LastAddressOfdV);
//     buf.put(d2psiM.first_address(),d2psiM.last_address());
    buf.put(myL.first_address(), myL.last_address());
    buf.put(FirstAddressOfG,LastAddressOfG);
    buf.put(LogValue);
    buf.put(PhaseValue);
  }
  else
  {
    buf.put(psiM.first_address(),psiM.last_address());
  }
}

void DiracDeterminantBase::copyFromDerivativeBuffer(ParticleSet& P, PooledData<RealType>& buf)
{
  if(DerivStorageType==0)
  {
    buf.get(psiM.first_address(),psiM.last_address());
//     buf.get(FirstAddressOfdV,LastAddressOfdV);
//     buf.get(d2psiM.first_address(),d2psiM.last_address());
    buf.get(myL.first_address(), myL.last_address());
    buf.get(FirstAddressOfG,LastAddressOfG);
    buf.get(LogValue);
    buf.get(PhaseValue);
  }
  else
  {
    buf.get(psiM.first_address(),psiM.last_address());
  }
  P.G += myG;
  P.L += myL;
}

DiracDeterminantBase::RealType DiracDeterminantBase::evaluateLogForDerivativeBuffer(ParticleSet& P, PooledData<RealType>& buf)
{
  myG=0.0;
  myL=0.0;
  //ValueType x=evaluate(P,myG,myL);
  LogValue=evaluateLog(P,myG,myL);
  P.G += myG;
  P.L += myL;
  return LogValue;
}

DiracDeterminantBase::RealType DiracDeterminantBase::evaluateLogFromDerivativeBuffer(ParticleSet& P, PooledData<RealType>& buf)
{
//    P.G += myG;
//    P.L += myL;
  return LogValue;
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
//       APP_ABORT("  DiracDeterminantBase::evaluate is disabled");
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
  //      cerr<<"I'm calling evaluate log"<<endl;
  SPOVGLTimer.start();
  Phi->evaluate(P, FirstIndex, LastIndex, psiM,dpsiM, d2psiM);
  SPOVGLTimer.stop();
  if(NumPtcls==1)
  {
    //CurrentDet=psiM(0,0);
    ValueType det=psiM(0,0);
    ValueType y=1.0/det;
    psiM(0,0)=y;
    GradType rv = y*dpsiM(0,0);
    G(FirstIndex) += rv;
    L(FirstIndex) += y*d2psiM(0,0) - dot(rv,rv);
    LogValue = evaluateLogAndPhase(det,PhaseValue);
  }
  else
  {
    InverseTimer.start();
    LogValue=InvertWithLog(psiM.data(),NumPtcls,NumOrbitals,WorkSpace.data(),Pivot.data(),PhaseValue);
    InverseTimer.stop();
    RatioTimer.start();
    for(int i=0, iat=FirstIndex; i<NumPtcls; i++, iat++)
    {
      GradType rv=simd::dot(psiM[i],dpsiM[i],NumOrbitals);
      ValueType lap=simd::dot(psiM[i],d2psiM[i],NumOrbitals);
      G(iat) += rv;
      L(iat) += lap - dot(rv,rv);
    }
    RatioTimer.stop();
  }
  psiM_temp = psiM;
  return LogValue;
}


void
DiracDeterminantBase::evaluateDerivatives(ParticleSet& P,
    const opt_variables_type& active,
    vector<RealType>& dlogpsi,
    vector<RealType>& dhpsioverpsi)
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
  : OrbitalBase(s), NP(0),Phi(s.Phi),FirstIndex(s.FirstIndex)
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
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/

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

#include "QMCWaveFunctions/Fermion/MultiDiracDeterminantBase.h"
#include "Numerics/DeterminantOperators.h"
#include "Numerics/OhmmsBlas.h"

namespace qmcplusplus {

  /** constructor
   *@param spos the single-particle orbital set
   *@param first index of the first particle
   */
  MultiDiracDeterminantBase::MultiDiracDeterminantBase(SPOSetBasePtr const &spos, int first): 
    NP(0), Phi(spos), FirstIndex(first) 
  {
    Optimizable=false;
    OrbitalName="MultiDiracDeterminantBase";
  }

  ///default destructor
  MultiDiracDeterminantBase::~MultiDiracDeterminantBase() {}

  MultiDiracDeterminantBase& MultiDiracDeterminantBase::operator=(const MultiDiracDeterminantBase& s) {
    NP=0;
    resize(s.NumPtcls, s.NumOrbitals);
    return *this;
  }

  /** set the index of the first particle in the determinant and reset the size of the determinant
   *@param first index of first particle
   *@param nel number of particles in the determinant
   */
  void MultiDiracDeterminantBase::set(int first, int nel) {
    FirstIndex = first;
    resize(nel,nel);
  }


  ///reset the size: with the number of particles and number of orbtials
  void MultiDiracDeterminantBase::resize(int nel, int morb) {
    int norb=nel;
    if(norb <= 0) norb = nel; // for morb == -1 (default)
    psiM.resize(nel,norb);

    dpsiM.resize(nel,norb);
    d2psiM.resize(nel,norb);
    psiM_temp.resize(nel,norb);
    dpsiM_temp.resize(nel,norb);
    d2psiM_temp.resize(nel,norb);
    psiMinv.resize(nel,norb);

    WorkSpace.resize(nel);
    Pivot.resize(nel);
    LastIndex = FirstIndex + nel;
    NumPtcls=nel;
    NumOrbitals=norb;
    NumOrbitals_unoccupied=morb-norb;
    NumOrbitals_total=morb;
    psiV.resize(NumOrbitals_total);
    psiM_actual.resize(NumOrbitals_total,nel); //WARNING: Orbitals are first!
    Excitations.set_excitations(10,12,5,9,coefs);
    // For forces
    grad_source_psiM.resize(nel,norb);
    grad_lapl_source_psiM.resize(nel,norb);
    grad_grad_source_psiM.resize(nel,norb);
    phi_alpha_Minv.resize(nel,norb);
    grad_phi_Minv.resize(nel,norb);
    lapl_phi_Minv.resize(nel,norb);
    grad_phi_alpha_Minv.resize(nel,norb);
  }

  MultiDiracDeterminantBase::RealType 
    MultiDiracDeterminantBase::registerData(ParticleSet& P, PooledData<RealType>& buf) 
    {

    if(NP == 0) {//first time, allocate once
      //int norb = cols();
      dpsiV.resize(NumOrbitals_total);
      d2psiV.resize(NumOrbitals_total);
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

  MultiDiracDeterminantBase::RealType MultiDiracDeterminantBase::updateBuffer(ParticleSet& P, 
      PooledData<RealType>& buf, bool fromscratch) 
  {
    ///needs to be fixed to work with tranposed psiM
    myG=0.0;
    myL=0.0;

    if(fromscratch)
      LogValue=evaluateLog(P,myG,myL);
    else
    {
      if(UpdateMode == ORB_PBYP_RATIO) 
	Phi->evaluate(P, FirstIndex, LastIndex, psiM_temp,dpsiM, d2psiM);    

      if(NumPtcls==1) {
        ValueType y=1.0/psiM_temp(0,0);
        psiM(0,0)=y;
        GradType rv = y*dpsiM(0,0);
        myG(FirstIndex) += rv;
        myL(FirstIndex) += y*d2psiM(0,0) - dot(rv,rv);
      } else {
        const ValueType* restrict yptr=psiM.data();
        const ValueType* restrict d2yptr=d2psiM.data();
        const GradType* restrict dyptr=dpsiM.data();
        for(int i=0, iat=FirstIndex; i<NumPtcls; i++, iat++) 
        {
          GradType rv;
          ValueType lap=0.0;
          for(int j=0; j<NumOrbitals; j++,yptr++) {
            rv += *yptr * *dyptr++;
            lap += *yptr * *d2yptr++;
          }
          myG(iat) += rv;
          myL(iat) += lap - dot(rv,rv);
        }
      }
    }

    P.G += myG;
    P.L += myL;

    buf.put(psiM.first_address(),psiM.last_address());
    buf.put(FirstAddressOfdV,LastAddressOfdV);
    buf.put(d2psiM.first_address(),d2psiM.last_address());
    buf.put(myL.first_address(), myL.last_address());
    buf.put(FirstAddressOfG,LastAddressOfG);
    buf.put(LogValue);
    buf.put(PhaseValue);

    return LogValue;
  }

  void MultiDiracDeterminantBase::copyFromBuffer(ParticleSet& P, PooledData<RealType>& buf) {

    buf.get(psiM.first_address(),psiM.last_address());
    buf.get(psiM_actual.first_address(),psiM_actual.last_address());
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
    psiM_temp = psiM;
    dpsiM_temp = dpsiM;
    d2psiM_temp = d2psiM;
  }

  /** dump the inverse to the buffer
  */
  void MultiDiracDeterminantBase::dumpToBuffer(ParticleSet& P, PooledData<RealType>& buf) {
    APP_ABORT("MultiDiracDeterminantBase::dumpToBuffer");
    buf.add(psiM.first_address(),psiM.last_address());
  }

  /** copy the inverse from the buffer
  */
  void MultiDiracDeterminantBase::dumpFromBuffer(ParticleSet& P, PooledData<RealType>& buf) {
    APP_ABORT("MultiDiracDeterminantBase::dumpFromBuffer");
    buf.get(psiM.first_address(),psiM.last_address());
  }

  


  /** return the ratio only for the  iat-th partcle move
   * @param P current configuration
   * @param iat the particle thas is being moved
   */
  MultiDiracDeterminantBase::ValueType MultiDiracDeterminantBase::ratio(ParticleSet& P, int iat) 
  {
    UpdateMode=ORB_PBYP_RATIO;
    WorkingIndex = iat-FirstIndex;
    Phi->evaluate(P, iat, psiV);
#ifdef DIRAC_USE_BLAS
    return curRatio = BLAS::dot(NumOrbitals,psiM[iat-FirstIndex],&psiV[0]);
#else
    return curRatio = DetRatio(psiM, psiV.begin(),iat-FirstIndex);
#endif
  }

  MultiDiracDeterminantBase::GradType 
    MultiDiracDeterminantBase::evalGrad(ParticleSet& P, int iat)
  {
    const ValueType* restrict yptr=psiM[iat-FirstIndex];
    const GradType* restrict dyptr=dpsiM[iat-FirstIndex];
    GradType rv;
    for(int j=0; j<NumOrbitals; ++j) rv += (*yptr++) *(*dyptr++);
    return rv;
  }

  MultiDiracDeterminantBase::GradType 
    MultiDiracDeterminantBase::evalGradSource(ParticleSet& P, ParticleSet& source,
					 int iat)
  {
    Phi->evaluateGradSource (P, FirstIndex, LastIndex,
			     source, iat, grad_source_psiM);
      
//     Phi->evaluate(P, FirstIndex, LastIndex, psiM, dpsiM, d2psiM);
//     LogValue=InvertWithLog(psiM.data(),NumPtcls,NumOrbitals,WorkSpace.data(),Pivot.data(),PhaseValue);

    const ValueType* restrict yptr=psiM[0];
    const GradType* restrict dyptr=grad_source_psiM[0];
    GradType rv(0.0,0.0,0.0);
    for (int i=0; i<NumPtcls; i++)
      for(int j=0; j<NumOrbitals; j++) 
	//rv += (*yptr++) *(*dyptr++);
	rv += grad_source_psiM(i,j) * psiM(i,j);
    // HACK HACK
    //return (grad_source_psiM(1,3));
    return rv;
  }

  MultiDiracDeterminantBase::GradType
  MultiDiracDeterminantBase::evalGradSource
  (ParticleSet& P, ParticleSet& source,int iat,
   TinyVector<ParticleSet::ParticleGradient_t, OHMMS_DIM> &grad_grad,
   TinyVector<ParticleSet::ParticleLaplacian_t,OHMMS_DIM> &lapl_grad)
  {
    Phi->evaluateGradSource (P, FirstIndex, LastIndex, source, iat, 
			     grad_source_psiM, grad_grad_source_psiM, 
			     grad_lapl_source_psiM);
    Phi->evaluate(P, FirstIndex, LastIndex, psiM, dpsiM, d2psiM);
    LogValue=InvertWithLog(psiM.data(),NumPtcls,NumOrbitals,
    			   WorkSpace.data(),Pivot.data(),PhaseValue);

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
    for (int ptcl=0;ptcl<NumPtcls;ptcl++){
      ValueType Grad2_psi_over_psi(0.0);
      GradType Grad_psi_over_psi(0.0);
      HessType Grad_psi_alpha_over_psi(0.0);
      HessType one_row_change(0.0);
      HessType two_row_change(0.0);
      GradType one_row_change_l(0.0);
      GradType two_row_change_l(0.0);

      for (int el_dim=0;el_dim<3;el_dim++){
	for (int orbital=0;orbital<NumOrbitals;orbital++){
	  Grad_psi_over_psi(el_dim)+=Grad_phi(ptcl,orbital)(el_dim)*psiM(ptcl,orbital);
	  if (el_dim==0)
	    Grad2_psi_over_psi+=Grad2_phi(ptcl,orbital)*psiM(ptcl,orbital);
	}
	for (int dim=0;dim<3;dim++){
	  one_row_change(dim,el_dim)=0.0;
	  for (int orbital=0;orbital<NumOrbitals;orbital++){
	    one_row_change(dim,el_dim)+=Grad_phi_alpha(ptcl,orbital)(dim,el_dim)*psiM(ptcl,orbital);
	    if (el_dim==0)
	      one_row_change_l(dim)+=Grad2_phi_alpha(ptcl,orbital)(dim)*psiM(ptcl,orbital);
	  }
	  for (int ptcl2=0;ptcl2<NumPtcls;ptcl2++){
	    if (ptcl!=ptcl2){
	      toDet=0.0;
	      toDet_l=0.0;
	      for (int orbital=0;orbital<NumOrbitals;orbital++){
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
      for (int dim=0;dim<3;dim++){
	lapl_grad(dim)(ptcl)=0.0;
	lapl_grad(dim)(ptcl)+=one_row_change_l(dim)+two_row_change_l(dim)- Psi_alpha_over_psi(dim)*Grad2_psi_over_psi;
	for (int el_dim=0;el_dim<3;el_dim++){
	  lapl_grad(dim)(ptcl)-= 2.0*Grad_psi_alpha_over_psi(dim,el_dim)*Grad_psi_over_psi(el_dim);
	  lapl_grad(dim)(ptcl)+= 2.0*Psi_alpha_over_psi(dim)*(Grad_psi_over_psi(el_dim)*Grad_psi_over_psi(el_dim));
	}
      }
      
    }
    outfile.close();
    return Psi_alpha_over_psi;
  }


  MultiDiracDeterminantBase::GradType
  MultiDiracDeterminantBase::evalGradSourcep
  (ParticleSet& P, ParticleSet& source,int iat,
   TinyVector<ParticleSet::ParticleGradient_t, OHMMS_DIM> &grad_grad,
   TinyVector<ParticleSet::ParticleLaplacian_t,OHMMS_DIM> &lapl_grad)
  {
    Phi->evaluateGradSource (P, FirstIndex, LastIndex, source, iat,
			     grad_source_psiM, grad_grad_source_psiM,
			     grad_lapl_source_psiM);

    // HACK HACK HACK
    // Phi->evaluate(P, FirstIndex, LastIndex, psiM, dpsiM, d2psiM);
    // LogValue=InvertWithLog(psiM.data(),NumPtcls,NumOrbitals,
    // 			   WorkSpace.data(),Pivot.data(),PhaseValue);

    // Compute matrices
    phi_alpha_Minv = 0.0;    grad_phi_Minv = 0.0;
    lapl_phi_Minv = 0.0;     grad_phi_alpha_Minv = 0.0;
    for (int i=0; i<NumPtcls; i++)
      for (int j=0; j<NumOrbitals; j++) {
	lapl_phi_Minv(i,j) = 0.0;
	for (int k=0; k<NumOrbitals; k++)
	  lapl_phi_Minv(i,j) += d2psiM(i,k)*psiM(j,k);
      }
    for (int dim=0; dim<3; dim++) {
      for (int i=0; i<NumPtcls; i++)
	for (int j=0; j<NumOrbitals; j++) {
	  for (int k=0; k<NumOrbitals; k++) {
	    phi_alpha_Minv(i,j)[dim] += grad_source_psiM(i,k)[dim] * psiM(j,k);
	    grad_phi_Minv(i,j)[dim] += dpsiM(i,k)[dim] * psiM(j,k);
	    for (int dim_el=0; dim_el<3; dim_el++)
	      grad_phi_alpha_Minv(i,j)(dim, dim_el) +=
		grad_grad_source_psiM(i,k)(dim,dim_el)*psiM(j,k);
	  }
	}
    }

    GradType gradPsi;
    for(int i=0, iel=FirstIndex; i<NumPtcls; i++, iel++) {
      HessType dval (0.0);
      GradType d2val(0.0);
      for (int dim=0; dim<3; dim++)
	for (int dim_el=0; dim_el<3; dim_el++)
	  dval(dim,dim_el) = grad_phi_alpha_Minv(i,i)(dim,dim_el);
      for(int j=0; j<NumOrbitals; j++) {
	gradPsi += grad_source_psiM(i,j) * psiM(i,j);
	for (int dim=0; dim<3; dim++)
	  for (int k=0; k<3; k++)
	    dval(dim,k) -= phi_alpha_Minv(j,i)[dim]*grad_phi_Minv(i,j)[k];
      }
      for (int dim=0; dim<OHMMS_DIM; dim++) {
	for (int k=0; k<3; k++)
	  grad_grad[dim][iel][k] += dval(dim,k);
	for (int j=0; j<NumOrbitals; j++) {
	  // First term, eq 9
	  lapl_grad[dim][iel] += grad_lapl_source_psiM(i,j)[dim] *
	    psiM(i,j);
	  // Second term, eq 9
	  if (j == i)
	    for (int dim_el=0; dim_el<3; dim_el++)
	      lapl_grad[dim][iel] -= 2.0 * grad_phi_alpha_Minv(j,i)(dim,dim_el)
		* grad_phi_Minv(i,j)[dim_el];
	  // Third term, eq 9
	  // First term, eq 10
	  lapl_grad[dim][iel] -= phi_alpha_Minv(j,i)[dim]*lapl_phi_Minv(i,j);
	  // Second term, eq 11
	  for (int dim_el=0; dim_el<3; dim_el++)
	    lapl_grad[dim][iel] += 2.0*phi_alpha_Minv(j,i)[dim] *
	      grad_phi_Minv(i,i)[dim_el]*grad_phi_Minv(i,j)[dim_el];
	}
      }
    }
    return gradPsi;
  }

    MultiDiracDeterminantBase::ValueType 
      MultiDiracDeterminantBase::ratioGrad(ParticleSet& P, int iat, GradType& grad_iat)
  {
    Phi->evaluate(P, iat, psiV, dpsiV, d2psiV);
    WorkingIndex = iat-FirstIndex;

    UpdateMode=ORB_PBYP_PARTIAL;
    
    const ValueType* restrict vptr = psiV.data();
    const ValueType* restrict yptr = psiM[WorkingIndex];
    const GradType* restrict dyptr = dpsiV.data();
    
    GradType rv;
    curRatio = 0.0;
    for(int j=0; j<NumOrbitals; ++j) {
      rv       += yptr[j] * dyptr[j];
      curRatio += yptr[j] * vptr[j];
    }
    grad_iat += (1.0/curRatio) * rv;
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
  MultiDiracDeterminantBase::ValueType MultiDiracDeterminantBase::ratio(ParticleSet& P, int iat,
      ParticleSet::ParticleGradient_t& dG, 
      ParticleSet::ParticleLaplacian_t& dL) 
  {
    UpdateMode=ORB_PBYP_ALL;
    Phi->evaluate(P, iat, psiV, dpsiV, d2psiV);
    WorkingIndex = iat-FirstIndex;

    //psiM_temp = psiM;
#ifdef DIRAC_USE_BLAS
    curRatio = BLAS::dot(NumOrbitals,psiM_temp[WorkingIndex],&psiV[0]);
#else
    curRatio= DetRatio(psiM_temp, psiV.begin(),WorkingIndex);
#endif

    if(abs(curRatio)<numeric_limits<RealType>::epsilon()) 
    {
      UpdateMode=ORB_PBYP_RATIO; //singularity! do not update inverse 
      return 0.0;
    }

    //update psiM_temp with the row substituted
    DetUpdate(psiM_temp,psiV,workV1,workV2,WorkingIndex,curRatio);

    //update dpsiM_temp and d2psiM_temp 
    for(int j=0; j<NumOrbitals; j++) {
      dpsiM_temp(WorkingIndex,j)=dpsiV[j];
      d2psiM_temp(WorkingIndex,j)=d2psiV[j];
    }

    int kat=FirstIndex;

    const ValueType* restrict yptr=psiM_temp.data();
    const ValueType* restrict d2yptr=d2psiM_temp.data();
    const GradType* restrict dyptr=dpsiM_temp.data();
    for(int i=0; i<NumPtcls; i++,kat++) {
      //This mimics gemm with loop optimization
      GradType rv;
      ValueType lap=0.0;
      for(int j=0; j<NumOrbitals; j++,yptr++) {
        rv += *yptr * *dyptr++;
        lap += *yptr * *d2yptr++;
      }

      //using inline dot functions
      //GradType rv=dot(psiM_temp[i],dpsiM_temp[i],NumOrbitals);
      //ValueType lap=dot(psiM_temp[i],d2psiM_temp[i],NumOrbitals);

      //Old index: This is not pretty
      //GradType rv =psiM_temp(i,0)*dpsiM_temp(i,0);
      //ValueType lap=psiM_temp(i,0)*d2psiM_temp(i,0);
      //for(int j=1; j<NumOrbitals; j++) {
      //  rv += psiM_temp(i,j)*dpsiM_temp(i,j);
      //  lap += psiM_temp(i,j)*d2psiM_temp(i,j);
      //}
      lap -= dot(rv,rv);
      dG[kat] += rv - myG[kat];  myG_temp[kat]=rv;
      dL[kat] += lap -myL[kat];  myL_temp[kat]=lap;
    }

    return curRatio;
  }



  /** return the ratio
   * psiM_temp must be equal to psiM upon entering!
   * @param P current configuration
   * @param iat particle whose position is moved
   * @param dG differential Gradients
   * @param dL differential Laplacians
   *
   * Data member *_temp contain the data assuming that the move is accepted
   * and are used to evaluate differential Gradients and Laplacians.
   */
  MultiDiracDeterminantBase::ValueType MultiDiracDeterminantBase::ratio_multi(ParticleSet& P, int iat)
  {
    UpdateMode=ORB_PBYP_ALL;
    WorkingIndex = iat-FirstIndex;
    Phi->evaluate(P, iat, psiV);
    psiMInv=psiM_temp;
    MatrixOperators::transpose(psiMInv);
`
    Excitations.BuildDotProducts(psiMInv,psiM_actual);
    ValueType oldVal=1.0; //eventually store previous valuess
    int coefIndex=0;
    Excitations.CalcSingleExcitations(coefs,oldVal,coefIndex);
    Excitations.CalcDoubleExcitations(coefs,oldVal,coefIndex);


    ValueType gs_ratio=0.0;
    for (int orbital=0;orbital<NumOrbitals;orbital++)
      gs_ratio+=psiM_temp(WorkingIndex,orbital)*psiV[orbital];

    //relying here on the fact that DetUpdate only uses the first (psiM_temp.cols()) rows of psiV!
    //update psiM_temp with the row substituted
    DetUpdate(psiM_temp,psiV,workV1,workV2,WorkingIndex,gs_ratio);
    psiMInv=psiM_temp;
    MatrixOperators::transpose(psiMInv);
   //currnetly psiM sending incorrectly transposed!

    for (int orbital=0;orbital<NumOrbitals;orbital++)
      psiM_actual(orbital,WorkingIndex)=psiV(orbital);
    Excitations.BuildDotProducts(psiMInv,psiM_actual);
    ValueType val=1.0;
    coefIndex=0;
    Excitations.CalcSingleExcitations(coefs,val,coefIndex);
    Excitations.CalcDoubleExcitations(coefs,val,coefIndex);
    curRatio=(gs_ratio*val)/oldVal;
//     int kat=FirstIndex;

//     const ValueType* restrict yptr=psiM_temp.data();
//     const ValueType* restrict d2yptr=d2psiM_temp.data();
//     const GradType* restrict dyptr=dpsiM_temp.data();
//     for(int i=0; i<NumPtcls; i++,kat++) {
//       //This mimics gemm with loop optimization
//       GradType rv;
//       ValueType lap=0.0;
//       for(int j=0; j<NumOrbitals; j++,yptr++) {
//         rv += *yptr * *dyptr++;
//         lap += *yptr * *d2yptr++;
//       }

//       //using inline dot functions
//       //GradType rv=dot(psiM_temp[i],dpsiM_temp[i],NumOrbitals);
//       //ValueType lap=dot(psiM_temp[i],d2psiM_temp[i],NumOrbitals);

//       //Old index: This is not pretty
//       //GradType rv =psiM_temp(i,0)*dpsiM_temp(i,0);
//       //ValueType lap=psiM_temp(i,0)*d2psiM_temp(i,0);
//       //for(int j=1; j<NumOrbitals; j++) {
//       //  rv += psiM_temp(i,j)*dpsiM_temp(i,j);
//       //  lap += psiM_temp(i,j)*d2psiM_temp(i,j);
//       //}
//       lap -= dot(rv,rv);
//       dG[kat] += rv - myG[kat];  myG_temp[kat]=rv;
//       dL[kat] += lap -myL[kat];  myL_temp[kat]=lap;
//     }

    return curRatio;
  }




  MultiDiracDeterminantBase::ValueType MultiDiracDeterminantBase::logRatio(ParticleSet& P, int iat,
      ParticleSet::ParticleGradient_t& dG, 
      ParticleSet::ParticleLaplacian_t& dL) {
    APP_ABORT("  logRatio is not allowed");
    //THIS SHOULD NOT BE CALLED
    ValueType r=ratio(P,iat,dG,dL);
    return LogValue = evaluateLogAndPhase(r,PhaseValue);
  }


  /** move was accepted, update the real container
  */
  void MultiDiracDeterminantBase::acceptMove(ParticleSet& P, int iat) 
  {
    PhaseValue += evaluatePhase(curRatio);
    LogValue +=std::log(std::abs(curRatio));
    switch(UpdateMode)
    {
      case ORB_PBYP_RATIO:
        DetUpdate(psiM,psiV,workV1,workV2,WorkingIndex,curRatio);
        break;
      case ORB_PBYP_PARTIAL:
        //psiM = psiM_temp;
	DetUpdate(psiM,psiV,workV1,workV2,WorkingIndex,curRatio);
        std::copy(dpsiV.begin(),dpsiV.end(),dpsiM[WorkingIndex]);
        std::copy(d2psiV.begin(),d2psiV.end(),d2psiM[WorkingIndex]);

        //////////////////////////////////////
        ////THIS WILL BE REMOVED. ONLY FOR DEBUG DUE TO WAVEFUNCTIONTEST
        //myG = myG_temp;
        //myL = myL_temp;
        ///////////////////////

        break;
      default:
        myG = myG_temp;
        myL = myL_temp;
        psiM = psiM_temp;
        std::copy(dpsiV.begin(),dpsiV.end(),dpsiM[WorkingIndex]);
        std::copy(d2psiV.begin(),d2psiV.end(),d2psiM[WorkingIndex]);
        break;
    }

    curRatio=1.0;
  }

  /** move was rejected. copy the real container to the temporary to move on
  */
  void MultiDiracDeterminantBase::restore(int iat) {
    if(UpdateMode == ORB_PBYP_ALL) {
      psiM_temp = psiM;
      std::copy(dpsiM[WorkingIndex],dpsiM[WorkingIndex+1],dpsiM_temp[WorkingIndex]);
      std::copy(d2psiM[WorkingIndex],d2psiM[WorkingIndex+1],d2psiM_temp[WorkingIndex]);
    }
    curRatio=1.0;
  }

  void MultiDiracDeterminantBase::update(ParticleSet& P, 
      ParticleSet::ParticleGradient_t& dG, 
      ParticleSet::ParticleLaplacian_t& dL,
      int iat) {
    DetUpdate(psiM,psiV,workV1,workV2,WorkingIndex,curRatio);
    for(int j=0; j<NumOrbitals; j++) {
      dpsiM(WorkingIndex,j)=dpsiV[j];
      d2psiM(WorkingIndex,j)=d2psiV[j];
    }

    int kat=FirstIndex;
    for(int i=0; i<NumPtcls; i++,kat++) {
      GradType rv=dot(psiM[i],dpsiM[i],NumOrbitals);
      ValueType lap=dot(psiM[i],d2psiM[i],NumOrbitals);
      //GradType rv =psiM(i,0)*dpsiM(i,0);
      //ValueType lap=psiM(i,0)*d2psiM(i,0);
      //for(int j=1; j<NumOrbitals; j++) {
      //  rv += psiM(i,j)*dpsiM(i,j);
      //  lap += psiM(i,j)*d2psiM(i,j);
      //}
      lap -= dot(rv,rv);
      dG[kat] += rv - myG[kat]; myG[kat]=rv;
      dL[kat] += lap -myL[kat]; myL[kat]=lap;
    }

    PhaseValue += evaluatePhase(curRatio);
    LogValue +=std::log(std::abs(curRatio));
    curRatio=1.0;
  }

  MultiDiracDeterminantBase::RealType 
    MultiDiracDeterminantBase::evaluateLog(ParticleSet& P, PooledData<RealType>& buf) 
    {
      buf.put(psiM.first_address(),psiM.last_address());
      buf.put(psiM_actual.first_address(),psiM_actual.last_address());
      buf.put(FirstAddressOfdV,LastAddressOfdV);
      buf.put(d2psiM.first_address(),d2psiM.last_address());
      buf.put(myL.first_address(), myL.last_address());
      buf.put(FirstAddressOfG,LastAddressOfG);
      buf.put(LogValue);
      buf.put(PhaseValue);
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
  MultiDiracDeterminantBase::ValueType
    MultiDiracDeterminantBase::evaluate(ParticleSet& P, 
        ParticleSet::ParticleGradient_t& G, 
        ParticleSet::ParticleLaplacian_t& L){

      APP_ABORT("  MultiDiracDeterminantBase::evaluate is distabled");

      Phi->evaluate(P, FirstIndex, LastIndex, psiM,dpsiM, d2psiM);

      ValueType CurrentDet;
      if(NumPtcls==1) {
        CurrentDet=psiM(0,0);
        ValueType y=1.0/CurrentDet;
        psiM(0,0)=y;
        GradType rv = y*dpsiM(0,0);
        G(FirstIndex) += rv;
        L(FirstIndex) += y*d2psiM(0,0) - dot(rv,rv);
      } else {
        CurrentDet = Invert(psiM.data(),NumPtcls,NumOrbitals, WorkSpace.data(), Pivot.data());
        //CurrentDet = Invert(psiM.data(),NumPtcls,NumOrbitals);
        
        const ValueType* restrict yptr=psiM.data();
        const ValueType* restrict d2yptr=d2psiM.data();
        const GradType* restrict dyptr=dpsiM.data();
        for(int i=0, iat=FirstIndex; i<NumPtcls; i++, iat++) {
          GradType rv;
          ValueType lap=0.0;
          for(int j=0; j<NumOrbitals; j++,yptr++) {
            rv += *yptr * *dyptr++;
            lap += *yptr * *d2yptr++;
          }
          //Old index
          //    GradType rv = psiM(i,0)*dpsiM(i,0);
          //    ValueType lap=psiM(i,0)*d2psiM(i,0);
          //    for(int j=1; j<NumOrbitals; j++) {
          //      rv += psiM(i,j)*dpsiM(i,j);
          //      lap += psiM(i,j)*d2psiM(i,j);
          //    }
          G(iat) += rv;
          L(iat) += lap - dot(rv,rv);
        }
      }
      return CurrentDet;
    }


  MultiDiracDeterminantBase::RealType
    MultiDiracDeterminantBase::evaluateLog(ParticleSet& P, 
        ParticleSet::ParticleGradient_t& G, 
        ParticleSet::ParticleLaplacian_t& L)
    {
      cerr<<"I'm calling evaluate log"<<endl;
      //      Phi->evaluate(P, FirstIndex, LastIndex, psiM,dpsiM, d2psiM);
      Phi->evaluate(P, FirstIndex, LastIndex, psiM_actual,dpsiM, d2psiM);

      if(NumPtcls==1) 
      {
	APP_ABORT("Evaluate Log with 1 particle in MultiMultiDiracDeterminantBase is potentially dangerous");
        //CurrentDet=psiM(0,0);
        ValueType det=psiM(0,0);
        ValueType y=1.0/det;
        psiM(0,0)=y;
        GradType rv = y*dpsiM(0,0);
        G(FirstIndex) += rv;
        L(FirstIndex) += y*d2psiM(0,0) - dot(rv,rv);
        LogValue = evaluateLogAndPhase(det,PhaseValue);
      } else {

	///Copying ground state data to psiM
	for (int orbital=0;orbital<NumOrbitals;orbital++)
	  for (int ptcl=0;ptcl<NumPtcls;ptcl++) 
	    psiM(orbital,ptcl)=psiM_actual(orbital,ptcl); //Note: Because psiM is typically an inverse, it's typically ptcl,orbital
	//done copying
	LogValue=InvertWithLog(psiM.data(),NumPtcls,NumOrbitals,WorkSpace.data(),Pivot.data(),PhaseValue);
        const ValueType* restrict yptr=psiM.data();
        const ValueType* restrict d2yptr=d2psiM.data();
        const GradType* restrict dyptr=dpsiM.data();
        for(int i=0, iat=FirstIndex; i<NumPtcls; i++, iat++) {
          GradType rv;
          ValueType lap=0.0;
          for(int j=0; j<NumOrbitals; j++,yptr++) {
            rv += *yptr * *dyptr++;
            lap += *yptr * *d2yptr++;
          }
          G(iat) += rv;
          L(iat) += lap - dot(rv,rv);
        }
      }
      psiM_temp = psiM;
      return LogValue;
    }

  OrbitalBasePtr MultiDiracDeterminantBase::makeClone(ParticleSet& tqp) const
  {
    APP_ABORT(" Cannot use MultiDiracDeterminantBase::makeClone");
    return 0;
    //SPOSetBase* sposclone=Phi->makeClone();
    //MultiDiracDeterminantBase* dclone= new MultiDiracDeterminantBase(sposclone);
    //dclone->set(FirstIndex,LastIndex-FirstIndex);
    //dclone->resetTargetParticleSet(tqp);
    //return dclone;
  }

  MultiDiracDeterminantBase::MultiDiracDeterminantBase(const MultiDiracDeterminantBase& s): 
    OrbitalBase(s), NP(0),Phi(s.Phi),FirstIndex(s.FirstIndex)
  {
    this->resize(s.NumPtcls,s.NumOrbitals);
  }

  SPOSetBasePtr  MultiDiracDeterminantBase::clonePhi() const
  {
    return Phi->makeClone();
  }

}
/***************************************************************************
 * $RCSfile$   $Author: kesler $
 * $Revision: 3582 $   $Date: 2009-02-20 18:07:41 -0600 (Fri, 20 Feb 2009) $
 * $Id: MultiDiracDeterminantBase.cpp 3582 2009-02-21 00:07:41Z kesler $ 
 ***************************************************************************/

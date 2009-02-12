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

namespace qmcplusplus {

  /** constructor
   *@param spos the single-particle orbital set
   *@param first index of the first particle
   */
  DiracDeterminantBase::DiracDeterminantBase(SPOSetBasePtr const &spos, int first): 
    NP(0), Phi(spos), FirstIndex(first) 
  {
    Optimizable=false;
    OrbitalName="DiracDeterminantBase";
  }

  ///default destructor
  DiracDeterminantBase::~DiracDeterminantBase() {}

  DiracDeterminantBase& DiracDeterminantBase::operator=(const DiracDeterminantBase& s) {
    NP=0;
    resize(s.NumPtcls, s.NumOrbitals);
    return *this;
  }

  /** set the index of the first particle in the determinant and reset the size of the determinant
   *@param first index of first particle
   *@param nel number of particles in the determinant
   */
  void DiracDeterminantBase::set(int first, int nel) {
    FirstIndex = first;
    resize(nel,nel);
  }


  ///reset the size: with the number of particles and number of orbtials
  void DiracDeterminantBase::resize(int nel, int morb) {
    int norb=morb;
    if(norb <= 0) norb = nel; // for morb == -1 (default)
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
  }

  DiracDeterminantBase::RealType 
    DiracDeterminantBase::registerData(ParticleSet& P, PooledData<RealType>& buf) 
    {

    if(NP == 0) {//first time, allocate once
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

  void DiracDeterminantBase::copyFromBuffer(ParticleSet& P, PooledData<RealType>& buf) {

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
    psiM_temp = psiM;
    dpsiM_temp = dpsiM;
    d2psiM_temp = d2psiM;
  }

  /** dump the inverse to the buffer
  */
  void DiracDeterminantBase::dumpToBuffer(ParticleSet& P, PooledData<RealType>& buf) {
    APP_ABORT("DiracDeterminantBase::dumpToBuffer");
    buf.add(psiM.first_address(),psiM.last_address());
  }

  /** copy the inverse from the buffer
  */
  void DiracDeterminantBase::dumpFromBuffer(ParticleSet& P, PooledData<RealType>& buf) {
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
    Phi->evaluate(P, iat, psiV);
#ifdef DIRAC_USE_BLAS
    return curRatio = BLAS::dot(NumOrbitals,psiM[iat-FirstIndex],&psiV[0]);
#else
    return curRatio = DetRatio(psiM, psiV.begin(),iat-FirstIndex);
#endif
  }

  DiracDeterminantBase::GradType 
    DiracDeterminantBase::evalGrad(ParticleSet& P, int iat)
  {
    cerr<<"Calling evalGrad"<<endl;
    const ValueType* restrict yptr=psiM[iat-FirstIndex];
    const GradType* restrict dyptr=dpsiM[iat-FirstIndex];
    GradType rv;
    for(int j=0; j<NumOrbitals; ++j) rv += (*yptr++) *(*dyptr++);
    return rv;
  }

  DiracDeterminantBase::GradType 
    DiracDeterminantBase::evalGradSource(ParticleSet& P, ParticleSet& source,
					 int iat)
  {
    Phi->evaluateGradSource (P, FirstIndex, LastIndex,
			     source, iat, grad_source_psiM);
      
    Phi->evaluate(P, FirstIndex, LastIndex, psiM, dpsiM, d2psiM);
    LogValue=InvertWithLog(psiM.data(),NumPtcls,NumOrbitals,WorkSpace.data(),Pivot.data(),PhaseValue);

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


  DiracDeterminantBase::GradType 
  DiracDeterminantBase::evalGradSource
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

    TinyVector<ValueMatrix_t,3> phi_alpha_Minv, grad_phi_Minv;
    for (int dim=0; dim<3; dim++) {
      phi_alpha_Minv[dim].resize(NumPtcls,NumOrbitals);
      grad_phi_Minv[dim].resize(NumPtcls,NumOrbitals);
      for (int i=0; i<NumPtcls; i++)
	for (int j=0; j<NumOrbitals; j++) {
	  phi_alpha_Minv[dim](i,j) = 0.0;
	  grad_phi_Minv[dim](i,j)  = 0.0;
	  for (int k=0; k<NumOrbitals; k++) {
	    phi_alpha_Minv[dim](i,j) += grad_source_psiM(i,k)[dim]
	      * psiM(j,k);
	    grad_phi_Minv[dim](i,j) += dpsiM(i,k)[dim] * psiM(j,k);
	  }
	}
    }

    // TinyVector<ValueMatrix_t,OHMMS_DIM> psi_alpha_inv;
    // for (int i=0; i<3; i++)
    //   psi_alpha_inv[i].resize(NumOrbitals,NumPtcls);

    // for (int dim=0; dim < 3; dim++) {
    //   for (int i=0; i<NumOrbitals; i++)
    // 	for (int j=0; j<NumPtcls; j++)
    // 	  psi_alpha_inv[dim](i,j) = grad_source_psiM(j,i)[dim];
    //   InvertWithLog(psi_alpha_inv[dim].data(),NumPtcls,NumOrbitals,
    // 		    WorkSpace.data(),Pivot.data(),PhaseValue);
    // }

    GradType gradPsi;

    const ValueType* restrict    invptr = psiM[0];
    const GradType*  restrict     yptr  = grad_source_psiM[0];
    const HessType*  restrict    dyptr  = grad_grad_source_psiM[0];
    const GradType*  restrict   d2yptr  = grad_lapl_source_psiM[0];
    const GradType*  restrict  dvalptr  = dpsiM.data();
    const ValueType* restrict d2valptr  = d2psiM.data();

    // for(int i=0; i<NumPtcls; i++) 
    //   for(int j=0; j<NumOrbitals; j++)
    // 	gradPsi += *invptr++ * *yptr++;
    gradPsi = evalGradSource(P, source, iat);
    invptr = psiM[0];

    // for(int i=0, iel=FirstIndex; i<NumPtcls; i++, iel++) {
    //   HessType dval (0.0);
    //   GradType d2val(0.0);
    //   GradType gv(0.0);
    //   ValueType lv(0.0);
    //   GradType dlog_psi(0.0);
    //   for(int j=0; j<NumOrbitals; j++,invptr++) {
    // 	dlog_psi = dlog_psi + *invptr *   *yptr++;

    // 	// dval     = dval     + *invptr *  *dyptr++;
    // 	for (int dim=0; dim<3; dim++)
    // 	  for (int k=0; k<3; k++)
    // 	    dval(dim,k) += psi_alpha_inv[dim](i,j) * (*dyptr)(dim,k);
    // 	dyptr++;

    // 	d2val    = d2val    + *invptr * *d2yptr++;
    // 	gv += *invptr *  *dvalptr++;
    // 	lv += *invptr * *d2valptr++;
    //   }
    //   for (int dim=0; dim<OHMMS_DIM; dim++) {
    // 	grad_grad[dim][iel] += gradPsi[dim] * (dval.getRow(dim) - gv);
    // 	lapl_grad[dim][iel] += d2val[dim];
    //   }
    // }

    for(int i=0, iel=FirstIndex; i<NumPtcls; i++, iel++) {
      HessType dval (0.0);
      GradType d2val(0.0);
      GradType gv(0.0);
      ValueType lv(0.0);
      GradType dlog_psi(0.0);
      for(int j=0; j<NumOrbitals; j++,invptr++) {
	for (int dim=0; dim<3; dim++)
	  for (int k=0; k<3; k++)
	    //dval(dim,k) += psi_alpha_inv[dim](i,j) * grad_grad_source_psiM(i,j)(dim,k);
	    dval(dim,k) += psiM(i,j) * grad_grad_source_psiM(i,j)(dim,k) - 
	      phi_alpha_Minv[dim](j,i)*grad_phi_Minv[k](i,j);

	gv += psiM(i,j) *  dpsiM(i,j);
	dlog_psi += psiM(i,j) * grad_source_psiM(i,j);
      }
      for (int dim=0; dim<OHMMS_DIM; dim++) {
	for (int k=0; k<3; k++) {
	  //grad_grad[dim][iel][k] += gradPsi[dim] * (dval(dim,k) /*- gv[k]*/);
	  // grad_grad[dim][iel][k] = ((i == 1) ? 1.0 : 0.0) *grad_grad_source_psiM(1,3)(dim,k);
	  grad_grad[dim][iel][k] += dval(dim,k);// - gradPsi[dim]*gv[k];
	  //	  grad_grad[dim][iel][k] += ComputeExtraTerms(i,dim,k);  ///BRYAN added extram terms!
	}
	//lapl_grad[dim][iel] += d2val[dim];
      }
    }


    // -0.0061939011

    return gradPsi;
  }

//   double DiracDeterminantBase::ComputeExtraTerms(int ptcl_gradient,int elDim,int ionDim)
//   {
//     ValueMatrix_t det;
//     det.resize(2,2);
//     ValueType toAdd(0);
//     for (int orbital=0;orbital<NumOrbitals;orbital++)
//       det(0,0)+=psiM(ptcl_gradient,orbital)*dpsiM(ptcl_gradient,orbital)[elDim]; //   grad_source_psiM(ptcl_gradient,orbital);
//     for (int ptcl=0;ptcl<NumPtcls;ptcl++){
//       if (ptcl!=ptcl_gradient){
// 	for (int orbital=0;orbital<NumOrbitals;orbital++)
// 	  det(0,1)+=psiM(ptcl_gradient,orbital)*grad_source_psiM(ptcl,orbital)[ionDim];
// 	for (int orbital=0;orbital<NumOrbitals;orbital++)
// 	  det(1,0)+=psiM(ptcl,orbital)*dpsiM(ptcl_gradient,orbital)[elDim];
// 	for (int orbital=0;orbital<NumOrbitals;orbital++)
// 	  det(1,1)+=psiM(ptcl,orbital)*grad_source_psiM(ptcl,orbital)[ionDim];
// 	toAdd+=det(0,0)*det(1,1)-det(1,0)*det(0,1);
//       }
//     }
    
//     return toAdd;
//   }


    DiracDeterminantBase::ValueType 
      DiracDeterminantBase::ratioGrad(ParticleSet& P, int iat, GradType& grad_iat)
  {
    cerr<<"Calling ratio grad"<<endl;
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
  DiracDeterminantBase::ValueType DiracDeterminantBase::ratio(ParticleSet& P, int iat,
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

  DiracDeterminantBase::ValueType DiracDeterminantBase::logRatio(ParticleSet& P, int iat,
      ParticleSet::ParticleGradient_t& dG, 
      ParticleSet::ParticleLaplacian_t& dL) {
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
  void DiracDeterminantBase::restore(int iat) {
    if(UpdateMode == ORB_PBYP_ALL) {
      psiM_temp = psiM;
      std::copy(dpsiM[WorkingIndex],dpsiM[WorkingIndex+1],dpsiM_temp[WorkingIndex]);
      std::copy(d2psiM[WorkingIndex],d2psiM[WorkingIndex+1],d2psiM_temp[WorkingIndex]);
    }
    curRatio=1.0;
  }

  void DiracDeterminantBase::update(ParticleSet& P, 
      ParticleSet::ParticleGradient_t& dG, 
      ParticleSet::ParticleLaplacian_t& dL,
      int iat) {
    cerr<<"UPDATE BEING CALLED"<<endl;
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

  DiracDeterminantBase::RealType 
    DiracDeterminantBase::evaluateLog(ParticleSet& P, PooledData<RealType>& buf) 
    {
      cerr<<"EVALUATE LOG BEING CALLED"<<endl;
      buf.put(psiM.first_address(),psiM.last_address());
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
  DiracDeterminantBase::ValueType
    DiracDeterminantBase::evaluate(ParticleSet& P, 
        ParticleSet::ParticleGradient_t& G, 
        ParticleSet::ParticleLaplacian_t& L){

      APP_ABORT("  DiracDeterminantBase::evaluate is distabled");

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


  DiracDeterminantBase::RealType
    DiracDeterminantBase::evaluateLog(ParticleSet& P, 
        ParticleSet::ParticleGradient_t& G, 
        ParticleSet::ParticleLaplacian_t& L)
    {
      cerr<<"EVALUATELOG SECOND IS BEING CALLED"<<endl;
      Phi->evaluate(P, FirstIndex, LastIndex, psiM,dpsiM, d2psiM);

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
      } else {
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

  OrbitalBasePtr DiracDeterminantBase::makeClone(ParticleSet& tqp) const
  {
    APP_ABORT(" Cannot use DiracDeterminantBase::makeClone");
    return 0;
    //SPOSetBase* sposclone=Phi->makeClone();
    //DiracDeterminantBase* dclone= new DiracDeterminantBase(sposclone);
    //dclone->set(FirstIndex,LastIndex-FirstIndex);
    //dclone->resetTargetParticleSet(tqp);
    //return dclone;
  }

  DiracDeterminantBase::DiracDeterminantBase(const DiracDeterminantBase& s): 
    OrbitalBase(s), NP(0),Phi(s.Phi),FirstIndex(s.FirstIndex)
  {
    this->resize(s.NumPtcls,s.NumOrbitals);
  }

  SPOSetBasePtr  DiracDeterminantBase::clonePhi() const
  {
    return Phi->makeClone();
  }

}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

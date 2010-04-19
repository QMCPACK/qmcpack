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

#include "QMCWaveFunctions/Fermion/DiracDeterminantWithBackflow.h"
#include "QMCWaveFunctions/Fermion/BackflowTransformation.h"
#include "Numerics/DeterminantOperators.h"
#include "Numerics/OhmmsBlas.h"
#include "Numerics/MatrixOperators.h"

namespace qmcplusplus {

  /** constructor
   *@param spos the single-particle orbital set
   *@param first index of the first particle
   */
  DiracDeterminantWithBackflow::DiracDeterminantWithBackflow(SPOSetBasePtr const &spos, BackflowTransformation * BF, int first): DiracDeterminantBase(spos,first)
  {
    Optimizable=true;
    OrbitalName="DiracDeterminantWithBackflow";
    registerTimers();
    BFTrans=BF;
  }

  ///default destructor
  DiracDeterminantWithBackflow::~DiracDeterminantWithBackflow() {}

  DiracDeterminantWithBackflow& DiracDeterminantWithBackflow::operator=(const DiracDeterminantWithBackflow& s) {
    NP=0;
    resize(s.NumPtcls, s.NumOrbitals);
    return *this;
  }


  ///reset the size: with the number of particles and number of orbtials
  void DiracDeterminantWithBackflow::resize(int nel, int morb) {
    int norb=morb;
    if(norb <= 0) norb = nel; // for morb == -1 (default)
    psiM.resize(nel,norb);
    dpsiM.resize(nel,norb);
    //d2psiM.resize(nel,norb);
    //psiM_temp.resize(nel,norb);
    dpsiM_temp.resize(nel,norb);
    //d2psiM_temp.resize(nel,norb);
    psiMinv.resize(nel,norb);
    //psiV.resize(norb);
    WorkSpace.resize(nel);
    Pivot.resize(nel);
    LastIndex = FirstIndex + nel;
    NumPtcls=nel;
    NumOrbitals=norb;
    grad_grad_psiM.resize(nel,norb);
    Qmat.resize(nel);

    // For forces
    /*  not used
    grad_source_psiM.resize(nel,norb);
    grad_lapl_source_psiM.resize(nel,norb);
    grad_grad_source_psiM.resize(nel,norb);
    phi_alpha_Minv.resize(nel,norb);
    grad_phi_Minv.resize(nel,norb);
    lapl_phi_Minv.resize(nel,norb);
    grad_phi_alpha_Minv.resize(nel,norb);
    */
  }

  DiracDeterminantWithBackflow::RealType 
    DiracDeterminantWithBackflow::registerData(ParticleSet& P, PooledData<RealType>& buf) 
    {

    if(NP == 0) {//first time, allocate once
      //int norb = cols();
      NP=P.getTotalNum();
      myG.resize(NP);
      myL.resize(NP);
      myG_temp.resize(NP);
      myL_temp.resize(NP);
      resize(NumPtcls,NumOrbitals);
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

  DiracDeterminantWithBackflow::RealType DiracDeterminantWithBackflow::updateBuffer(ParticleSet& P, 
      PooledData<RealType>& buf, bool fromscratch) 
  {
    myG=0.0;
    myL=0.0;

    UpdateTimer.start();

    LogValue=evaluateLog(P,myG,myL);

    P.G += myG;
    P.L += myL;

    //copy psiM to psiM_temp
    psiM_temp=psiM;

    buf.put(psiM.first_address(),psiM.last_address());
    buf.put(FirstAddressOfdV,LastAddressOfdV);
    buf.put(d2psiM.first_address(),d2psiM.last_address());
    buf.put(myL.first_address(), myL.last_address());
    buf.put(FirstAddressOfG,LastAddressOfG);
    buf.put(LogValue);
    buf.put(PhaseValue);

    UpdateTimer.stop();
    return LogValue;
  }

  void DiracDeterminantWithBackflow::copyFromBuffer(ParticleSet& P, PooledData<RealType>& buf) {

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
  }

  /** return the ratio only for the  iat-th partcle move
   * @param P current configuration
   * @param iat the particle thas is being moved
   */
  DiracDeterminantWithBackflow::ValueType DiracDeterminantWithBackflow::ratio(ParticleSet& P, int iat) 
  {
    RealType OldLog = LogValue;
    RealType OldPhase = PhaseValue;

    Phi->evaluate(BFTrans->QP, FirstIndex, LastIndex, psiM,dpsiM,grad_grad_psiM);
    //app_log() <<psiM <<endl;

    //std::copy(psiM.begin(),psiM.end(),psiMinv.begin());
    psiMinv=psiM;

    // invert backflow matrix
    InverseTimer.start();
    LogValue=InvertWithLog(psiMinv.data(),NumPtcls,NumOrbitals,WorkSpace.data(),Pivot.data(),PhaseValue);
    InverseTimer.stop();
   
    return LogValue/OldLog; 
  }

  DiracDeterminantWithBackflow::GradType 
    DiracDeterminantWithBackflow::evalGrad(ParticleSet& P, int iat)
  {
     APP_ABORT(" Need to implement DiracDeterminantWithBackflow::ratio(ParticleSet& P, int iat). \n");
     return 0.0;
  }

  DiracDeterminantWithBackflow::GradType 
    DiracDeterminantWithBackflow::evalGradSource(ParticleSet& P, ParticleSet& source,
					 int iat)
  {
     APP_ABORT(" Need to implement DiracDeterminantWithBackflow::evalGradSource() \n");
     return 0.0;
  }

  DiracDeterminantWithBackflow::GradType
  DiracDeterminantWithBackflow::evalGradSourcep
  (ParticleSet& P, ParticleSet& source,int iat,
   TinyVector<ParticleSet::ParticleGradient_t, OHMMS_DIM> &grad_grad,
   TinyVector<ParticleSet::ParticleLaplacian_t,OHMMS_DIM> &lapl_grad)
  {
     APP_ABORT(" Need to implement DiracDeterminantWithBackflow::evalGradSourcep() \n");
     return 0.0;
  }


  DiracDeterminantWithBackflow::GradType
  DiracDeterminantWithBackflow::evalGradSource
  (ParticleSet& P, ParticleSet& source,int iat,
   TinyVector<ParticleSet::ParticleGradient_t, OHMMS_DIM> &grad_grad,
   TinyVector<ParticleSet::ParticleLaplacian_t,OHMMS_DIM> &lapl_grad)
  {
     APP_ABORT(" Need to implement DiracDeterminantWithBackflow::evalGradSource() \n");
     return 0.0;
  }

    DiracDeterminantWithBackflow::ValueType 
      DiracDeterminantWithBackflow::ratioGrad(ParticleSet& P, int iat, GradType& grad_iat)
  {
     APP_ABORT(" Need to implement DiracDeterminantWithBackflow::ratioGrad() \n");
     return 0.0;
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
  DiracDeterminantWithBackflow::ValueType DiracDeterminantWithBackflow::ratio(ParticleSet& P, int iat,
      ParticleSet::ParticleGradient_t& dG, 
      ParticleSet::ParticleLaplacian_t& dL) 
  { 
     APP_ABORT(" Need to implement DiracDeterminantWithBackflow::ratio() \n");
     return 0.0;
  }


  DiracDeterminantWithBackflow::RealType
    DiracDeterminantWithBackflow::evaluateLog(ParticleSet& P,
        ParticleSet::ParticleGradient_t& G,
        ParticleSet::ParticleLaplacian_t& L)
  {
    // calculate backflow matrix, 1st and 2nd derivatives
// do I need no transpose here???
    Phi->evaluate(BFTrans->QP, FirstIndex, LastIndex, psiM,dpsiM,grad_grad_psiM);
    //app_log() <<psiM <<endl;

    //std::copy(psiM.begin(),psiM.end(),psiMinv.begin());
    psiMinv=psiM;

    // invert backflow matrix
    InverseTimer.start();
    LogValue=InvertWithLog(psiMinv.data(),NumPtcls,NumOrbitals,WorkSpace.data(),Pivot.data(),PhaseValue);
    InverseTimer.stop();

    // calculate F matrix (gradients wrt bf coordinates)
    // could use dgemv with increments of 3*nCols  
    for(int i=0; i<NumPtcls; i++)
    for(int j=0; j<NumPtcls; j++)
    {
       dpsiM_temp(i,j)=dot(psiMinv[i],dpsiM[j],NumOrbitals);
    }
    //for(int i=0, iat=FirstIndex; i<NumPtcls; i++, iat++)
    // G(iat) += dpsiM_temp(i,i);

    // calculate gradients and first piece of laplacians 
    GradType temp;
    ValueType temp2;
    int num = P.getTotalNum();
    for(int i=0; i<num; i++) {
      temp=0.0;
      temp2=0.0;
      for(int j=0; j<NumPtcls; j++) {
        HessType *ptr1 = &(BFTrans->Amat(i,FirstIndex+j)); 
        ValueType *ptr2 = dpsiM_temp(j,j).data(); 
        ValueType *ptr3 = BFTrans->Bmat_full(i,FirstIndex+j).data();
        for(int la=0; la<3; la++) { 
          temp2 += *(ptr3+la) * (*(ptr2+la));
          for(int lb=0; lb<3; lb++)  {
            temp[la] += (*ptr1)(la,lb) * (*(ptr2+lb));
          }   
        }
      } 
      //app_log() <<"i, Aii, inv(i,i), temp: "
      //          <<i <<"  "
      //          <<BFTrans->Amat(i,i) <<"  "
      //          <<dpsiM_temp(i,i) <<"  " 
      //          <<temp <<endl;
      G(i) += temp; 
      L(i) += temp2;
    }


    // calculate temporary matrix
   for(int j=0; j<NumPtcls; j++) {
     ValueType *ptr = Qmat[j].data();
     for(int la=0; la<9; la++)
       *(ptr+la) = 0.0;
     for(int k=0; k<NumPtcls; k++) {
       ValueType *ptr2 = grad_grad_psiM(j,k).data();
       for(int la=0; la<9; la++)
         *(ptr+la) += *(ptr2+la) * psiMinv(j,k);  
     }
   } 

    // calculate laplacians  // FIX FIX FIX
   for(int j=0; j<NumPtcls; j++) {
    HessType *ptr5 = &(Qmat(j));
    for(int k=0; k<NumPtcls; k++) {
      ValueType *ptr3 = dpsiM_temp(k,j).data();
      ValueType *ptr4 = dpsiM_temp(j,k).data();
      if(k == j) {
        for(int i=0; i<num; i++) {
          HessType *ptr1 = &(BFTrans->Amat(i,FirstIndex+j));
          HessType *ptr2 = &(BFTrans->Amat(i,FirstIndex+k));
          for(int lb=0; lb<3; lb++)
           for(int lc=0; lc<3; lc++) {
             ValueType temp=(*(ptr3+lb))*(*(ptr4+lc))-(*(ptr5))(lb,lc);
             for(int la=0; la<3; la++) 
               L(i) -= temp *  
                       ( (*ptr1)(la,lb) ) * 
                       ( (*ptr2)(la,lc) );
           }
        }
      } else { // k != j
        for(int i=0; i<num; i++) {
          HessType *ptr1 = &(BFTrans->Amat(i,FirstIndex+j));
          HessType *ptr2 = &(BFTrans->Amat(i,FirstIndex+k));
          for(int lb=0; lb<3; lb++)
           for(int lc=0; lc<3; lc++) {
             ValueType temp=(*(ptr3+lb))*(*(ptr4+lc));
             for(int la=0; la<3; la++) 
               L(i) -= temp *
                       ( (*ptr1)(la,lb) ) *
                       ( (*ptr2)(la,lc) );
           }
        }
      } // k==j
    }  // k
   }   // j

    return LogValue;
  }

  DiracDeterminantWithBackflow::ValueType DiracDeterminantWithBackflow::logRatio(ParticleSet& P, int iat,
      ParticleSet::ParticleGradient_t& dG, 
      ParticleSet::ParticleLaplacian_t& dL) {
    APP_ABORT("  logRatio is not allowed");
    //THIS SHOULD NOT BE CALLED
    ValueType r=ratio(P,iat,dG,dL);
    return LogValue = evaluateLogAndPhase(r,PhaseValue);
  }


  /** move was accepted, update the real container
  */
  void DiracDeterminantWithBackflow::acceptMove(ParticleSet& P, int iat) 
  {
     APP_ABORT(" Need to implement DiracDeterminantWithBackflow::acceptMove() \n");
  }

  /** move was rejected. Nothing to restore for now. 
  */
  void DiracDeterminantWithBackflow::restore(int iat) {
     curRatio=1.0;
  }

  void DiracDeterminantWithBackflow::update(ParticleSet& P, 
      ParticleSet::ParticleGradient_t& dG, 
      ParticleSet::ParticleLaplacian_t& dL,
      int iat) {
      
     APP_ABORT(" Need to implement DiracDeterminantWithBackflow::update() \n");

  }


  
  DiracDeterminantWithBackflow::RealType
    DiracDeterminantWithBackflow::evaluateLog(ParticleSet& P, PooledData<RealType>& buf)
    {

     APP_ABORT(" Need to implement DiracDeterminantWithBackflow::evaluateLog(PooldedData)() \n");
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
  DiracDeterminantWithBackflow::ValueType
    DiracDeterminantWithBackflow::evaluate(ParticleSet& P, 
        ParticleSet::ParticleGradient_t& G, 
        ParticleSet::ParticleLaplacian_t& L)
    {

//       APP_ABORT("  DiracDeterminantWithBackflow::evaluate is disabled");

      ValueType logval = evaluateLog(P, G, L);
      return std::cos(PhaseValue)*std::exp(logval);
    }

  void
  DiracDeterminantWithBackflow::evaluateDerivatives(ParticleSet& P,
                                            const opt_variables_type& active,
                                            vector<RealType>& dlogpsi,
                                            vector<RealType>& dhpsioverpsi)
  {

  }


  OrbitalBasePtr DiracDeterminantWithBackflow::makeClone(ParticleSet& tqp) const
  {
    APP_ABORT(" Illegal action. Cannot use DiracDeterminantWithBackflow::makeClone");
    return 0;
  }

  DiracDeterminantWithBackflow* DiracDeterminantWithBackflow::makeCopy(SPOSetBasePtr spo) const
  {
    BackflowTransformation *BF = BFTrans->makeClone(); 
    DiracDeterminantWithBackflow* dclone= new DiracDeterminantWithBackflow(spo,BF);
    dclone->set(FirstIndex,LastIndex-FirstIndex);
    return dclone;
  }

  DiracDeterminantWithBackflow::DiracDeterminantWithBackflow(const DiracDeterminantWithBackflow& s): 
    DiracDeterminantBase(s),BFTrans(s.BFTrans)
  {
    registerTimers();
    this->resize(s.NumPtcls,s.NumOrbitals);
  }

  //SPOSetBasePtr  DiracDeterminantWithBackflow::clonePhi() const
  //{
  //  return Phi->makelone();
  //}

}
/***************************************************************************
 * $RCSfile$   $Author: jeongnim.kim $
 * $Revision: 4710 $   $Date: 2010-03-07 18:05:22 -0600 (Sun, 07 Mar 2010) $
 * $Id: DiracDeterminantWithBackflow.cpp 4710 2010-03-08 00:05:22Z jeongnim.kim $ 
 ***************************************************************************/

//////////////////////////////////////////////////////////////////
// (c) Copyright 2006- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#include "QMCWaveFunctions/AGPDeterminant.h"
#include "Numerics/DeterminantOperators.h"
#include "Numerics/OhmmsBlas.h"
#include "Numerics/MatrixOperators.h"

namespace qmcplusplus {

  AGPDeterminant::AGPDeterminant(BasisSetType* bs): GeminalBasis(bs), NumPtcls(0){
  }
  AGPDeterminant::~AGPDeterminant() {}

  void AGPDeterminant::resize(int nup, int ndown) {
    BasisSize=GeminalBasis->size();

    if(NumPtcls == 0) { //use Numptcls to ensure resizing only once
      Lambda.resize(BasisSize,BasisSize);
      if(nup>ndown) {
        LambdaUP.resize(nup-ndown,BasisSize);
      }

      Nup=nup;
      Ndown=ndown;
      NumPtcls=nup+ndown;

      psiM.resize(nup,nup);
      psiM_temp.resize(nup,nup);
      psiMinv.resize(nup,nup);

      psiU.resize(nup);
      psiD.resize(ndown);
      dpsiU.resize(nup);
      dpsiD.resize(ndown);
      d2psiU.resize(nup);
      d2psiD.resize(ndown);
      WorkSpace.resize(nup);
      Pivot.resize(nup);
      phiT.resize(NumPtcls,BasisSize);
    }

    cout << "<<<< checking size nup, ndown, basis " << nup << " " << ndown 
      << " " << BasisSize << endl;
  }
    void AGPDeterminant::resetTargetParticleSet(ParticleSet& P) { 
      cout << "Reseting target particle " << P.getName() << endl;
      GeminalBasis->resetTargetParticleSet(P);
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
  AGPDeterminant::ValueType
  AGPDeterminant::evaluate(ParticleSet& P, 
      ParticleSet::ParticleGradient_t& G, 
      ParticleSet::ParticleLaplacian_t& L){

    GeminalBasis->evaluate(P);

    /* evaluate psi_up(iat)= \sum_{j} C_{ij} \phi_j^{u}(r_{iat}) 
     * psi_down(iat-Nup) =  \sum_{j} C_{ij} \phi_j^{d}(r_{iat})
     */
    MatrixOperators::product(GeminalBasis->Y, Lambda, phiT);

    //psiM=0.0;
    for(int u=0; u<Nup; u++) {
      //paired block
      for(int d=0, jat=Nup; d<Ndown; d++,jat++) {
        psiM(d,u) = BLAS::dot(BasisSize,phiT[u],GeminalBasis->y(jat));
      }
      //unpaired block Ndown x unpaired
      for(int d=Ndown,unpaired=0; d<Nup; d++,unpaired++) {
        psiM(d,u) = BLAS::dot(BasisSize,LambdaUP[unpaired],GeminalBasis->y(u));
      }
    }

    CurrentDet = Invert(psiM.data(),Nup,Nup,WorkSpace.data(),Pivot.data());

    for(int iat=0; iat<Nup; iat++) {
      GradType rv;
      ValueType lap=0;
      int jat=Nup;
      for(int d=0; d<Ndown; d++,jat++) {
        ValueType dfac=psiM(iat,d);
        rv += dfac*dot(phiT[jat],GeminalBasis->dy(iat),BasisSize);
        lap += dfac*dot(phiT[jat],GeminalBasis->d2y(iat),BasisSize);
      }
      for(int d=Ndown,unpaired=0; d<Nup; d++,unpaired++) {
        ValueType dfac=psiM(iat,d);
        rv += dfac*dot(LambdaUP[unpaired],GeminalBasis->dy(iat),BasisSize);
        lap += dfac*dot(LambdaUP[unpaired],GeminalBasis->d2y(iat),BasisSize);
      }
      G(iat) += rv;
      L(iat) += lap-dot(rv,rv);
    }

    for(int jat=Nup,d=0; jat<NumPtcls; jat++,d++) {
      GradType rv;
      ValueType lap=0;
      for(int u=0; u<Nup; u++) {
        ValueType dfac=psiM(u,d);
        rv += dfac*dot(phiT[u],GeminalBasis->dy(jat),BasisSize);
        lap += dfac*dot(phiT[u],GeminalBasis->d2y(jat),BasisSize);
      }
      G(jat) += rv;
      L(jat) += lap-dot(rv,rv);
    }

    return CurrentDet;
  }

  void
  AGPDeterminant::evaluateLogAndStore(ParticleSet& P) {

    GeminalBasis->evaluate(P);

    /* evaluate psi_up(iat)= \sum_{j} C_{ij} \phi_j^{u}(r_{iat}) 
     * psi_down(iat-Nup) =  \sum_{j} C_{ij} \phi_j^{d}(r_{iat})
     */
    MatrixOperators::product(GeminalBasis->Y, Lambda, phiT);

    for(int u=0; u<Nup; u++) {
      //paired block
      for(int d=0, jat=Nup; d<Ndown; d++,jat++) {
        psiM(d,u) = BLAS::dot(BasisSize,phiT[u],GeminalBasis->y(jat));
      }
      //unpaired block Ndown x unpaired
      for(int d=Ndown,unpaired=0; d<Nup; d++,unpaired++) {
        psiM(d,u) = BLAS::dot(BasisSize,LambdaUP[unpaired],GeminalBasis->y(u));
      }
    }

    CurrentDet = Invert(psiM.data(),Nup,Nup,WorkSpace.data(),Pivot.data());

    for(int iat=0; iat<Nup; iat++) {
      GradType rv;
      ValueType lap=0;
      int jat=Nup;
      for(int d=0; d<Ndown; d++,jat++) {
        ValueType dfac=psiM(iat,d);
        rv += dfac*dot(phiT[jat],GeminalBasis->dy(iat),BasisSize);
        lap += dfac*dot(phiT[jat],GeminalBasis->d2y(iat),BasisSize);
      }
      for(int d=Ndown,unpaired=0; d<Nup; d++,unpaired++) {
        ValueType dfac=psiM(iat,d);
        rv += dfac*dot(LambdaUP[unpaired],GeminalBasis->dy(iat),BasisSize);
        lap += dfac*dot(LambdaUP[unpaired],GeminalBasis->d2y(iat),BasisSize);
      }
      myG[iat]=rv;
      myL[iat]=lap-dot(rv,rv);
    }

    for(int jat=Nup,d=0; jat<NumPtcls; jat++,d++) {
      GradType rv;
      ValueType lap=0;
      for(int u=0; u<Nup; u++) {
        ValueType dfac=psiM(u,d);
        rv += dfac*dot(phiT[u],GeminalBasis->dy(jat),BasisSize);
        lap += dfac*dot(phiT[u],GeminalBasis->d2y(jat),BasisSize);
      }
      myG[jat]=rv;
      myL[jat]=lap-dot(rv,rv);
    }
  }

  AGPDeterminant::ValueType 
  AGPDeterminant::registerData(ParticleSet& P, PooledData<RealType>& buf) {

    if(myG.size() == 0) {
      myG.resize(NumPtcls);
      myL.resize(NumPtcls);
      myG_temp.resize(NumPtcls);
      myL_temp.resize(NumPtcls);
      dY.resize(NumPtcls,BasisSize);
      d2Y.resize(NumPtcls,BasisSize);
      phiTv.resize(BasisSize);
      dYv.resize(BasisSize);
      d2Yv.resize(BasisSize);
    }

    evaluateLogAndStore(P);

    dY = GeminalBasis->dY;
    d2Y = GeminalBasis->d2Y;

    P.G += myG;
    P.L += myL;

//      //add the data: determinant, inverse, gradient and laplacians
//      buf.add(psiM.begin(),psiM.end());
//      buf.add(FirstAddressOfdV,LastAddressOfdV);
//      buf.add(d2psiM.begin(),d2psiM.end());
//      buf.add(myL.begin(), myL.end());
//      buf.add(FirstAddressOfG,LastAddressOfG);
//      buf.add(CurrentDet);
//

      SignValue = (CurrentDet<0.0)?-1.0:1.0;
      return LogValue = log(abs(CurrentDet));
    }

  AGPDeterminant::ValueType 
  AGPDeterminant::updateBuffer(ParticleSet& P, PooledData<RealType>& buf) {
    evaluateLogAndStore(P);
    P.G += myG;
    P.L += myL;
//      myG=0.0;
//      myL=0.0;
//      ValueType x=evaluate(P,myG,myL); 
//      P.G += myG;
//      P.L += myL;
//      buf.put(psiM.begin(),psiM.end());
//      buf.put(FirstAddressOfdV,LastAddressOfdV);
//      buf.put(d2psiM.begin(),d2psiM.end());
//      buf.put(myL.begin(), myL.end());
//      buf.put(FirstAddressOfG,LastAddressOfG);
//      buf.put(CurrentDet);

      return CurrentDet;
    }

    void AGPDeterminant::copyFromBuffer(ParticleSet& P, PooledData<RealType>& buf) {
      //copy current inverse of the determinant
      psiM_temp = psiM;
//      buf.get(psiM.begin(),psiM.end());
//      buf.get(FirstAddressOfdV,LastAddressOfdV);
//      buf.get(d2psiM.begin(),d2psiM.end());
//      buf.get(myL.begin(), myL.end());
//      buf.get(FirstAddressOfG,LastAddressOfG);
//      buf.get(CurrentDet);
//
//      //re-evaluate it for testing
//      //Phi.evaluate(P, FirstIndex, LastIndex, psiM, dpsiM, d2psiM);
//      //CurrentDet = Invert(psiM.data(),NumPtcls,NumOrbitals);
//      //need extra copy for gradient/laplacian calculations without updating it
//      psiM_temp = psiM;
//      dpsiM_temp = dpsiM;
//      d2psiM_temp = d2psiM;
    }


    /** return the ratio only for the  iat-th partcle move
     * @param P current configuration
     * @param iat the particle thas is being moved
     */
    AGPDeterminant::ValueType 
    AGPDeterminant::ratio(ParticleSet& P, int iat) {

      GeminalBasis->evaluate(P,iat);

      const ValueType* restrict y_ptr=GeminalBasis->y(0);
      if(iat<Nup) {
        for(int d=0,jat=Nup; d<Ndown; d++,jat++) {
          psiU[d]=BLAS::dot(BasisSize,y_ptr,phiT[jat]);
        }
        //unpaired block Ndown x unpaired
        for(int d=Ndown,unpaired=0; d<Nup; d++,unpaired++) {
          psiU[d] = BLAS::dot(BasisSize,LambdaUP[unpaired],y_ptr);
        }
        return DetRatio(psiM, psiU.data(),iat);
      } else {
        for(int u=0; u<Ndown; u++) {
          psiD[u]=BLAS::dot(BasisSize,y_ptr,phiT[u]);
        }
        for(int u=Ndown; u<Nup; u++) {
          psiD[u]=0.0;
        }
        return DetRatioTranspose(psiM, psiD.data(),iat-Nup);
      }
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
    AGPDeterminant::ValueType 
    AGPDeterminant::ratio(ParticleSet& P, int iat,
		    ParticleSet::ParticleGradient_t& dG, 
		    ParticleSet::ParticleLaplacian_t& dL) {

      //copy the iat-row to temporary vectors, restored when rejected
      std::copy(phiT[iat],phiT[iat]+BasisSize,phiTv.begin());
      std::copy(dY[iat],dY[iat]+BasisSize, dYv.data());
      std::copy(d2Y[iat],d2Y[iat]+BasisSize, d2Yv.data());

      GeminalBasis->evaluateAll(P,iat);

      BLAS::gemv(Lambda.rows(),Lambda.cols(), Lambda.data(), GeminalBasis->y(0), 
          phiT[iat]);
      std::copy(GeminalBasis->dy(0),GeminalBasis->dy(0)+BasisSize, dY[iat]);
      std::copy(GeminalBasis->d2y(0),GeminalBasis->d2y(0)+BasisSize, d2Y[iat]);

      if(iat<Nup)  
        ratioUp(P,iat,dG,dL);
      else
        ratioDown(P,iat,dG,dL);

      for(int kat=0; kat<Nup; kat++) {
        GradType rv;
        ValueType lap=0;
        for(int d=0, jat=Nup; d<Ndown; d++,jat++) {
          ValueType dfac=psiM_temp(kat,d);
          rv += dfac*dot(phiT[jat],dY[kat],BasisSize);
          lap += dfac*dot(phiT[jat],d2Y[kat],BasisSize);
        }
        for(int d=Ndown,unpaired=0; d<Nup; d++,unpaired++) {
          ValueType dfac=psiM_temp(kat,d);
          rv += dfac*dot(LambdaUP[unpaired],dY[kat],BasisSize);
          lap += dfac*dot(LambdaUP[unpaired],d2Y[kat],BasisSize);
        }
        lap -= dot(rv,rv);
        dG[kat] += (rv-myG[kat]); myG_temp[kat]=rv;
        dL[kat] += (lap-myL[kat]); myL_temp[kat]=lap;
      }

      for(int jat=Nup,d=0; jat<NumPtcls; jat++,d++) {
        GradType rv;
        ValueType lap=0;
        for(int u=0; u<Nup; u++) {
          ValueType dfac=psiM_temp(u,d);
          rv += dfac*dot(phiT[u],dY[jat],BasisSize);
          lap += dfac*dot(phiT[u],d2Y[jat],BasisSize);
        }
        lap -= dot(rv,rv);
        dG[jat] +=  (rv-myG[jat]); myG_temp[jat]=rv;
        dL[jat] +=  (lap-myL[jat]); myL_temp[jat]=lap;
      }

      return curRatio;
    }

    void
    AGPDeterminant::ratioUp(ParticleSet& P, int iat,
		    ParticleSet::ParticleGradient_t& dG, 
		    ParticleSet::ParticleLaplacian_t& dL) {
      const ValueType* restrict y_ptr=GeminalBasis->y(0);
      for(int d=0,jat=Nup; d<Ndown; d++,jat++) {
        psiU[d]=BLAS::dot(BasisSize,y_ptr,phiT[jat]);
      }
      //unpaired block Ndown x unpaired
      for(int d=Ndown,unpaired=0; d<Nup; d++,unpaired++) {
        psiU[d] = BLAS::dot(BasisSize,LambdaUP[unpaired],y_ptr);
      }

      curRatio = DetRatio(psiM_temp, psiU.data(),iat);
      DetUpdate(psiM_temp,psiU,workV1,workV2,iat,curRatio);
    }

    void
    AGPDeterminant::ratioDown(ParticleSet& P, int iat,
		    ParticleSet::ParticleGradient_t& dG, 
		    ParticleSet::ParticleLaplacian_t& dL) {
      const ValueType* restrict y_ptr=GeminalBasis->y(0);
      int d=iat-Nup;
      for(int u=0; u<Ndown; u++) {
        psiD[u]=BLAS::dot(BasisSize,y_ptr,phiT[u]);
      }
      for(int u=Ndown; u<Nup; u++) {
        psiD[u]=0.0;
      }

      curRatio = DetRatioTranspose(psiM_temp, psiD.data(),d);
      DetUpdateTranspose(psiM_temp,psiD,workV1,workV2,d,curRatio);
    }


    /** move was accepted, update the real container
     */
    void AGPDeterminant::acceptMove(ParticleSet& P, int iat) {
      CurrentDet *= curRatio;
      myG = myG_temp;
      myL = myL_temp;
      psiM = psiM_temp;
      curRatio=1.0;
    }

    /** move was rejected. copy the real container to the temporary to move on
     */
    void AGPDeterminant::restore(int iat) {
      psiM_temp = psiM;
      std::copy(phiTv.begin(), phiTv.end(),phiT[iat]);
      std::copy(dYv.begin(), dYv.end(), dY[iat]);
      std::copy(d2Yv.begin(), d2Yv.end(), d2Y[iat]);
      curRatio=1.0;
    }

    
    void AGPDeterminant::update(ParticleSet& P, 
		ParticleSet::ParticleGradient_t& dG, 
		ParticleSet::ParticleLaplacian_t& dL,
		int iat) {

//      DetUpdate(psiM,psiV,workV1,workV2,WorkingIndex,curRatio);
//      for(int j=0; j<NumOrbitals; j++) {
//	dpsiM(WorkingIndex,j)=dpsiV[j];
//	d2psiM(WorkingIndex,j)=d2psiV[j];
//      }
//
//      int kat=FirstIndex;
//      for(int i=0; i<NumPtcls; i++,kat++) {
//	PosType rv =psiM(i,0)*dpsiM(i,0);
//	ValueType lap=psiM(i,0)*d2psiM(i,0);
//	for(int j=1; j<NumOrbitals; j++) {
//	  rv += psiM(i,j)*dpsiM(i,j);
//	  lap += psiM(i,j)*d2psiM(i,j);
//	}
//	lap -= dot(rv,rv);
//	dG[kat] += rv - myG[kat]; myG[kat]=rv;
//	dL[kat] += lap -myL[kat]; myL[kat]=lap;
//      }
//
//      //not very useful
//      CurrentDet *= curRatio;
//      curRatio=1.0;
//    }
//
//    ValueType evaluate(ParticleSet& P, PooledData<RealType>& buf) {
//
//      buf.put(psiM.begin(),psiM.end());
//      buf.put(FirstAddressOfdV,LastAddressOfdV);
//      buf.put(d2psiM.begin(),d2psiM.end());
//      buf.put(myL.begin(), myL.end());
//      buf.put(FirstAddressOfG,LastAddressOfG);
//      buf.put(CurrentDet);
//
    }

  AGPDeterminant::ValueType 
  AGPDeterminant::evaluate(ParticleSet& P, PooledData<RealType>& buf) {
    cout << " AGPDeterminant::evaluate(ParticleSet& P, PooledData<RealType>& buf) " << endl;
    return CurrentDet;
  }

  void AGPDeterminant::resizeByWalkers(int nwalkers) {
    //don't know what to do here
  }
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

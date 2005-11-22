//////////////////////////////////////////////////////////////////
// (c) Copyright 2005-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
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
#include "OhmmsPETE/OhmmsMatrix.h"
#include "QMCWaveFunctions/ThreeBodyGeminal.h"
#include "Particle/DistanceTable.h"
#include "Particle/DistanceTableData.h"
#include "Numerics/DeterminantOperators.h"
#include "Numerics/OhmmsBlas.h"
namespace qmcplusplus {

  ThreeBodyGeminal::ThreeBodyGeminal(ParticleSet& ions, ParticleSet& els): 
    CenterRef(ions), GeminalBasis(0) {
    d_table = DistanceTable::getTable(DistanceTable::add(ions,els));
    NumPtcls=els.getTotalNum();
  }

  ThreeBodyGeminal::~ThreeBodyGeminal() {
    //clean up

  }
    ///reset the value of all the Two-Body Jastrow functions
  void ThreeBodyGeminal::reset() {

  }

  OrbitalBase::ValueType 
  ThreeBodyGeminal::evaluateLog(ParticleSet& P,
		                 ParticleSet::ParticleGradient_t& G, 
		                 ParticleSet::ParticleLaplacian_t& L) {
    GeminalBasis->evaluate(P);

    //Rewrite with gemm
    for(int i=0; i<NumPtcls; i++) {
      for(int k=0; k<BasisSize; k++) {
        V(i,k) = BLAS::dot(BasisSize,GeminalBasis->Y[i],Lambda[k]);
      }
    }

    LogValue=ValueType();
    for(int i=0; i< NumPtcls-1; i++) {
      for(int j=i+1; j<NumPtcls; j++) {
        ValueType x= dot(V[j],GeminalBasis->Y[i],BasisSize);
        LogValue += x;
        Uk[i]+= x;
        Uk[j]+= x;
      }
    }

    for(int i=0; i<NumPtcls; i++)  {
      PosType grad=0.0;
      ValueType lap=0.0;
      for(int j=0; j<NumPtcls; j++) {
        if(j==i) continue;
        grad+=dot(V[j],GeminalBasis->dY[i],BasisSize);
        lap+=dot(V[j],GeminalBasis->d2Y[i],BasisSize);
      }
      G(i)+=0.5*grad;
      L(i)+=0.5*lap;
    }

    return LogValue;
  }

  OrbitalBase::ValueType 
  ThreeBodyGeminal::ratio(ParticleSet& P, int iat) {
    GeminalBasis->evaluate(P,iat);
    curVal=0.0;
    for(int j=0; j<NumPtcls; j++) {
      if(j == iat) continue;
      curVal+= dot(V[j],GeminalBasis->y(0),BasisSize);
    }
    //double counting
    curVal*=0.5;
    return exp(curVal-Uk[iat]);
  }

    /** later merge the loop */
  OrbitalBase::ValueType 
  ThreeBodyGeminal::ratio(ParticleSet& P, int iat,
		    ParticleSet::ParticleGradient_t& dG,
		    ParticleSet::ParticleLaplacian_t& dL) {
    return exp(logRatio(P,iat,dG,dL));
  }

    /** later merge the loop */
  OrbitalBase::ValueType 
  ThreeBodyGeminal::logRatio(ParticleSet& P, int iat,
		    ParticleSet::ParticleGradient_t& dG,
		    ParticleSet::ParticleLaplacian_t& dL) {

    GeminalBasis->evaluateAll(P,iat);

    curVal=0.0;
    curLap=0.0;
    curGrad=0.0;
    for(int j=0; j<NumPtcls; j++) {
      if(j == iat) continue;
      curVal+= dot(V[j],GeminalBasis->y(0),BasisSize);
      curGrad += dot(V[j],GeminalBasis->dy(0),BasisSize);
      curLap += dot(V[j],GeminalBasis->d2y(0),BasisSize);
    }
    
    curGrad*=0.5;
    dG[iat] += curGrad-dUk[iat];

    curLap*=0.5;
    dL[iat] += curLap-d2Uk[iat];

    curVal*=0.5;
    return curVal-Uk[iat];
  }

  void ThreeBodyGeminal::restore(int iat) {
    //nothing to do here
  }

  void ThreeBodyGeminal::update(ParticleSet& P, int iat) {
    Uk[iat]=curVal;
    dUk[iat]=curGrad;
    d2Uk[iat]=curLap;
    //change to gemv
    for(int k=0; k<BasisSize; k++) {
      V(iat,k) = BLAS::dot(BasisSize,GeminalBasis->y(0),Lambda[k]);
    }
  }

  void ThreeBodyGeminal::update(ParticleSet& P, 		
		       ParticleSet::ParticleGradient_t& dG, 
		       ParticleSet::ParticleLaplacian_t& dL,
		       int iat) {
    dG[iat]+=curGrad-dUk[iat]; 
    dL[iat]+=curLap-d2Uk[iat]; 
    update(P,iat);
  }

  OrbitalBase::ValueType 
  ThreeBodyGeminal::registerData(ParticleSet& P, PooledData<RealType>& buf) {

    LogValue=evaluateLog(P,P.G,P.L);
    buf.add(V.begin(), V.end());
    return LogValue;
  }

  OrbitalBase::ValueType 
  ThreeBodyGeminal::updateBuffer(ParticleSet& P, 
      PooledData<RealType>& buf) {
    LogValue=evaluateLog(P,P.G,P.L);
    buf.put(V.begin(), V.end());
    return LogValue;
  }
    
  void 
  ThreeBodyGeminal::copyFromBuffer(ParticleSet& P, 
      PooledData<RealType>& buf) {
    buf.get(V.begin(), V.end());
  }

  OrbitalBase::ValueType 
  ThreeBodyGeminal::evaluate(ParticleSet& P, PooledData<RealType>& buf) {
    LogValue=evaluateLog(P,P.G,P.L);
    buf.put(V.begin(), V.end());
    return exp(LogValue);
  }

  template<class T>
    inline bool
    putContent(Matrix<T>& a, xmlNodePtr cur){
      std::istringstream
        stream((const char*)
            (xmlNodeListGetString(cur->doc, cur->xmlChildrenNode, 1)));
      int i=0, ntot=a.size();
      while(!stream.eof() && i<ntot){ stream >> a(i++);}
      return true;
    }


  bool ThreeBodyGeminal::put(xmlNodePtr cur, VarRegistry<RealType>& varlist) {

    BasisSize = GeminalBasis->TotalBasis;

    app_log() << "  The number of Geminal functions "
      <<"for Three-body Jastrow " << BasisSize << endl;
    app_log() << "  The number of particles " << NumPtcls << endl;

    U.resize(NumPtcls,BasisSize);
    V.resize(NumPtcls,BasisSize);
    Lambda.resize(BasisSize,BasisSize);

    //assign the coefficients
    putContent(Lambda,cur);

    //symmetrize it
    //for(int ib=1; ib<BasisSize; ib++) {
    //  for(int jb=0; jb<ib; jb++) {
    //    Lambda(ib,jb) = Lambda(jb,ib);
    //  }
    //}

    Uk.resize(NumPtcls);
    dUk.resize(NumPtcls);
    d2Uk.resize(NumPtcls);

    GeminalBasis->resize(NumPtcls);

    return true;
  }
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/


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
#include "Numerics/MatrixOperators.h"
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

    //only symmetrize it
    for(int ib=0; ib<BasisSize-1; ib++) 
      for(int jb=ib+1; jb<BasisSize; jb++)
        Lambda(jb,ib)=Lambda(ib,jb);
  }

  OrbitalBase::ValueType 
  ThreeBodyGeminal::evaluateLog(ParticleSet& P,
		                 ParticleSet::ParticleGradient_t& G, 
		                 ParticleSet::ParticleLaplacian_t& L) {
    GeminalBasis->evaluate(P);

    MatrixOperators::product(GeminalBasis->Y, Lambda, V);

    //Rewrite with gemm
    //for(int i=0; i<NumPtcls; i++) {
    //  for(int k=0; k<BasisSize; k++) {
    //    V(i,k) = BLAS::dot(BasisSize,GeminalBasis->Y[i],Lambda[k]);
    //  }
    //}
    Uk=0.0;

    LogValue=ValueType();
    for(int i=0; i< NumPtcls-1; i++) {
      const RealType* restrict yptr=GeminalBasis->Y[i];
      for(int j=i+1; j<NumPtcls; j++) {
        ValueType x= dot(V[j],yptr,BasisSize);
        LogValue += x;
        Uk[i]+= x;
        Uk[j]+= x;
      }
    }

    for(int i=0; i<NumPtcls; i++)  {
      const PosType* restrict dptr=GeminalBasis->dY[i];
      const RealType* restrict d2ptr=GeminalBasis->d2Y[i];
      const RealType* restrict vptr=V[0];
      GradType grad(0.0);
      ValueType lap(0.0);
      for(int j=0; j<NumPtcls; j++, vptr+=BasisSize) {
        if(j!=i) {
          grad += dot(vptr,dptr,BasisSize);
          lap +=  dot(vptr,d2ptr,BasisSize);
        }
      }
      G(i)+=grad;
      L(i)+=lap;
    }

    return LogValue;
  }

  OrbitalBase::ValueType 
  ThreeBodyGeminal::ratio(ParticleSet& P, int iat) {
    GeminalBasis->evaluate(P,iat);
    diffVal=0.0;
    for(int j=0; j<NumPtcls; j++) {
      if(j == iat) continue;
      diffVal+= dot(V[j],GeminalBasis->y(0),BasisSize);
    }
    return exp(diffVal-Uk[iat]);
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

    const ValueType* restrict y_ptr=GeminalBasis->y(0);
    const GradType* restrict  dy_ptr=GeminalBasis->dy(0);
    const ValueType* restrict d2y_ptr=GeminalBasis->d2y(0);

    for(int k=0; k<BasisSize; k++) {
      curV[k] = BLAS::dot(BasisSize,y_ptr,Lambda[k]);
      delV[k] = curV[k]-V[iat][k];
    }

    diffVal=0.0;
    GradType dg_acc(0.0);
    ValueType dl_acc(0.0);
    const RealType* restrict vptr=V[0];
    for(int j=0; j<NumPtcls; j++, vptr+=BasisSize) {
      if(j == iat) {
        curLap[j]=0.0;
        curGrad[j]=0.0;
        tLap[j]=0.0;
        tGrad[j]=0.0;
      } else {
        diffVal+= (curVal[j]=dot(delV.data(),Y[j],BasisSize));
        dG[j] += (tGrad[j]=dot(delV.data(),dY[j],BasisSize));
        dL[j] += (tLap[j]=dot(delV.data(),d2Y[j],BasisSize));

        curGrad[j]= dot(vptr,dy_ptr,BasisSize);
        curLap[j] = dot(vptr,d2y_ptr,BasisSize);

        dg_acc += curGrad[j]-dUk(iat,j);
        dl_acc += curLap[j]-d2Uk(iat,j);
      }
    }
    
    dG[iat] += dg_acc;
    dL[iat] += dl_acc;

    curVal[iat]=diffVal;

    return diffVal;
  }

  void ThreeBodyGeminal::restore(int iat) {
    //nothing to do here
  }

  void ThreeBodyGeminal::acceptMove(ParticleSet& P, int iat) {

    //add the differential
    LogValue += diffVal;
    Uk+=curVal; //accumulate the differences

    dUk.replaceRow(curGrad.begin(),iat);
    d2Uk.replaceRow(curLap.begin(),iat);

    dUk.add2Column(tGrad.begin(),iat);
    d2Uk.add2Column(tLap.begin(),iat);

    Y.replaceRow(GeminalBasis->y(0),iat);
    dY.replaceRow(GeminalBasis->dy(0),iat);
    d2Y.replaceRow(GeminalBasis->d2y(0),iat);
    V.replaceRow(curV.begin(),iat);

  }

  void ThreeBodyGeminal::update(ParticleSet& P, 		
		       ParticleSet::ParticleGradient_t& dG, 
		       ParticleSet::ParticleLaplacian_t& dL,
		       int iat) {
    cout << "****  This is to be removed " << endl;
    //dG[iat]+=curGrad-dUk[iat]; 
    //dL[iat]+=curLap-d2Uk[iat]; 
    acceptMove(P,iat);
  }

  OrbitalBase::ValueType 
  ThreeBodyGeminal::registerData(ParticleSet& P, PooledData<RealType>& buf) {

    evaluateLogAndStore(P);
    FirstAddressOfdY=&(dY(0,0)[0]);
    LastAddressOfdY=FirstAddressOfdY+NumPtcls*BasisSize*DIM;

    FirstAddressOfgU=&(dUk(0,0)[0]);
    LastAddressOfgU = FirstAddressOfgU + NumPtcls*NumPtcls*DIM;

    buf.add(LogValue);
    buf.add(V.begin(), V.end());

    buf.add(Y.begin(), Y.end());
    buf.add(FirstAddressOfdY,LastAddressOfdY);
    buf.add(d2Y.begin(),d2Y.end());

    buf.add(Uk.begin(), Uk.end());
    buf.add(FirstAddressOfgU,LastAddressOfgU);
    buf.add(d2Uk.begin(), d2Uk.end());

    return LogValue;
  }

  void 
  ThreeBodyGeminal::evaluateLogAndStore(ParticleSet& P) {
    GeminalBasis->evaluate(P);

    MatrixOperators::product(GeminalBasis->Y, Lambda, V);
    
    Y=GeminalBasis->Y;
    dY=GeminalBasis->dY;
    d2Y=GeminalBasis->d2Y;

    Uk=0.0;
    LogValue=ValueType();
    for(int i=0; i< NumPtcls-1; i++) {
      const RealType* restrict yptr=GeminalBasis->Y[i];
      for(int j=i+1; j<NumPtcls; j++) {
        ValueType x= dot(V[j],yptr,BasisSize);
        LogValue += x;
        Uk[i]+= x;
        Uk[j]+= x;
      }
    }

    for(int i=0; i<NumPtcls; i++)  {
      const PosType* restrict dptr=GeminalBasis->dY[i];
      const RealType* restrict d2ptr=GeminalBasis->d2Y[i];
      const RealType* restrict vptr=V[0];
      GradType grad(0.0);
      ValueType lap(0.0);
      for(int j=0; j<NumPtcls; j++, vptr+=BasisSize) {
        if(j==i) {
          dUk(i,j) = 0.0;
          d2Uk(i,j)= 0.0;
        } else {
          grad+= (dUk(i,j) = dot(vptr,dptr,BasisSize));
          lap += (d2Uk(i,j)= dot(vptr,d2ptr,BasisSize));
        }
      }
      P.G(i)+=grad;
      P.L(i)+=lap;
    }
  }

  void 
  ThreeBodyGeminal::copyFromBuffer(ParticleSet& P, PooledData<RealType>& buf) {
    buf.get(LogValue);
    buf.get(V.begin(), V.end());

    buf.get(Y.begin(), Y.end());
    buf.get(FirstAddressOfdY,LastAddressOfdY);
    buf.get(d2Y.begin(),d2Y.end());

    buf.get(Uk.begin(), Uk.end());
    buf.get(FirstAddressOfgU,LastAddressOfgU);
    buf.get(d2Uk.begin(), d2Uk.end());
  }

  OrbitalBase::ValueType 
  ThreeBodyGeminal::evaluate(ParticleSet& P, PooledData<RealType>& buf) {
    buf.put(LogValue);
    buf.put(V.begin(), V.end());

    buf.put(Y.begin(), Y.end());
    buf.put(FirstAddressOfdY,LastAddressOfdY);
    buf.put(d2Y.begin(),d2Y.end());

    buf.put(Uk.begin(), Uk.end());
    buf.put(FirstAddressOfgU,LastAddressOfgU);
    buf.put(d2Uk.begin(), d2Uk.end());

    return exp(LogValue);
  }

  OrbitalBase::ValueType 
  ThreeBodyGeminal::updateBuffer(ParticleSet& P, PooledData<RealType>& buf) {
    evaluateLogAndStore(P);
    buf.put(LogValue);
    buf.put(V.begin(), V.end());

    buf.put(Y.begin(), Y.end());
    buf.put(FirstAddressOfdY,LastAddressOfdY);
    buf.put(d2Y.begin(),d2Y.end());

    buf.put(Uk.begin(), Uk.end());
    buf.put(FirstAddressOfgU,LastAddressOfgU);
    buf.put(d2Uk.begin(), d2Uk.end());
    return LogValue;
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

    Lambda.resize(BasisSize,BasisSize);

    char coeffname[128];
    string aname("j3g");

    const xmlChar* aptr=xmlGetProp(cur,(const xmlChar*)"name");
    if(aptr!=NULL) aname = (const char*)aptr;

    ////assign the coefficients
    //putContent(Lambda,cur);

    ////symmetrize it
    //for(int ib=0; ib<BasisSize; ib++) {
    //  sprintf(coeffname,"%s_%d_%d",aname.c_str(),ib,ib);
    //  RealType* lptr=Lambda.data()+ib*BasisSize+ib;
    //  varlist.add(coeffname,lptr);
    //  for(int jb=ib+1; jb<BasisSize; jb++) {
    //    Lambda(jb,ib) = Lambda(ib,jb);

    //    ++lptr;
    //    sprintf(coeffname,"%s_%d_%d",aname.c_str(),ib,jb);
    //    varlist.add(coeffname,lptr);
    //  }
    //}
    int offset=1;
    xmlNodePtr tcur=cur->xmlChildrenNode;
    while(tcur != NULL) {
      if(xmlStrEqual(tcur->name,(const xmlChar*)"lambda")) {
        int i=atoi((const char*)(xmlGetProp(tcur,(const xmlChar*)"i")));
        int j=atoi((const char*)(xmlGetProp(tcur,(const xmlChar*)"j")));
        double c=atof((const char*)(xmlGetProp(tcur,(const xmlChar*)"c")));
        Lambda(i-offset,j-offset)=c;
        if(i != j) {
          Lambda(j-offset,i-offset)=c;
        }
        sprintf(coeffname,"%s_%d_%d",aname.c_str(),i,j);
        varlist.add(coeffname,Lambda.data()+i*BasisSize+j);
      }
      tcur=tcur->next;
    }

    V.resize(NumPtcls,BasisSize);
    Y.resize(NumPtcls,BasisSize);
    dY.resize(NumPtcls,BasisSize);
    d2Y.resize(NumPtcls,BasisSize);

    curGrad.resize(NumPtcls);
    curLap.resize(NumPtcls);
    curVal.resize(NumPtcls);

    tGrad.resize(NumPtcls);
    tLap.resize(NumPtcls);
    curV.resize(BasisSize);
    delV.resize(BasisSize);

    Uk.resize(NumPtcls);
    dUk.resize(NumPtcls,NumPtcls);
    d2Uk.resize(NumPtcls,NumPtcls);


    //app_log() << "  Three-body Geminal coefficients " << endl;
    //app_log() << Lambda << endl;

    GeminalBasis->resize(NumPtcls);

    return true;
  }
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/


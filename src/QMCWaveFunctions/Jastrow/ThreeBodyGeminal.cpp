//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#include "OhmmsPETE/OhmmsMatrix.h"
#include "QMCWaveFunctions/Jastrow/ThreeBodyGeminal.h"
#include "Particle/DistanceTable.h"
#include "Particle/DistanceTableData.h"
#include "Numerics/DeterminantOperators.h"
#include "Numerics/OhmmsBlas.h"
#include "Numerics/LibxmlNumericIO.h"
#include "Utilities/RandomGenerator.h"
#include "Numerics/MatrixOperators.h"
#include "OhmmsData/AttributeSet.h"

namespace qmcplusplus
{

ThreeBodyGeminal::ThreeBodyGeminal(const ParticleSet& ions, ParticleSet& els):
  CenterRef(ions), GeminalBasis(0), IndexOffSet(1), ID_Lambda("j3")
{
  d_table = DistanceTable::add(ions,els,DT_AOS);
  NumPtcls=els.getTotalNum();
  NormFac=1.0/static_cast<RealType>(NumPtcls*NumPtcls);
  Optimizable=true;
}

ThreeBodyGeminal::~ThreeBodyGeminal()
{
  //clean up
}

void ThreeBodyGeminal::resetTargetParticleSet(ParticleSet& P)
{
  d_table = DistanceTable::add(CenterRef,P,DT_AOS);
  GeminalBasis->resetTargetParticleSet(P);
}

void ThreeBodyGeminal::checkInVariables(opt_variables_type& active)
{
  active.insertFrom(myVars);
  int ncur=active.size();
  GeminalBasis->checkInVariables(active);
  if(ncur!= active.size())
    GeminalBasis->checkInVariables(myVars);
}

void ThreeBodyGeminal::checkOutVariables(const opt_variables_type& active)
{
  myVars.getIndex(active);
  GeminalBasis->checkOutVariables(active);
  app_log() << "<j3-variables>"<< std::endl;
  myVars.print(app_log());
  app_log() << "</j3-variables>"<< std::endl;
}

///reset the value of all the Two-Body Jastrow functions
void ThreeBodyGeminal::resetParameters(const opt_variables_type& active)
{
  int ii=0; //counter for (i,j) pairs for i<j<BasisSize
  int aii=0;//counter for the variables that are meant to be optimized
  for(int ib=0; ib<BasisSize; ib++)
  {
    if(FreeLambda[ii++])
    {
      int loc=myVars.where(aii++);
      if(loc>=0)
        Lambda(ib,ib)=active[loc];
    }
    for(int jb=ib+1; jb<BasisSize; jb++)
    {
      if(FreeLambda[ii++])
      {
        int loc=myVars.where(aii++);
        if(loc>=0)
          Lambda(jb,ib)=Lambda(ib,jb)=active[loc];
      }
    }
  }
  GeminalBasis->resetParameters(active);
  for(int i=0; i<myVars.size(); ++i)
    if(myVars.where(i)>=0)
      myVars[i]=active[myVars.where(i)];
  //app_log() << "ThreeBodyGeminal::resetParameters" << std::endl;
  //app_log() << Lambda << std::endl;
}

void ThreeBodyGeminal::reportStatus(std::ostream& os)
{
  myVars.print(os);
}

OrbitalBase::RealType
ThreeBodyGeminal::evaluateLog(ParticleSet& P,
                              ParticleSet::ParticleGradient_t& G,
                              ParticleSet::ParticleLaplacian_t& L)
{
  //this is necessary to handle ratio for pseudopotentials
  evaluateLogAndStore(P);
  ////GeminalBasis->evaluate(P);
  //GeminalBasis->evaluateForWalkerMove(P);
  //
  //MatrixOperators::product(GeminalBasis->Y, Lambda, V);
  ////Rewrite with gemm
  ////for(int i=0; i<NumPtcls; i++) {
  ////  for(int k=0; k<BasisSize; k++) {
  ////    V(i,k) = BLAS::dot(BasisSize,GeminalBasis->Y[i],Lambda[k]);
  ////  }
  ////}
  //Uk=0.0;
  //LogValue=ValueType();
  //for(int i=0; i< NumPtcls-1; i++) {
  //  const RealType* restrict yptr=GeminalBasis->Y[i];
  //  for(int j=i+1; j<NumPtcls; j++) {
  //    RealType x= dot(V[j],yptr,BasisSize);
  //    LogValue += x;
  //    Uk[i]+= x;
  //    Uk[j]+= x;
  //  }
  //}
  //for(int i=0; i<NumPtcls; i++)  {
  //  const PosType* restrict dptr=GeminalBasis->dY[i];
  //  const RealType* restrict d2ptr=GeminalBasis->d2Y[i];
  //  const RealType* restrict vptr=V[0];
  //  GradType grad(0.0);
  //  ValueType lap(0.0);
  //  for(int j=0; j<NumPtcls; j++, vptr+=BasisSize) {
  //    if(j!=i) {
  //      grad += dot(vptr,dptr,BasisSize);
  //      lap +=  dot(vptr,d2ptr,BasisSize);
  //    }
  //  }
  //  G(i)+=grad;
  //  L(i)+=lap;
  //}
  return LogValue;
}

OrbitalBase::ValueType
ThreeBodyGeminal::ratio(ParticleSet& P, int iat)
{
  UpdateMode=ORB_PBYP_RATIO;
  GeminalBasis->evaluateForPtclMove(P,iat);
  const BasisSetType::RealType* restrict y_ptr=GeminalBasis->Phi.data();
  for(int k=0; k<BasisSize; k++)
  {
    //curV[k] = BLAS::dot(BasisSize,y_ptr,Lambda[k]);
    curV[k] = simd::dot(y_ptr,Lambda[k],BasisSize);
    delV[k] = curV[k]-V[iat][k];
  }
  diffVal=0.0;
  const RealType* restrict vptr=V[0];
  for(int j=0; j<NumPtcls; j++, vptr+=BasisSize)
  {
    if(j==iat)
      continue;
    diffVal+= (curVal[j]=simd::dot(delV.data(),Y[j],BasisSize));
  }
  curVal[iat]=diffVal;
  return std::exp(diffVal);
}

void ThreeBodyGeminal::restore(int iat)
{
  //nothing to do here
}

void ThreeBodyGeminal::acceptMove(ParticleSet& P, int iat)
{
  //add the differential
  LogValue += diffVal;
  Uk+=curVal; //accumulate the differences
  if(UpdateMode == ORB_PBYP_RATIO)
  {
    Y.replaceRow(GeminalBasis->Phi.data(),iat);
    V.replaceRow(curV.begin(),iat);
  }
  else
  {
    dUk.replaceRow(curGrad.begin(),iat);
    d2Uk.replaceRow(curLap.begin(),iat);
    dUk.add2Column(tGrad.begin(),iat);
    d2Uk.add2Column(tLap.begin(),iat);
    //Y.replaceRow(GeminalBasis->y(0),iat);
    //dY.replaceRow(GeminalBasis->dy(0),iat);
    //d2Y.replaceRow(GeminalBasis->d2y(0),iat);
    Y.replaceRow(GeminalBasis->Phi.data(),iat);
    dY.replaceRow(GeminalBasis->dPhi.data(),iat);
    d2Y.replaceRow(GeminalBasis->d2Phi.data(),iat);
    V.replaceRow(curV.begin(),iat);
  }
}

void
ThreeBodyGeminal::registerData(ParticleSet& P, WFBufferType& buf)
{
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
}

void
ThreeBodyGeminal::evaluateLogAndStore(ParticleSet& P)
{
  GeminalBasis->evaluateForWalkerMove(P);
  MatrixOperators::product(GeminalBasis->Y, Lambda, V);
  Y=GeminalBasis->Y;
  dY=GeminalBasis->dY;
  d2Y=GeminalBasis->d2Y;
  Uk=0.0;
  LogValue=RealType();
  for(int i=0; i< NumPtcls-1; i++)
  {
    const RealType* restrict yptr=GeminalBasis->Y[i];
    for(int j=i+1; j<NumPtcls; j++)
    {
      RealType x= simd::dot(V[j],yptr,BasisSize);
      LogValue += x;
      Uk[i]+= x;
      Uk[j]+= x;
    }
  }
  for(int i=0; i<NumPtcls; i++)
  {
    const PosType* restrict dptr=GeminalBasis->dY[i];
    const RealType* restrict d2ptr=GeminalBasis->d2Y[i];
    const RealType* restrict vptr=V[0];
    BasisSetType::GradType grad(0.0);
    BasisSetType::ValueType lap(0.0);
    for(int j=0; j<NumPtcls; j++, vptr+=BasisSize)
    {
      if(j==i)
      {
        dUk(i,j) = 0.0;
        d2Uk(i,j)= 0.0;
      }
      else
      {
        grad+= (dUk(i,j) = simd::dot(vptr,dptr,BasisSize));
        lap += (d2Uk(i,j)= simd::dot(vptr,d2ptr,BasisSize));
      }
    }
    P.G[i]+=grad;
    P.L[i]+=lap;
  }
}

void
ThreeBodyGeminal::copyFromBuffer(ParticleSet& P, WFBufferType& buf)
{
  buf.get(LogValue);
  buf.get(V.begin(), V.end());
  buf.get(Y.begin(), Y.end());
  buf.get(FirstAddressOfdY,LastAddressOfdY);
  buf.get(d2Y.begin(),d2Y.end());
  buf.get(Uk.begin(), Uk.end());
  buf.get(FirstAddressOfgU,LastAddressOfgU);
  buf.get(d2Uk.begin(), d2Uk.end());
}

OrbitalBase::RealType
ThreeBodyGeminal::updateBuffer(ParticleSet& P, WFBufferType& buf,
                               bool fromscratch)
{
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

bool ThreeBodyGeminal::put(xmlNodePtr cur)
{
  //BasisSize = GeminalBasis->TotalBasis;
  BasisSize = GeminalBasis->getBasisSetSize();
  app_log() << "  The number of Geminal functions "
            <<"for Three-body Jastrow " << BasisSize << std::endl;
  app_log() << "  The number of particles " << NumPtcls << std::endl;
  Lambda.resize(BasisSize,BasisSize);
  //disable lambda's so that element-by-element input can be handled
  FreeLambda.resize(BasisSize*(BasisSize+1)/2);
  FreeLambda=false;
  //zero is default
  Lambda=0.0;
  //for(int ib=0; ib<BasisSize; ib++) Lambda(ib,ib)=NormFac;
  //for(int ib=0; ib<BasisSize; ib++)
  //  for(int jb=ib; jb<BasisSize; ++jb)
  //  {
  //    Lambda(ib,jb)=Random();
  //    if(jb!=ib) Lambda(jb,ib)=Lambda(ib,jb);
  //  }
  if(cur == NULL)
  {
    FreeLambda=true;
  }
  else
  {
    //read from an input nodes
    std::string aname("j3");
    std::string datatype("no");
    int sizeIn(0);
    IndexOffSet=1;
    OhmmsAttributeSet attrib;
    attrib.add(aname,"id");
    attrib.add(sizeIn,"size");
    attrib.add(aname,"name");
    attrib.add(datatype,"type");
    attrib.add(IndexOffSet,"offset");
    attrib.put(cur);
    ID_Lambda=aname;
    if(datatype.find("rray")<datatype.size())
    {
      if (sizeIn==Lambda.rows())
      {
        putContent(Lambda,cur);
      }
      FreeLambda=true;
      //addOptimizables(varlist);
      //symmetrize it
      //for(int ib=0; ib<BasisSize; ib++) {
      //  sprintf(coeffname,"%s_%d_%d",aname.c_str(),ib+IndexOffSet,ib+IndexOffSet);
      //  varlist[coeffname]=Lambda(ib,ib);
      //  for(int jb=ib+1; jb<BasisSize; jb++) {
      //    sprintf(coeffname,"%s_%d_%d",aname.c_str(),ib+IndexOffSet,jb+IndexOffSet);
      //    Lambda(jb,ib) = Lambda(ib,jb);
      //    varlist[coeffname]=Lambda(ib,jb);
      //  }
      //}
    }
    else
    {
      xmlNodePtr tcur=cur->xmlChildrenNode;
      while(tcur != NULL)
      {
        if(xmlStrEqual(tcur->name,(const xmlChar*)"lambda"))
        {
          int iIn=atoi((const char*)(xmlGetProp(tcur,(const xmlChar*)"i")));
          int jIn=atoi((const char*)(xmlGetProp(tcur,(const xmlChar*)"j")));
          int i=iIn-IndexOffSet;
          int j=jIn-IndexOffSet;
          double c=atof((const char*)(xmlGetProp(tcur,(const xmlChar*)"c")));
          Lambda(i,j)=c;
          FreeLambda[i*BasisSize+j]=true;
          if(i != j)
            Lambda(j,i)=c;
          //sprintf(coeffname,"%s_%d_%d",aname.c_str(),iIn,jIn);
          //varlist[coeffname]=c;
        }
        tcur=tcur->next;
      }
    }
  }
  //myVars are set
  myVars.clear();
  std::string coeffname;
  int ii=0;
  for(int ib=0; ib<BasisSize; ib++)
  {
    if(FreeLambda[ii++])
    {
      coeffname=ID_Lambda+"_"+std::to_string(ib)+"_"+std::to_string(ib);
      myVars.insert(coeffname,Lambda(ib,ib));
    }
    for(int jb=ib+1; jb<BasisSize; jb++)
    {
      if(FreeLambda[ii++])
      {
        coeffname=ID_Lambda+"_"+std::to_string(ib)+"_"+std::to_string(jb);
        myVars.insert(coeffname,Lambda(ib,jb));
      }
    }
  }
  //app_log() << "  Lambda Variables " << std::endl;
  //myVars.print(app_log());
  //app_log() << std::endl;
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
  //app_log() << "  Three-body Geminal coefficients " << std::endl;
  //app_log() << Lambda << std::endl;
  //GeminalBasis->resize(NumPtcls);
  return true;
}

OrbitalBasePtr ThreeBodyGeminal::makeClone(ParticleSet& tqp) const
{
  ThreeBodyGeminal* myclone=new ThreeBodyGeminal(CenterRef,tqp);
  myclone->GeminalBasis=GeminalBasis->makeClone();
  myclone->GeminalBasis->resetTargetParticleSet(tqp);
  myclone->put(NULL);
  myclone->myVars=myVars;
  myclone->Lambda=Lambda;
  myclone->ID_Lambda=ID_Lambda;
  myclone->FreeLambda=FreeLambda;
  return myclone;
}

//void ThreeBodyGeminal::addOptimizables(OptimizableSetType& varlist)
//{
//  char coeffname[16];
//  for(int ib=0; ib<BasisSize; ib++) {
//    sprintf(coeffname,"%s_%d_%d",ID_Lambda.c_str(),ib+IndexOffSet,ib+IndexOffSet);
//    varlist[coeffname]=Lambda(ib,ib);
//    for(int jb=ib+1; jb<BasisSize; jb++) {
//      sprintf(coeffname,"%s_%d_%d",ID_Lambda.c_str(),ib+IndexOffSet,jb+IndexOffSet);
//      Lambda(jb,ib) = Lambda(ib,jb);
//      varlist[coeffname]=Lambda(ib,jb);
//    }
//  }
//}
}


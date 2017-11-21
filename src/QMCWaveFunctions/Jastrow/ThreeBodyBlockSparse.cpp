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
#include "QMCWaveFunctions/Jastrow/ThreeBodyBlockSparse.h"
#include "Particle/DistanceTable.h"
#include "Particle/DistanceTableData.h"
#include "Numerics/DeterminantOperators.h"
#include "Numerics/OhmmsBlas.h"
#include "Numerics/BlockMatrixFunctions.h"
#include "Numerics/LibxmlNumericIO.h"
#include "Utilities/IteratorUtility.h"
#include "Numerics/MatrixOperators.h"
#include "OhmmsData/AttributeSet.h"

namespace qmcplusplus
{

ThreeBodyBlockSparse::ThreeBodyBlockSparse(const ParticleSet& ions, ParticleSet& els):
  CenterRef(ions), GeminalBasis(0), IndexOffSet(1), ID_Lambda("j3"), SameBlocksForGroup(true)
{
  myTableIndex=els.addTable(ions,DT_AOS);
  //d_table = DistanceTable::add(ions,els);
  NumPtcls=els.getTotalNum();
  Optimizable=true;
}

ThreeBodyBlockSparse::~ThreeBodyBlockSparse()
{
  delete_iter(LambdaBlocks.begin(), LambdaBlocks.end());
}

OrbitalBase::RealType
ThreeBodyBlockSparse::evaluateLog(ParticleSet& P,
                                  ParticleSet::ParticleGradient_t& G,
                                  ParticleSet::ParticleLaplacian_t& L)
{
  evaluateLogAndStore(P);
  //GeminalBasis->evaluateForWalkerMove(P);
  ////this could be better but it is used only sparsely
  ////MatrixOperators::product(GeminalBasis->Y, Lambda, V);
  //MatrixOperators::product(GeminalBasis->Y, Lambda, V, BlockOffset);
  //Uk=0.0;
  //LogValue=ValueType();
  //for(int i=0; i< NumPtcls-1; i++) {
  //  const BasisSetType::RealType* restrict yptr=GeminalBasis->Y[i];
  //  for(int j=i+1; j<NumPtcls; j++) {
  //    RealType x= dot(V[j],yptr,BasisSize);
  //    LogValue += x;
  //    Uk[i]+= x;
  //    Uk[j]+= x;
  //  }
  //}
  //for(int i=0; i<NumPtcls; i++)  {
  //  const BasisSetType::GradType* restrict dptr=GeminalBasis->dY[i];
  //  const BasisSetType::ValueType* restrict d2ptr=GeminalBasis->d2Y[i];
  //  const BasisSetType::ValueType* restrict vptr=V[0];
  //  BasisSetType::GradType grad(0.0);
  //  BasisSetType::ValueType lap(0.0);
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
ThreeBodyBlockSparse::ratio(ParticleSet& P, int iat)
{
  UpdateMode=ORB_PBYP_RATIO;
  GeminalBasis->evaluateForPtclMove(P,iat);
  const BasisSetType::ValueType* restrict y_ptr=GeminalBasis->Phi.data();
  //This is the only difference from ThreeBodyGeminal
  RealType* restrict cv_ptr = curV.data();
  for(int b=0; b<Blocks.size(); b++)
  {
    int nb=Blocks[b];
    if(nb)
    {
      //GEMV<RealType,0>::apply(LambdaBlocks[b]->data(),y_ptr,cv_ptr,nb,nb);
      GEMV<RealType,0>::apply(LambdaBlocks[BlockID[b]]->data(),y_ptr,cv_ptr,nb,nb);
      for(int ib=0,k=BlockOffset[b]; ib<nb; k++,ib++)
        delV[k] = (*cv_ptr++)-V[iat][k];
      y_ptr+=nb;
    }
  }
  diffVal=0.0;
  const RealType* restrict vptr=V[0];
  for(int j=0; j<NumPtcls; j++, vptr+=BasisSize)
  {
    if(j == iat)
      continue;
    diffVal+= (curVal[j]=simd::dot(delV.data(),Y[j],BasisSize));
  }
  curVal[iat]=diffVal;
  return std::exp(diffVal);
}

void ThreeBodyBlockSparse::restore(int iat)
{
  //nothing to do here
}

void ThreeBodyBlockSparse::acceptMove(ParticleSet& P, int iat)
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
ThreeBodyBlockSparse::registerData(ParticleSet& P, WFBufferType& buf)
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
ThreeBodyBlockSparse::evaluateLogAndStore(ParticleSet& P)
{
  GeminalBasis->evaluateForWalkerMove(P);
  //this could be better but it is used only sparsely
  //MatrixOperators::product(GeminalBasis->Y, Lambda, V);
  MatrixOperators::product(GeminalBasis->Y, Lambda, V, BlockOffset);
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
    const BasisSetType::GradType* restrict dptr=GeminalBasis->dY[i];
    const BasisSetType::ValueType* restrict d2ptr=GeminalBasis->d2Y[i];
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
ThreeBodyBlockSparse::copyFromBuffer(ParticleSet& P, WFBufferType& buf)
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
ThreeBodyBlockSparse::updateBuffer(ParticleSet& P, WFBufferType& buf,
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

void ThreeBodyBlockSparse::resetTargetParticleSet(ParticleSet& P)
{
  //d_table = DistanceTable::add(CenterRef,P);
  GeminalBasis->resetTargetParticleSet(P);
}

void ThreeBodyBlockSparse::checkInVariables(opt_variables_type& active)
{
  active.insertFrom(myVars);
  int ncur=active.size();
  GeminalBasis->checkInVariables(active);
  if(active.size()>ncur)
    GeminalBasis->checkInVariables(myVars);
}

void ThreeBodyBlockSparse::checkOutVariables(const opt_variables_type& active)
{
  myVars.getIndex(active);
  GeminalBasis->checkOutVariables(active);
  Optimizable=myVars.is_optimizable();
  app_log() << "<j3-variables>"<< std::endl;
  myVars.print(app_log());
  app_log() << "</j3-variables>"<< std::endl;
}

///reset the value of all the Two-Body Jastrow functions
void ThreeBodyBlockSparse::resetParameters(const opt_variables_type& active)
{
  int ii=0;
  for(int b=0; b<LambdaBlocks.size(); ++b)
  {
    Matrix<RealType>& m(*LambdaBlocks[b]);
    for(int ib=0; ib<m.rows(); ++ib)
    {
      int loc=myVars.where(ii++);
      if(loc>=0)
        m(ib,ib)=active[loc];
    }
    for(int ib=0; ib<m.rows(); ++ib)
      for(int jb=ib+1; jb<m.rows(); ++jb)
      {
        int loc=myVars.where(ii++);
        if(loc>=0)
          m(jb,ib)=m(ib,jb)=active[loc];
      }
  }
  //this is not  necessary
  if(SameBlocksForGroup)
    checkLambda();
  GeminalBasis->resetParameters(active);
  for(int i=0; i<myVars.size(); ++i)
    if(myVars.where(i)>=0)
      myVars[i]=active[myVars.where(i)];
  //app_log() << "ThreeBodyBlockSparse::resetParameters " << std::endl << Lambda << std::endl;
  //for(int b=0; b<LambdaBlocks.size();++b)
  //{
  //  app_log() << "3Blody Block " << b << std::endl << *LambdaBlocks[b] << std::endl;
  //}
}

void ThreeBodyBlockSparse::reportStatus(std::ostream& os)
{
  myVars.print(os);
}

bool ThreeBodyBlockSparse::put(xmlNodePtr cur)
{
  //BasisSize = GeminalBasis->TotalBasis;
  BasisSize = GeminalBasis->getBasisSetSize();
  app_log() << "  The number of Geminal functions "
            <<"for Three-body Jastrow " << BasisSize << std::endl;
  app_log() << "  The number of particles " << NumPtcls << std::endl;
  Lambda.resize(BasisSize,BasisSize);
  FreeLambda.resize(BasisSize,BasisSize);
  FreeLambda=false;
  //identity is the default
  for(int ib=0; ib<BasisSize; ib++)
    Lambda(ib,ib)=1.0;
  if(cur == NULL)
  {
    FreeLambda=true;
  }
  else
  {
    //read from an input nodes
    char coeffname[16];
    std::string aname("j3");
    std::string datatype("lambda");
    std::string sameblocks("yes");
    IndexOffSet=0;//defaults are c
    OhmmsAttributeSet attrib;
    attrib.add(aname,"id");
    attrib.add(aname,"name");
    attrib.add(datatype,"type");
    attrib.add(IndexOffSet,"offset");
    attrib.add(sameblocks,"sameBlocksForGroup");
    attrib.put(cur);
    SameBlocksForGroup = (sameblocks == "yes");
    ID_Lambda=aname;
    if(datatype.find("rray")<datatype.size())
    {
      putContent(Lambda,cur);
      FreeLambda=true;
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
          FreeLambda(i,j)=true;
          if(i != j)
            Lambda(j,i)=c;
          //sprintf(coeffname,"%s_%d_%d",aname.c_str(),iIn,jIn);
          //varlist[coeffname]=c;
        }
        tcur=tcur->next;
      }
    }
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
  return true;
}


void ThreeBodyBlockSparse::setBlocks(const std::vector<int>& blockspergroup)
{
  //copyc this
  BlocksPerGroup=blockspergroup;
  //test with three blocks
  Blocks.resize(CenterRef.getTotalNum());
  for(int i=0; i<CenterRef.getTotalNum(); i++)
    Blocks[i]=blockspergroup[CenterRef.GroupID[i]];
  BlockOffset.resize(Blocks.size()+1,0);
  for(int i=0; i<Blocks.size(); i++)
    BlockOffset[i+1]=BlockOffset[i]+Blocks[i];
  BlockID.resize(Blocks.size(),-1);
  if(LambdaBlocks.empty())
  {
    std::vector<int> firstK,lastK;
    if(SameBlocksForGroup)
    {
      std::vector<int> need2copy(CenterRef.getSpeciesSet().getTotalNum(),-1);
      firstK.resize(need2copy.size(),-1);
      lastK.resize(need2copy.size(),-1);
      LambdaBlocks.resize(need2copy.size(),0);
      for(int b=0; b<Blocks.size(); b++)
      {
        if(Blocks[b] == 0)
          continue;
        int gid=need2copy[CenterRef.GroupID[b]];
        if(gid<0)//assign the current block index
        {
          gid=CenterRef.GroupID[b];
          need2copy[gid] = gid;
          LambdaBlocks[gid]= new Matrix<RealType>(Blocks[b],Blocks[b]);
          //populate LambdaBlocks with input Lambda
          Matrix<RealType>& m(*LambdaBlocks[gid]);
          firstK[gid]=BlockOffset[b];
          lastK[gid]=BlockOffset[b+1];
          for(int k=BlockOffset[b],ib=0; k<BlockOffset[b+1]; ++k,++ib)
            for(int kp=BlockOffset[b],jb=0; kp<BlockOffset[b+1]; ++kp,++jb)
              m(ib,jb)=Lambda(k,kp);
        }
        BlockID[b]=gid;
      }
    }
    else
    {
      LambdaBlocks.resize(Blocks.size(),0);
      firstK.resize(Blocks.size(),-1);
      lastK.resize(Blocks.size(),-1);
      for(int b=0; b<Blocks.size(); b++)
      {
        if(Blocks[b] ==0)
          continue;
        BlockID[b]=b;
        LambdaBlocks[b]=new Matrix<RealType>(Blocks[b],Blocks[b]);
        Matrix<RealType>& m(*LambdaBlocks[b]);
        firstK[b]=BlockOffset[b];
        lastK[b]=BlockOffset[b+1];
        for(int k=BlockOffset[b],ib=0; k<BlockOffset[b+1]; ++k,++ib)
          for(int kp=BlockOffset[b],jb=0; kp<BlockOffset[b+1]; ++kp,++jb)
            m(ib,jb)=Lambda(k,kp);
      }
    }
    myVars.clear();
    for(int b=0; b<LambdaBlocks.size(); ++b)
    {
      int first=firstK[b];
      int last=lastK[b];
      Matrix<RealType>& m(*LambdaBlocks[b]);
      for(int k=first,ib=0; k<last; ++k,++ib)
      {
        std::ostringstream sstr;
        sstr << ID_Lambda<<"_"<<k<< "_"<<k;
        myVars.insert(sstr.str(),m(ib,ib),true);
      }
      for(int k=first,ib=0; k<last; ++k,++ib)
        for(int kp=k+1,jb=ib+1; kp<last; ++kp,++jb)
        {
          std::ostringstream sstr;
          sstr << ID_Lambda<<"_"<<k<< "_"<<kp;
          myVars.insert(sstr.str(),m(ib,jb),true);
        }
    }
  }
  if(SameBlocksForGroup)
    checkLambda();
  app_log() << "  ThreeBodyBlockSparse handles optimizable variables " << std::endl;
  myVars.print(app_log());
}

/** set dependent Lambda
 */
void ThreeBodyBlockSparse::checkLambda()
{
  for(int b=0; b<Blocks.size(); b++)
  {
    if(Blocks[b] ==0)
      continue;
    Matrix<RealType>& m(*LambdaBlocks[BlockID[b]]);
    for(int k=BlockOffset[b],ib=0; k<BlockOffset[b+1]; k++,ib++)
      for(int kp=BlockOffset[b],jb=0; kp<BlockOffset[b+1]; kp++,jb++)
        Lambda(k,kp)=m(ib,jb);
  }
}

OrbitalBasePtr ThreeBodyBlockSparse::makeClone(ParticleSet& tqp) const
{
  ThreeBodyBlockSparse* myclone=new ThreeBodyBlockSparse(CenterRef,tqp);
  myclone->GeminalBasis=GeminalBasis->makeClone();
  myclone->GeminalBasis->resetTargetParticleSet(tqp);
  myclone->put(NULL);//to resize data
  myclone->Lambda=Lambda;
  myclone->FreeLambda=FreeLambda;
  myclone->ID_Lambda=ID_Lambda;
  myclone->setBlocks(BlocksPerGroup);
  return myclone;
}

//void ThreeBodyBlockSparse::addOptimizables(OptimizableSetType& varlist)
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


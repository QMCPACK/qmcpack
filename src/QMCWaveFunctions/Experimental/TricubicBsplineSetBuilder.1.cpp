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
    
    
#include "QMCWaveFunctions/PlaneWave/PWParameterSet.h"
#include "QMCWaveFunctions/TricubicBsplineSetBuilder.h"
#include "OhmmsData/AttributeSet.h"
#include "Numerics/OhmmsBlas.h"
#include "Message/OpenMP.h"
#include "Message/CommOperators.h"
#include "Utilities/ProgressReportEngine.h"
namespace qmcplusplus
{

template<typename T1, typename T2>
struct bspline_data_check
{
  enum {is_same=0};
};

template<typename T1>
struct bspline_data_check<T1,T1>
{
  enum {is_same=1};
};

template<typename Tin, typename Tout>
void
TricubicBsplineSetBuilder::readData(const char* hroot, const std::vector<int>& occSet, int spinIndex, int degeneracy)
{
  ReportEngine PRE(ClassName,"readData");
  bool truncate = (myParam->Rcut>0.0);
  int norb=occSet.size()/degeneracy;
  typedef Array<Tin,DIM> InArrayType;
  InArrayType inTemp(BoxGrid[0],BoxGrid[1],BoxGrid[2]);
  StorageType inData(BoxGrid[0],BoxGrid[1],BoxGrid[2]);
  getCenterAndOrigin(hroot,occSet,norb);
  if(bspline_data_check<Tin,Tout>::is_same)
    PRE<<"  No data conversion is needed.\n";
  else
    PRE<< "  Data conversion is done.\n";
  for(int iorb=0; iorb<norb; iorb++)
  {
    std::ostringstream wnshort;
    std::string eigvName=myParam->getEigVectorName(hroot,occSet[iorb]/degeneracy,spinIndex);
    wnshort<<curH5Fname << "#"<<occSet[iorb]/degeneracy << "#" << spinIndex;
    std::map<std::string,RSOType*>::iterator it(BigDataSet.find(wnshort.str()));
    if(it == BigDataSet.end())
    {
      if(ReportLevel)
      {
        PRE << "   Reading spline function " << eigvName << " (" << wnshort.str()  << ")";
        if(truncate)
          PRE<< "\n     center=" << activeBasis->Centers[iorb];
        PRE<<"\n";
      }
      StorageType* newP=new StorageType;
      BigDataSet[wnshort.str()]=new RSOType(activeBasis->Centers[iorb],activeBasis->Origins[iorb],newP);
      if(bspline_data_check<Tin,Tout>::is_same)
      {
        if(is_manager())
        {
          HDFAttribIO<StorageType> dummy(inData);
          dummy.read(hfileID,eigvName.c_str());
        }
        myComm->bcast(inData);
      }
      else
      {
        if(is_manager())
        {
          HDFAttribIO<InArrayType> dummy(inTemp);
          dummy.read(hfileID,eigvName.c_str());
        }
        myComm->bcast(inTemp);
        BLAS::copy(inTemp.size(),inTemp.data(),inData.data());
      }
      activeBasis->add(iorb,inData,newP);
    }
    else
    {
      if(ReportLevel)
      {
        PRE << "   Reusing spline function " << eigvName << " (" << wnshort.str()  << ")";
        if(truncate)
          PRE<< "\n     center=" << activeBasis->Centers[iorb];
        PRE<<"\n";
      }
      activeBasis->add(iorb,(*it).second->Coeffs);
    }
  }
}

template<typename Tin, typename Tout>
void TricubicBsplineSetBuilder::readDataOMP(const char* hroot, const std::vector<int>& occSet, int spinIndex, int degeneracy)
{
  ReportEngine PRE(ClassName,"readDataOMP");
  if(BigDataSet.size()) //data exist, use the standard one
  {
    readData<Tin,Tout>(hroot,occSet,spinIndex,degeneracy);
    return;
  }
  int norb=occSet.size()/degeneracy;
  std::vector<StorageType*> bsset(norb);
  std::vector<int> odist(omp_get_max_threads()+1);
  FairDivideLow(norb,omp_get_max_threads(),odist);
  #pragma omp parallel
  {
    typedef Array<Tin,DIM> InArrayType;
    InArrayType inTemp(BoxGrid[0],BoxGrid[1],BoxGrid[2]);
    StorageType inData(BoxGrid[0],BoxGrid[1],BoxGrid[2]);
    int ip=omp_get_thread_num();
    for(int iorb=odist[ip]; iorb<odist[ip+1]; iorb++)
    {
      #pragma omp critical
      {
        std::string eigvName=myParam->getEigVectorName(hroot,occSet[iorb]/degeneracy,spinIndex);
        PRE << "   Reading spline function " << ip << " "  << eigvName <<'\n';
        if(bspline_data_check<Tin,Tout>::is_same)
        {
          HDFAttribIO<StorageType> dummy(inData);
          dummy.read(hfileID,eigvName.c_str());
        }
        else
        {
          HDFAttribIO<InArrayType> dummy(inTemp);
          dummy.read(hfileID,eigvName.c_str());
          BLAS::copy(inTemp.size(),inTemp.data(),inData.data());
        }
      }
      StorageType* newP =new StorageType;
      bsset[iorb]=newP;
      activeBasis->add(iorb,inData,newP);
    }
  }
  getCenterAndOrigin(hroot,occSet,norb);
  for(int iorb=0; iorb<norb; iorb++)
  {
    std::ostringstream wnshort;
    wnshort<<curH5Fname << "#"<<occSet[iorb]/degeneracy << "#" << spinIndex;
    BigDataSet[wnshort.str()]=new RSOType(activeBasis->Centers[iorb],activeBasis->Origins[iorb],bsset[iorb]);
  }
}

void TricubicBsplineSetBuilder::setBsplineBasisSet(xmlNodePtr cur)
{
  ReportEngine PRE(ClassName,"setBsplineBasisSet");
  int norb(0);
  int degeneracy(1);
  std::string ompIn("no");
  OhmmsAttributeSet aAttrib;
  aAttrib.add(norb,"orbitals");
  aAttrib.add(norb,"size");
  aAttrib.add(ompIn,"omp");
  //aAttrib.add(degeneracy,"degeneracy");
  aAttrib.put(cur);
  std::vector<int> occSet(norb);
  for(int i=0; i<norb; i++)
    occSet[i]=i;
  //set the root name, e.g., /eigenstates/twist0
  char hroot[128];
  sprintf(hroot,"/%s/%s%d",
          myParam->eigTag.c_str(), myParam->twistTag.c_str(), myParam->twistIndex);
  int spinIndex=0;
  cur=cur->children;
  while(cur != NULL)
  {
    std::string cname((const char*)(cur->name));
    if(cname == "occupation")
    {
      std::string occ_mode("ground");
      OhmmsAttributeSet oAttrib;
      oAttrib.add(occ_mode,"mode");
      oAttrib.add(spinIndex,"spindataset");
      oAttrib.put(cur);
      //Do nothing if mode == ground
      if(occ_mode == "excited")
      {
        std::vector<int> occ_in, occRemoved;
        putContent(occ_in,cur);
        for(int k=0; k<occ_in.size(); k++)
        {
          if(occ_in[k]<0)
            occRemoved.push_back(-occ_in[k]-1);
        }
        int kpopd=0;
        for(int k=0; k<occ_in.size(); k++)
        {
          if(occ_in[k]>0)
            occSet[occRemoved[kpopd++]]=occ_in[k]-1;
        }
      }
      const xmlChar* h5path = xmlGetProp(cur,(const xmlChar*)"h5path");
      if(h5path != NULL)
        sprintf(hroot,"%s",(const char*)h5path);
    }
    cur=cur->next;
  }
  //calculate degeneracy with BoxDup
  degeneracy=BoxDup[0]*BoxDup[1]*BoxDup[2];
  //always translate grid
  if(degeneracy>1)
    TranslateGrid=true;
  bool in_complex = myParam->getEigVectorType(hfileID);
  //disable readDataOMP since mpi cannot handle communications within threads
  //if(ompIn == "no")
  //{
  if(in_complex)
    readData<ComplexType,ValueType>(hroot,occSet,spinIndex,degeneracy);
  else
    readData<RealType,ValueType>(hroot,occSet,spinIndex,degeneracy);
  //}
  //else
  //{
  //  if(in_complex)
  //    readDataOMP<ComplexType,ValueType>(hroot,occSet,spinIndex,degeneracy);
  //  else
  //    readDataOMP<RealType,ValueType>(hroot,occSet,spinIndex,degeneracy);
  //}
  if(degeneracy>1)
  {
    if(myParam->Rcut>0.0)//duplicate only the localized orbitals
    {
      PRE<<"Tiling localized orbitals by " << BoxDup << "\n";
      activeBasis->tileOrbitals(BoxDup);
    }
    else
    {
      //fatal
      PRE.error("Tiling is not allowed with non-localized orbitals.",true);
    }
  }
}

void TricubicBsplineSetBuilder::getCenterAndOrigin(const char* hroot, const std::vector<int>& occSet, int norb)
{
  if(myParam->Rcut>0.0)
  {
    if(myComm->rank()==0)
    {
      for(int iorb=0; iorb<norb; iorb++)
      {
        std::string centerName=myParam->getCenterName(hroot,occSet[iorb]);
        HDFAttribIO<PosType > cdummy(activeBasis->Centers[iorb]);
        cdummy.read(hfileID,centerName.c_str());
      }
    }
    myComm->bcast(activeBasis->Centers);
  }
  if(FloatingGrid)
  {
    if(myComm->rank() ==0)
    {
      for(int iorb=0; iorb<norb; iorb++)
      {
        std::string originName=myParam->getOriginName(hroot,occSet[iorb]);
        HDFAttribIO<PosType > cdummy(activeBasis->Origins[iorb]);
        cdummy.read(hfileID,originName.c_str());
      }
    }
    myComm->bcast(activeBasis->Origins);
  }
}

//  void TricubicBsplineSetBuilder::readComplex2RealData(const char* hroot, const std::vector<int>& occSet,
//      int spinIndex, int degeneracy)
//  {
//    bool truncate = (myParam->Rcut>0.0);
//    //PosType center(0.0),origin(0.0);
//    int norb=occSet.size()/degeneracy;
//    StorageType inData(BoxGrid[0],BoxGrid[1],BoxGrid[2]);
//    Array<ComplexType,3> inTemp(BoxGrid[0],BoxGrid[1],BoxGrid[2]);
//
//    //first obtain centers and origins
//    getCenterAndOrigin(hroot,occSet,norb);
//
//    for(int iorb=0; iorb<norb; iorb++)
//    {
//      std::ostringstream wnshort;
//      std::string eigvName=myParam->getEigVectorName(hroot,occSet[iorb],spinIndex);
//      wnshort<<curH5Fname << "#"<<occSet[iorb] << "#" << spinIndex;
//      std::map<std::string,RSOType*>::iterator it(BigDataSet.find(wnshort.str()));
//      if(it == BigDataSet.end()) {
//        if(print_log)
//        {
//          app_log() << "   Reading spline function " << eigvName << " (" << wnshort.str()  << ")"  << std::endl;
//          if(truncate)
//            app_log() << "     center=" << activeBasis->Centers[iorb] << std::endl;
//        }
//        StorageType* newP =new StorageType;
//        HDFAttribIO<Array<ComplexType,3> > dummy(inTemp);
//        dummy.read(hfileID,eigvName.c_str());
//        BLAS::copy(inTemp.size(),inTemp.data(),inData.data());
//        BigDataSet[wnshort.str()]=new RSOType(activeBasis->Centers[iorb],activeBasis->Origins[iorb],newP);
//        activeBasis->add(iorb,inData,newP);
//      } else {
//        if(print_log)
//        {
//          app_log() << "   Reusing spline function " << eigvName << " (" << wnshort.str()  << ")" << std::endl;
//          if(truncate)
//            app_log() << "     center=" << activeBasis->Centers[iorb] << std::endl;
//        }
//        activeBasis->add(iorb,(*it).second->Coeffs);
//      }
//    }
//
//  }
//
//  void TricubicBsplineSetBuilder::readComplex2RealDataOMP(const char* hroot, const std::vector<int>& occSet,
//      int spinIndex, int degeneracy)
//  {
//
//    if(BigDataSet.size()) //data exist, use the standard one
//    {
//      readComplex2RealData(hroot,occSet,spinIndex,degeneracy);
//      return;
//    }
//
//    int norb=occSet.size()/degeneracy;
//    std::vector<StorageType*> bsset(norb);
//    std::vector<int> odist(omp_get_max_threads()+1);
//    FairDivideLow(norb,omp_get_max_threads(),odist);
//
//#pragma omp parallel
//    {
//      bool truncate = (myParam->Rcut>0.0);
//      StorageType inData(BoxGrid[0],BoxGrid[1],BoxGrid[2]);
//      Array<ComplexType,3> inTemp(BoxGrid[0],BoxGrid[1],BoxGrid[2]);
//
//      int ip=omp_get_thread_num();
//      for(int iorb=odist[ip]; iorb<odist[ip+1]; iorb++)
//      {
//#pragma omp critical
//        {
//          std::string eigvName=myParam->getEigVectorName(hroot,occSet[iorb]/degeneracy,spinIndex);
//          app_log() << "   Reading spline function " << ip << " "  << eigvName << std::endl;
//          HDFAttribIO<Array<ComplexType,3> > dummy(inTemp);
//          dummy.read(hfileID,eigvName.c_str());
//        }
//        BLAS::copy(inTemp.size(),inTemp.data(),inData.data());
//        StorageType* newP =new StorageType;
//        bsset[iorb]=newP;
//        activeBasis->add(iorb,inData,newP);
//      }
//    }
//
//    //first obtain centers and origins
//    getCenterAndOrigin(hroot,occSet,norb);
//    for(int iorb=0; iorb<norb; iorb++)
//    {
//      std::ostringstream wnshort;
//      wnshort<<curH5Fname << "#"<<occSet[iorb]/degeneracy << "#" << spinIndex;
//      BigDataSet[wnshort.str()]=new RSOType(activeBasis->Centers[iorb],activeBasis->Origins[iorb],bsset[iorb]);
//    }
//  }
//
//  void TricubicBsplineSetBuilder::readData(const char* hroot, const std::vector<int>& occSet,
//      int spinIndex, int degeneracy)
//  {
//    bool truncate = (myParam->Rcut>0.0);
//    int norb=occSet.size()/degeneracy;
//    StorageType inData(BoxGrid[0],BoxGrid[1],BoxGrid[2]);
//
//    getCenterAndOrigin(hroot,occSet,norb);
//    for(int iorb=0; iorb<norb; iorb++)
//    {
//      std::ostringstream wnshort;
//      std::string eigvName=myParam->getEigVectorName(hroot,occSet[iorb]/degeneracy,spinIndex);
//      wnshort<<curH5Fname << "#"<<occSet[iorb]/degeneracy << "#" << spinIndex;
//      std::map<std::string,RSOType*>::iterator it(BigDataSet.find(wnshort.str()));
//      if(it == BigDataSet.end())
//      {
//        if(print_log)
//        {
//          app_log() << "   Reading spline function " << eigvName << " (" << wnshort.str()  << ")"  << std::endl;
//          if(truncate)
//            app_log() << "     center=" << activeBasis->Centers[iorb] << std::endl;
//        }
//        StorageType* newP=new StorageType;
//        HDFAttribIO<StorageType> dummy(inData);
//        dummy.read(hfileID,eigvName.c_str());
//        BigDataSet[wnshort.str()]=new RSOType(activeBasis->Centers[iorb],activeBasis->Origins[iorb],newP);
//        activeBasis->add(iorb,inData,newP);
//      }
//      else
//      {
//        if(print_log)
//        {
//          app_log() << "   Reusing spline function " << eigvName << " (" << wnshort.str()  << ")" << std::endl;
//          if(truncate)
//            app_log() << "     center=" << activeBasis->Centers[iorb] << std::endl;
//        }
//        activeBasis->add(iorb,(*it).second->Coeffs);
//      }
//    }
//  }
//
//  void TricubicBsplineSetBuilder::readDataOMP(const char* hroot, const std::vector<int>& occSet,
//      int spinIndex, int degeneracy)
//  {
//
//    if(BigDataSet.size()) //data exist, use the standard one
//    {
//      readData(hroot,occSet,spinIndex,degeneracy);
//      return;
//    }
//
//    int norb=occSet.size()/degeneracy;
//    std::vector<StorageType*> bsset(norb);
//    std::vector<int> odist(omp_get_max_threads()+1);
//    FairDivideLow(norb,omp_get_max_threads(),odist);
//
//#pragma omp parallel
//    {
//      bool truncate = (myParam->Rcut>0.0);
//      StorageType inData(BoxGrid[0],BoxGrid[1],BoxGrid[2]);
//      int ip=omp_get_thread_num();
//      for(int iorb=odist[ip]; iorb<odist[ip+1]; iorb++)
//      {
//#pragma omp critical
//        {
//          std::string eigvName=myParam->getEigVectorName(hroot,occSet[iorb]/degeneracy,spinIndex);
//          app_log() << "   Reading spline function " << ip << " "  << eigvName << std::endl;
//          HDFAttribIO<StorageType> dummy(inData);
//          dummy.read(hfileID,eigvName.c_str());
//        }
//        StorageType* newP =new StorageType;
//        bsset[iorb]=newP;
//        activeBasis->add(iorb,inData,newP);
//      }
//    }
//
//    getCenterAndOrigin(hroot,occSet,norb);
//    for(int iorb=0; iorb<norb; iorb++)
//    {
//      std::ostringstream wnshort;
//      wnshort<<curH5Fname << "#"<<occSet[iorb]/degeneracy << "#" << spinIndex;
//      BigDataSet[wnshort.str()]=new RSOType(activeBasis->Centers[iorb],activeBasis->Origins[iorb],bsset[iorb]);
//    }
//  }
}

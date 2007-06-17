//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
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
#include "QMCWaveFunctions/PlaneWave/PWParameterSet.h"
#include "QMCWaveFunctions/TricubicBsplineSetBuilder.h"
#include "OhmmsData/AttributeSet.h"
#include "Numerics/OhmmsBlas.h"
#include "Message/Communicate.h"
#include "Message/OpenMP.h"
namespace qmcplusplus {

  void TricubicBsplineSetBuilder::setBsplineBasisSet(xmlNodePtr cur)
  {
    int norb(0);
    int degeneracy(1);
    string ompIn("no");
    OhmmsAttributeSet aAttrib;
    aAttrib.add(norb,"orbitals"); aAttrib.add(norb,"size");
    aAttrib.add(ompIn,"omp"); 
    //aAttrib.add(degeneracy,"degeneracy");
    aAttrib.put(cur);

    vector<int> occSet(norb);
    for(int i=0; i<norb; i++) occSet[i]=i;

    //set the root name, e.g., /eigenstates/twist0
    char hroot[128];
    sprintf(hroot,"/%s/%s%d",
        myParam->eigTag.c_str(), myParam->twistTag.c_str(), myParam->twistIndex);

    int spinIndex=0;
    cur=cur->children;
    while(cur != NULL) {
      string cname((const char*)(cur->name));
      if(cname == "occupation") {
        string occ_mode("ground");
        OhmmsAttributeSet oAttrib;
        oAttrib.add(occ_mode,"mode");
        oAttrib.add(spinIndex,"spindataset");
        oAttrib.put(cur);
        //Do nothing if mode == ground
        if(occ_mode == "excited") {
          vector<int> occ_in, occRemoved;
          putContent(occ_in,cur);
          for(int k=0; k<occ_in.size(); k++) {
            if(occ_in[k]<0) 
              occRemoved.push_back(-occ_in[k]-1);
          }
          int kpopd=0;
          for(int k=0; k<occ_in.size(); k++) {
            if(occ_in[k]>0) 
              occSet[occRemoved[kpopd++]]=occ_in[k]-1;
          }
        }
        const xmlChar* h5path = xmlGetProp(cur,(const xmlChar*)"h5path");
        if(h5path != NULL) sprintf(hroot,"%s",(const char*)h5path);
      }
      cur=cur->next;
    }

    //calculate degeneracy with BoxDup
    degeneracy=BoxDup[0]*BoxDup[1]*BoxDup[2];

    bool complex2real = myParam->hasComplexData(hfileID);
#if defined(QMC_COMPLEX)
    if(!complex2real) 
    {
      app_error() << "  Real wavefunctions cannot be used with QMC_COMPLEX=1" << endl;
      abort(); //FIXABORT
    }
#endif

    if(complex2real)
    {
      if(ompIn == "yes")
        readComplex2RealDataOMP(hroot,occSet,spinIndex,degeneracy);
      else
        readComplex2RealData(hroot,occSet,spinIndex,degeneracy);
    }
    else
      readData(hroot,occSet,spinIndex,degeneracy);

    if(degeneracy>1)
    { 
      if(myParam->Rcut>0.0)//duplicate only the localized orbitals
      {
        if(print_log)
          app_log() << "  Tiling localized orbitals by " << BoxDup << endl;
        activeBasis->tileOrbitals(BoxDup);
      }
      else
      {
        app_error() << "  Tiling is not allowed with non-localized orbitals." << endl;
        abort(); //FIXABORT
      }
    }
  }

  void TricubicBsplineSetBuilder::readComplex2RealData(const char* hroot, const vector<int>& occSet,
      int spinIndex, int degeneracy)
  {
    bool truncate = (myParam->Rcut>0.0);
    PosType center(0.0);
    int norb=occSet.size()/degeneracy;
    StorageType inData(BoxGrid[0],BoxGrid[1],BoxGrid[2]);
    Array<ComplexType,3> inTemp(BoxGrid[0],BoxGrid[1],BoxGrid[2]);
    for(int iorb=0; iorb<norb; iorb++) 
    {
      if(truncate)
      {
        string centerName=myParam->getCenterName(hroot,occSet[iorb]);
        HDFAttribIO<PosType > cdummy(center);
        cdummy.read(hfileID,centerName.c_str());
      }
      ostringstream wnshort;
      string eigvName=myParam->getEigVectorName(hroot,occSet[iorb],spinIndex);
      wnshort<<curH5Fname << "#"<<occSet[iorb] << "#" << spinIndex;
      map<string,StorageType*>::iterator it(BigDataSet.find(wnshort.str()));
      if(it == BigDataSet.end()) {
        if(print_log)
        {
          app_log() << "   Reading spline function " << eigvName << " (" << wnshort.str()  << ")"  << endl;
          if(truncate) app_log() << "     center=" << center << endl;
        }
        StorageType* newP =new StorageType;
        HDFAttribIO<Array<ComplexType,3> > dummy(inTemp);
        dummy.read(hfileID,eigvName.c_str());
        BLAS::copy(inTemp.size(),inTemp.data(),inData.data());
        BigDataSet[wnshort.str()]=newP;
        activeBasis->add(iorb,center,inData,newP);
      } else {
        if(print_log)
        {
          app_log() << "   Reusing spline function " << eigvName << " (" << wnshort.str()  << ")" << endl;
          if(truncate) app_log() << "     center=" << center << endl;
        }
        activeBasis->add(iorb,center,(*it).second);
      }
    } 

  }

  void TricubicBsplineSetBuilder::readComplex2RealDataOMP(const char* hroot, const vector<int>& occSet,
      int spinIndex, int degeneracy)
  {

    if(BigDataSet.size()) //data exist, use the standard one
    {
      readComplex2RealData(hroot,occSet,spinIndex,degeneracy);
      return;
    }

    int norb=occSet.size()/degeneracy;
    vector<StorageType*> bsset(norb);
    vector<int> odist(omp_get_max_threads()+1);
    FairDivideLow(norb,omp_get_max_threads(),odist);

#pragma omp parallel
    {
      bool truncate = (myParam->Rcut>0.0);
      PosType center;
      StorageType inData(BoxGrid[0],BoxGrid[1],BoxGrid[2]);
      Array<ComplexType,3> inTemp(BoxGrid[0],BoxGrid[1],BoxGrid[2]);

      int ip=omp_get_thread_num();
      for(int iorb=odist[ip]; iorb<odist[ip+1]; iorb++)
      {
#pragma omp critical
        {
          string eigvName=myParam->getEigVectorName(hroot,occSet[iorb]/degeneracy,spinIndex);
          app_log() << "   Reading spline function " << ip << " "  << eigvName << endl;
          HDFAttribIO<Array<ComplexType,3> > dummy(inTemp);
          dummy.read(hfileID,eigvName.c_str());
        }
        BLAS::copy(inTemp.size(),inTemp.data(),inData.data());
        StorageType* newP =new StorageType;
        bsset[iorb]=newP;
        activeBasis->add(iorb,center,inData,newP);
      }
    }

    bool localize=myParam->Rcut>0.0;
    for(int iorb=0; iorb<norb; iorb++) 
    {
      if(localize)
      {
        string centerName=myParam->getCenterName(hroot,occSet[iorb]);
        HDFAttribIO<PosType > cdummy(activeBasis->Centers[iorb]);
        cdummy.read(hfileID,centerName.c_str());
      }
      ostringstream wnshort;
      wnshort<<curH5Fname << "#"<<occSet[iorb]/degeneracy << "#" << spinIndex;
      BigDataSet[wnshort.str()]=bsset[iorb];
    } 
  }

  void TricubicBsplineSetBuilder::readData(const char* hroot, const vector<int>& occSet,
      int spinIndex, int degeneracy)
  {
    bool truncate = (myParam->Rcut>0.0);
    int norb=occSet.size()/degeneracy;
    StorageType inData(BoxGrid[0],BoxGrid[1],BoxGrid[2]);

    PosType center(0.0);
    for(int iorb=0; iorb<norb; iorb++) 
    {
      if(truncate)
      {
        string centerName=myParam->getCenterName(hroot,occSet[iorb]);
        HDFAttribIO<PosType > cdummy(activeBasis->Centers[iorb]);
        cdummy.read(hfileID,centerName.c_str());
      }
      ostringstream wnshort;
      string eigvName=myParam->getEigVectorName(hroot,occSet[iorb]/degeneracy,spinIndex);
      wnshort<<curH5Fname << "#"<<occSet[iorb]/degeneracy << "#" << spinIndex;
      map<string,StorageType*>::iterator it(BigDataSet.find(wnshort.str()));
      if(it == BigDataSet.end()) 
      {
        if(print_log)
        {
          app_log() << "   Reading spline function " << eigvName << " (" << wnshort.str()  << ")"  << endl;
          if(truncate) app_log() << "     center=" << center << endl;
        }
        StorageType* newP=new StorageType;
        HDFAttribIO<StorageType> dummy(inData);
        dummy.read(hfileID,eigvName.c_str());
        BigDataSet[wnshort.str()]=newP;
        activeBasis->add(iorb,center,inData,newP);
      } 
      else 
      {
        if(print_log)
        {
          app_log() << "   Reusing spline function " << eigvName << " (" << wnshort.str()  << ")" << endl;
          if(truncate) app_log() << "     center=" << center << endl;
        }
        activeBasis->add(iorb,center,(*it).second);
      }
    } 
  }

  void TricubicBsplineSetBuilder::readDataOMP(const char* hroot, const vector<int>& occSet,
      int spinIndex, int degeneracy)
  {

    if(BigDataSet.size()) //data exist, use the standard one
    {
      readData(hroot,occSet,spinIndex,degeneracy);
      return;
    }

    int norb=occSet.size()/degeneracy;
    vector<StorageType*> bsset(norb);
    vector<int> odist(omp_get_max_threads()+1);
    FairDivideLow(norb,omp_get_max_threads(),odist);

#pragma omp parallel
    {
      bool truncate = (myParam->Rcut>0.0);
      PosType center;
      StorageType inData(BoxGrid[0],BoxGrid[1],BoxGrid[2]);
      int ip=omp_get_thread_num();
      for(int iorb=odist[ip]; iorb<odist[ip+1]; iorb++)
      {
#pragma omp critical
        {
          string eigvName=myParam->getEigVectorName(hroot,occSet[iorb]/degeneracy,spinIndex);
          app_log() << "   Reading spline function " << ip << " "  << eigvName << endl;
          HDFAttribIO<StorageType> dummy(inData);
          dummy.read(hfileID,eigvName.c_str());
        }
        StorageType* newP =new StorageType;
        bsset[iorb]=newP;
        activeBasis->add(iorb,center,inData,newP);
      }
    }

    bool localize=myParam->Rcut>0.0;
    for(int iorb=0; iorb<norb; iorb++) 
    {
      if(localize)
      {
        string centerName=myParam->getCenterName(hroot,occSet[iorb]);
        HDFAttribIO<PosType > cdummy(activeBasis->Centers[iorb]);
        cdummy.read(hfileID,centerName.c_str());
      }
      ostringstream wnshort;
      wnshort<<curH5Fname << "#"<<occSet[iorb]/degeneracy << "#" << spinIndex;
      BigDataSet[wnshort.str()]=bsset[iorb];
    } 
  }
}
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 2082 $   $Date: 2007-06-15 16:25:49 -0500 (Fri, 15 Jun 2007) $
 * $Id: TricubicBsplineSetBuilder.cpp 2082 2007-06-15 21:25:49Z jnkim $ 
 ***************************************************************************/

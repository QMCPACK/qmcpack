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
    
    
    
    
    
    
    
    



#include "Estimators/PairCorrEstimator.h"
#include "Particle/DistanceTable.h"
#include "Particle/DistanceTableData.h"

namespace qmcplusplus
{

/** HDF5 output engine  which appends the current data to a dataset
 *
 * Raise the dimension of the working data (ref) by one. The first dimension
 * is set to H5S_UNLIMITED. Each write call results in extending
 * the HDF5 dataset by the size of the ref.
 * @todo use template specialization and look for a better name of the class.
 */
struct AppendData
{
  typedef Matrix<double> ArrayType_t;
  hid_t datasetID;
  hid_t dataspaceID;
  hsize_t stride;
  hsize_t maxdims[3];
  hsize_t curdims[3];
  hsize_t dims[3];
  hsize_t offset[3];
  ArrayType_t&  ref;

  AppendData(ArrayType_t& a, int nthreads=1, int ip=0):
    datasetID(-1), stride(nthreads), ref(a)
  {
    maxdims[0] = H5S_UNLIMITED;
    maxdims[1] = a.rows();
    maxdims[2] = a.cols();
    curdims[0] = 1;
    curdims[1] = a.rows();
    curdims[2] = a.cols();
    dims[0] = 1;
    dims[1] = a.rows();
    dims[2] = a.cols();
    offset[0]=ip;
    offset[1]=0;
    offset[2]=0;
  }

  ~AppendData()
  {
    if(datasetID>-1)
    {
      H5Sclose(dataspaceID);
      H5Dclose(datasetID);
    }
  }

  inline void write(hid_t grp, const char* name)
  {
    const hsize_t RANK=3;
    if(datasetID<0)
    {
      dataspaceID=H5Screate_simple(RANK, dims, maxdims);
      hid_t p = H5Pcreate (H5P_DATASET_CREATE);
      H5Pset_chunk(p,RANK,dims);
      datasetID= H5Dcreate(grp,name,H5T_NATIVE_DOUBLE,dataspaceID,p);
      hid_t memspace = H5Screate_simple(RANK, dims, NULL);
      hid_t ret = H5Dwrite(datasetID, H5T_NATIVE_DOUBLE, memspace, dataspaceID, H5P_DEFAULT, ref.data());
      H5Sclose(memspace);
      H5Pclose(p);
    }
    else
    {
      H5Dextend(datasetID,curdims);
      H5Sset_extent_simple(dataspaceID,RANK,curdims,maxdims);
      H5Sselect_hyperslab(dataspaceID, H5S_SELECT_SET, offset, NULL, dims, NULL);
      hid_t memspace = H5Screate_simple(RANK, dims, NULL);
      hid_t ret = H5Dwrite(datasetID, H5T_NATIVE_DOUBLE, memspace, dataspaceID, H5P_DEFAULT, ref.data());
      H5Sclose(memspace);
    }
    curdims[0]++;
    offset[0]+=stride;
  }

};

PairCorrEstimator::PairCorrEstimator(ParticleSet& source):
  Symmetric(true),sourcePtcl(source)
{
  myTable = DistanceTable::add(source);
  int ns=sourcePtcl.groups();
  std::vector<int> mask(ns*ns,-1);
  int ij=0;
  for(int i=0; i<ns; i++)
    for(int j=i; j<ns; j++,ij++)
    {
      mask[j+i*ns]=ij;
#if PRINT_DEBUG
      char fname[32];
      sprintf(fname,"gofr.%s_%d_%d.dat",myTable->Name.c_str(),i,j);
      fout.push_back(new std::ofstream(fname));
      fout[ij]->setf(std::ios::scientific, std::ios::floatfield);
      fout[ij]->precision(5);
#endif
    }
  NumPairTypes=ij;
  Centers=sourcePtcl.getTotalNum();
  PairID.resize(myTable->getTotNadj());
  for(int iat=0; iat<Centers; iat++)
  {
    for(int nn=myTable->M[iat]; nn<myTable->M[iat+1]; nn++)
    {
      PairID[nn]=mask[myTable->PairID[nn]];
    }
  }
  setBound(10,0.1);
}

PairCorrEstimator::PairCorrEstimator(const ParticleSet& source, ParticleSet& target):
  Symmetric(false),sourcePtcl(source)
{
  myTable = DistanceTable::add(source,target);
  NumPairTypes=sourcePtcl.getSpeciesSet().getTotalNum();
#if PRINT_DEBUG
  for(int i=0; i<NumPairTypes; i++)
  {
    char fname[32];
    sprintf(fname,"gofr.%s_%s.dat",myTable->Name.c_str(),
            sourcePtcl.getSpeciesSet().speciesName[i].c_str());
    fout.push_back(new std::ofstream(fname));
    fout[i]->setf(std::ios::scientific, std::ios::floatfield);
    fout[i]->precision(5);
  }
#endif
  Centers=sourcePtcl.getTotalNum();
  PairID.resize(myTable->getTotNadj());
  for(int iat=0; iat<Centers; iat++)
  {
    for(int nn=myTable->M[iat]; nn<myTable->M[iat+1]; nn++)
    {
      PairID[nn]=sourcePtcl.GroupID[iat];
    }
  }
  setBound(10,0.1);
}

PairCorrEstimator::~PairCorrEstimator()
{
}


void PairCorrEstimator::resetTargetParticleSet(ParticleSet& p)
{
  if(Symmetric)
    myTable=DistanceTable::add(p);
  else
    myTable=DistanceTable::add(sourcePtcl,p);
}

void PairCorrEstimator::open(hid_t hroot)
{
  if(GroupID<0)
  {
    Title="pc_"+myTable->Name;
    GroupID = H5Gcreate(hroot,Title.c_str(),0);
    v_h = new AppendData(gofr);
    v2_h = new AppendData(gofr2);
  }
}

void PairCorrEstimator::close()
{
  if(GroupID>-1)
  {
    delete v_h;
    delete v2_h;
    H5Gclose(GroupID);
    GroupID=-1;
  }
}

void PairCorrEstimator::startAccumulate()
{
  gofrInst=0.0;
}

/** accumulate the observables */
void PairCorrEstimator::accumulate(ParticleSet& p, RealType wgt)
{
  for(int iat=0; iat<Centers; iat++)
  {
    for(int nn=myTable->M[iat]; nn<myTable->M[iat+1]; nn++)
    {
      if(myTable->r(nn)>=Dmax)
        continue;
      gofrInst(PairID[nn],DeltaInv*myTable->r(nn))+=wgt;
    }
  }
}

/** add gofrInst which contains sum over walkers */
void PairCorrEstimator::stopAccumulate(RealType wgtinv)
{
  //gofr += wgtinv*gofrInst;
  collect(gofrInst.begin(),gofrInst.end(),gofr.begin(),gofr2.begin(),wgtinv);
}

/** save the block average */
void PairCorrEstimator::stopBlock(RealType wgtnorm, RealType errnorm)
{
  v_h->write(GroupID,"v");
  v2_h->write(GroupID,"v2");
#if PRINT_DEBUG
  for(int i=0; i<gofr.size1(); ++i)
  {
    RealType r=0.0;
    for(int j=0; j<gofr.size2(); ++j, r+=Delta)
    {
      RealType avg=gofr(i,j)*wgtnorm;
      *fout[i] << std::setw(15) << r << std::setw(15) << avg
               << std::setw(15) << (gofr2(i,j)*wgtnorm-avg*avg)*errnorm
               << std::endl;
    }
    *fout[i] << std::endl;
  }
#endif
}

void PairCorrEstimator::startBlock(int steps)
{
  gofrInst=0.0;
  gofr=0.0;
  gofr2=0.0;
}

void PairCorrEstimator::setBound(RealType rmax, RealType dr)
{
  Dmax=rmax;
  Delta=dr;
  DeltaInv=1.0/dr;
  int n=(Dmax)/Delta+1;
  gofrInst.resize(NumPairTypes,n);
  gofr.resize(NumPairTypes,n);
  gofr2.resize(NumPairTypes,n);
  gofrerr.resize(NumPairTypes,n);
}

//  /////////////////////////////////////////////////////
//  PairCorrEstimator::PairCorrEstimator(ParticleSet& source):
//    Symmetric(true),sourcePtcl(source), fout(0)
//  {
//    myTable = DistanceTable::add(source);
//    setBound(10,0.1);
//  }
//
//  PairCorrEstimator::PairCorrEstimator(const ParticleSet& source, ParticleSet& target):
//    Symmetric(false),sourcePtcl(source), fout(0)
//  {
//    myTable = DistanceTable::add(source,target);
//    setBound(10,0.1);
//  }
//
//  PairCorrEstimator::~PairCorrEstimator()
//  {
//  }
//
//  void PairCorrEstimator::resetTargetParticleSet(ParticleSet& p)
//  {
//    if(Symmetric)
//      myTable=DistanceTable::add(p);
//    else
//      myTable=DistanceTable::add(sourcePtcl,p);
//  }
//
//  void PairCorrEstimator::startAccumulate()
//  {
//    dCInst=0.0;
//  }
//
//  /** accumulate the observables */
//  void PairCorrEstimator::accumulate(ParticleSet& p, RealType wgt)
//  {
//    for(int i=0; i<myTable->getTotNadj(); i++)
//    {
//      if(myTable->r(i)<Dmax) dCInst[DeltaInv*myTable->r(i)]+=wgt;
//    }
//  }
//
//  /** reweight of the current cummulative  values */
//  void PairCorrEstimator::stopAccumulate(RealType wgtinv)
//  {
//    //RealType norm=wgtinv/static_cast<RealType>(myTable->getTotNadj());
//    //dCBlock += dCInst*norm;
//    dCBlock += dCInst*wgtinv;
//  }
//
//  void PairCorrEstimator::stopBlock(RealType wgtinv)
//  {
//    if(!fout)
//    {
//      char fname[32];
//      sprintf(fname,"gofr.%s.dat",myTable->Name.c_str());
//      fout = new std::ofstream(fname);
//      fout->setf(std::ios::scientific, std::ios::floatfield);
//      fout->precision(5);
//    }
//    RealType r=0.0;
//    for(int i=0; i<dCBlock.size(); i++, r+=Delta)
//      *fout << std::setw(15) << r << std::setw(15) << wgtinv*dCBlock[i] << std::endl;
//    *fout << std::endl;
//  }
//
//  void PairCorrEstimator::startBlock(int steps)
//  {
//    dCInst=0.0;
//    dCBlock=0.0;
//  }
//
//  void PairCorrEstimator::setBound(RealType rmax, RealType dr)
//  {
//    Dmax=rmax;
//    Delta=dr;
//    DeltaInv=1.0/dr;
//    int n=(Dmax)/Delta+1;
//    dCInst.resize(n);
//    dCBlock.resize(n);
//  }
}


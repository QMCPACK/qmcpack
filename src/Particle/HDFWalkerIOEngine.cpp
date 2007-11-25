//////////////////////////////////////////////////////////////////
// (c) Copyright 2004- by Jeongnim Kim and Jordan Vincent
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
/** @file HDFWalkerIOEngine.cpp
 * @breif definition  of HDFWalkerIOEngine
 */
#include "Particle/HDFWalkerIOEngine.h"
#include "Particle/MCWalkerConfiguration.h"
#include "Utilities/OhmmsInfo.h"
#include "Message/Communicate.h"

namespace qmcplusplus 
{

  HDFWalkerIOEngine::HDFWalkerIOEngine(MCWalkerConfiguration& a, bool reuse):
    W(a), replace(reuse) { }

  void HDFWalkerIOEngine::write(hid_t grp, const char* name) {
    typedef MCWalkerConfiguration::PosType PosType;
    Matrix<PosType> tp(W.getActiveWalkers(),W.R.size());
    //store walkers in a temporary array and pack them
    int item(0);
    MCWalkerConfiguration::const_iterator it(W.begin());
    MCWalkerConfiguration::const_iterator it_end(W.end());
    while(it != it_end) 
    {
      std::copy((*it)->R.begin(),(*it)->R.end(),tp[item++]);
      ++it; 
    }

    //now save to grp as a named object
    const int rank = 3;
    hsize_t dims[rank], maxdims[rank];
    dims[0] = tp.rows(); dims[1] = tp.cols(); dims[2] = OHMMS_DIM;
    maxdims[0] = H5S_UNLIMITED; maxdims[1] = tp.cols(); maxdims[2] = OHMMS_DIM;
    if(replace)
    {
      hid_t dataset =  H5Dopen(grp, name);
      //herr_t ret=H5Dextend(dataset,dims); //use dynamic resizing
      hid_t memspace = H5Screate_simple(rank, dims, NULL);
      herr_t ret = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, memspace, H5S_ALL, H5P_DEFAULT,tp.data());
      H5Sclose(memspace);
      H5Dclose(dataset);
    }
    else
    {//
      hid_t dataspace  = H5Screate_simple(rank, dims, maxdims);
      hid_t p = H5Pcreate (H5P_DATASET_CREATE);
      H5Pset_chunk(p,rank,dims);
      hid_t dataset =  H5Dcreate(grp, name, H5T_NATIVE_DOUBLE, dataspace, p);
      hid_t memspace = H5Screate_simple(rank, dims, NULL);
      hid_t ret = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace, H5P_DEFAULT,tp.data());
      H5Sclose(memspace);
      H5Sclose(dataspace);
      H5Dclose(dataset);
      //H5Pclose(p);
    }
  }

  void HDFWalkerIOEngine::writeAll(hid_t grp, const char* name, Communicate* myComm) 
  {
    int mynode=myComm->rank();
    int nwloc=W.getActiveWalkers();
    if(nwloc != W.WalkerOffsets[mynode+1]-W.WalkerOffsets[mynode])
    {
      app_error() << " Fatal Error: inconsident number of walkers per node " << endl;
      abort();//FIXABORT
    }

    typedef MCWalkerConfiguration::PosType PosType;
    Matrix<PosType> tp(nwloc,W.R.size());
    //store walkers in a temporary array and pack them
    int item(0);
    MCWalkerConfiguration::const_iterator it(W.begin());
    MCWalkerConfiguration::const_iterator it_end(W.end());
    while(it != it_end) 
    {
      std::copy((*it)->R.begin(),(*it)->R.end(),tp[item++]);
      ++it; 
    }

    const int RANK = 3;
    hsize_t offset[RANK],gcount[RANK],count[RANK],stride[]={1,1,1};
    count[0] = nwloc; count[1] = W.getTotalNum(); count[2] = OHMMS_DIM;
    gcount[0] = W.getGlobalNumWalkers(); gcount[1] = W.getTotalNum(); gcount[2] = OHMMS_DIM;
    offset[0] = W.WalkerOffsets[mynode]; offset[1] = 0; offset[2] = 0;

    hid_t dset_id;
    if(replace)
    {
      dset_id=H5Dopen(grp,name);
    }
    else
    { //create global dataset
      hid_t sid1  = H5Screate_simple(RANK,gcount,NULL);
      dset_id=H5Dcreate(grp,name,H5T_NATIVE_DOUBLE,sid1,H5P_DEFAULT);
      H5Sclose(sid1);
    }

    hid_t memspace=H5Screate_simple(RANK,count,NULL);
    hid_t filespace=H5Dget_space(dset_id);
    herr_t ret=H5Sselect_hyperslab(filespace,H5S_SELECT_SET,offset,stride,count,NULL);

    //transfer method is provided by the users 
    ret = H5Dwrite(dset_id,H5T_NATIVE_DOUBLE,memspace,filespace,xfer_plist,tp.data());

    H5Sclose(filespace);
    H5Sclose(memspace);
    H5Dclose(dset_id);
  }

  void HDFWalkerIOEngine::read(hid_t grp, const char* name) 
  {
    const int RANK = 3;
    hsize_t count[RANK], offset[]={0,0,0};

    hid_t dataset = H5Dopen(grp,name);
    hid_t dataspace = H5Dget_space(dataset);
    int rank = H5Sget_simple_extent_ndims(dataspace);
    int status_n = H5Sget_simple_extent_dims(dataspace, count, NULL);

    vector<MCWalkerConfiguration::PosType> posIn(count[0]*count[1]);

    hid_t memspace = H5Screate_simple(RANK, count, NULL);
    herr_t status = H5Sselect_hyperslab(dataspace,H5S_SELECT_SET, offset,NULL,count,NULL);
    status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace, H5P_DEFAULT, &(posIn[0][0]));

    H5Sclose(dataspace);
    H5Sclose(memspace);
    H5Dclose(dataset);

    int curWalker = W.getActiveWalkers();
    if(curWalker) {
      W.createWalkers(count[0]);
    } else {
      W.resize(count[0],count[1]);
    }

    MCWalkerConfiguration::iterator it = W.begin()+curWalker; 
    int ii=0;
    for(int iw=0; iw<count[0]; iw++) {
      //std::copy(Post_temp[iw],Post_temp[iw+1], (*it)->R.begin());
      for(int iat=0; iat < count[1]; iat++,ii++){
        (*it)->R(iat) = posIn[ii];
      }
      ++it;
    }
  }

  void HDFWalkerIOEngine::readAll(hid_t grp, const char* name, Communicate* comm) 
  {
    int mynode=comm->rank();
    int nprocs=comm->size();
    vector<int> woffset(nprocs+1,0);

    const int RANK = 3;
    hsize_t offset[]={1,1,1};
    hsize_t gcount[RANK],count[RANK];
    hsize_t stride[]={1,1,1};

    hid_t dataset = H5Dopen(grp,name);
    hid_t dataspace = H5Dget_space(dataset);
    int rank_n = H5Sget_simple_extent_ndims(dataspace);
    int status_n = H5Sget_simple_extent_dims(dataspace, gcount, NULL);

    //assign offsets and size
    FairDivideLow(gcount[0],nprocs,woffset);
    offset[0]=woffset[mynode]; offset[1] = 0; offset[2] = 0;
    count[0]=woffset[mynode+1]-woffset[mynode];
    count[1]=gcount[1];
    count[2]=gcount[2];

    app_log() << "   Initial walker distribution: ";
    std::copy(woffset.begin(),woffset.end(),ostream_iterator<int>(app_log()," "));
    app_log() << endl;

    vector<MCWalkerConfiguration::PosType> posIn(count[0]*count[1]);

    hid_t memspace = H5Screate_simple(RANK, count, NULL);
    herr_t status = H5Sselect_hyperslab(dataspace,H5S_SELECT_SET, offset,NULL,count,NULL);

#if defined(H5_HAVE_PARALLEL)
    hid_t  xfer_plist = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(xfer_plist,H5FD_MPIO_COLLECTIVE);
#else
    hid_t  xfer_plist =  H5P_DEFAULT;
#endif
    //status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace, H5P_DEFAULT, &(posIn[0][0]));
    status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace, xfer_plist, &(posIn[0][0]));

    H5Sclose(dataspace);
    H5Sclose(memspace);
    H5Dclose(dataset);

#if defined(H5_HAVE_PARALLEL)
    H5Pclose(xfer_plist);
#endif

    int curWalker = W.getActiveWalkers();
    if(curWalker) {
      W.createWalkers(count[0]);
    } else {
      W.resize(count[0],count[1]);
    }

    MCWalkerConfiguration::iterator it = W.begin()+curWalker; 
    int ii=0;
    for(int iw=0; iw<count[0]; iw++) {
      //std::copy(Post_temp[iw],Post_temp[iw+1], (*it)->R.begin());
      for(int iat=0; iat < count[1]; iat++,ii++){
        (*it)->R(iat) = posIn[ii];
      }
      ++it;
    }
  }
}
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 2303 $   $Date: 2007-11-19 13:25:50 -0600 (Mon, 19 Nov 2007) $
 * $Id: HDFWalkerOutput.cpp 2303 2007-11-19 19:25:50Z jnkim $ 
 ***************************************************************************/

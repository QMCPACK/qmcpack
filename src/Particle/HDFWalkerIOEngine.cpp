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
    
    



/** @file HDFWalkerIOEngine.cpp
 * @brief definition  of HDFWalkerIOEngine
 */
#include "Particle/HDFWalkerIOEngine.h"
#include "Particle/MCWalkerConfiguration.h"
#include "Message/Communicate.h"

namespace qmcplusplus
{

HDFWalkerIOEngine::HDFWalkerIOEngine(MCWalkerConfiguration& a, bool reuse):
  W(a), replace(reuse) { }

void HDFWalkerIOEngine::read(hid_t grp, const char* name)
{
  const int RANK = 3;
  hsize_t count[RANK], offset[]= {0,0,0};
  hid_t dataset = H5Dopen(grp,name);
  hid_t dataspace = H5Dget_space(dataset);
  int rank = H5Sget_simple_extent_ndims(dataspace);
  int status_n = H5Sget_simple_extent_dims(dataspace, count, NULL);
  std::vector<MCWalkerConfiguration::PosType> posIn(count[0]*count[1]);
  hid_t memspace = H5Screate_simple(RANK, count, NULL);
  herr_t status = H5Sselect_hyperslab(dataspace,H5S_SELECT_SET, offset,NULL,count,NULL);
  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace, H5P_DEFAULT, &(posIn[0][0]));
  H5Sclose(dataspace);
  H5Sclose(memspace);
  H5Dclose(dataset);
  int curWalker = W.getActiveWalkers();
  if(curWalker)
  {
    W.createWalkers(count[0]);
  }
  else
  {
    W.resize(count[0],count[1]);
  }
  MCWalkerConfiguration::iterator it = W.begin()+curWalker;
  int ii=0;
  for(int iw=0; iw<count[0]; iw++)
  {
    //std::copy(Post_temp[iw],Post_temp[iw+1], (*it)->R.begin());
    for(int iat=0; iat < count[1]; iat++,ii++)
    {
      (*it)->R[iat] = posIn[ii];
    }
    ++it;
  }
}

void HDFWalkerIOEngine::readAll(hid_t grp, const char* name, Communicate* comm)
{
  int mynode=comm->rank();
  int nprocs=comm->size();
  std::vector<int> woffset(nprocs+1,0);
  const int RANK = 3;
  hsize_t offset[]= {1,1,1};
  hsize_t gcount[RANK],count[RANK];
  hsize_t stride[]= {1,1,1};
  hid_t dataset = H5Dopen(grp,name);
  hid_t dataspace = H5Dget_space(dataset);
  int rank_n = H5Sget_simple_extent_ndims(dataspace);
  int status_n = H5Sget_simple_extent_dims(dataspace, gcount, NULL);
  //assign offsets and size
  FairDivideLow(gcount[0],nprocs,woffset);
  offset[0]=woffset[mynode];
  offset[1] = 0;
  offset[2] = 0;
  count[0]=woffset[mynode+1]-woffset[mynode];
  count[1]=gcount[1];
  count[2]=gcount[2];
  app_log() << "   Initial walker distribution: ";
  copy(woffset.begin(),woffset.end(),std::ostream_iterator<int>(app_log()," "));
  app_log() << std::endl;
  std::vector<MCWalkerConfiguration::PosType> posIn(count[0]*count[1]);
  hid_t memspace = H5Screate_simple(RANK, count, NULL);
  herr_t status = H5Sselect_hyperslab(dataspace,H5S_SELECT_SET, offset,NULL,count,NULL);
#if defined(H5_HAVE_PARALLEL)
  xfer_plist = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(xfer_plist,H5FD_MPIO_COLLECTIVE);
#else
  xfer_plist =  H5P_DEFAULT;
#endif
  hid_t type_id=get_h5_datatype(posIn[0][0]);
  status = H5Dread(dataset, type_id, memspace, dataspace, xfer_plist, &(posIn[0][0]));
  H5Sclose(dataspace);
  H5Sclose(memspace);
  H5Dclose(dataset);
#if defined(H5_HAVE_PARALLEL)
  H5Pclose(xfer_plist);
#endif
  int curWalker = W.getActiveWalkers();
  if(curWalker)
  {
    W.createWalkers(count[0]);
  }
  else
  {
    W.resize(count[0],count[1]);
  }
  MCWalkerConfiguration::iterator it = W.begin()+curWalker;
  int ii=0;
  for(int iw=0; iw<count[0]; iw++)
  {
    //std::copy(Post_temp[iw],Post_temp[iw+1], (*it)->R.begin());
    for(int iat=0; iat < count[1]; iat++,ii++)
    {
      (*it)->R[iat] = posIn[ii];
    }
    ++it;
  }
}
}

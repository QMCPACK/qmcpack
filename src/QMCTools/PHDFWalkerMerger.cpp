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
    
    



#include <vector>
#include <iostream>
#include "QMCTools/HDFWalkerMerger.h"
#include "Utilities/RandomGeneratorIO.h"
#include "Utilities/UtilityFunctions.h"
#include "Message/Communicate.h"
#include "Message/CommOperators.h"

using namespace qmcplusplus;


HDFWalkerMerger::HDFWalkerMerger(const std::string& aroot, int ncopy)
  : NumCopy(ncopy), FileRoot(aroot)
{
}

HDFWalkerMerger::~HDFWalkerMerger() { }

void HDFWalkerMerger::setCommunicator(Communicate* c)
{
  myComm=c;
}

void HDFWalkerMerger::init()
{
  char fname[128];
  char GrpName[128];
  int nprocs=myComm->ncontexts();
  int mynode=myComm->mycontext();
  hsize_t dimin[3], summary_size=3;
  std::vector<int> nw_pernode(nprocs,0);
  std::vector<int> nc_offset(nprocs+1,0);
  FairDivideLow(NumCopy,nprocs,nc_offset);
  std::cout << "Number of copies " << NumCopy << std::endl;
  std::cout << "Number of proces " << nprocs << std::endl;
  for(int i=0; i<=nprocs; i++)
    std::cout << nc_offset[i] << std::endl;
  numWalkersIn.push_back(new std::vector<hsize_t>);
  //determine the minimum configuration and walkers for each config
  int min_config=1;
  for(int ip=nc_offset[mynode]; ip<nc_offset[mynode+1]; ip++)
  {
    int nconf=1;
    sprintf(fname,"%s.p%03d.config.h5",FileRoot.c_str(),ip);
    hid_t hfile = H5Fopen(fname,H5F_ACC_RDWR,H5P_DEFAULT);
    hid_t h_config = H5Gopen(hfile,"config_collection");
    hid_t h1=H5Dopen(h_config,"NumOfConfigurations");
    if(h1>-1)
    {
      H5Dread(h1, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&(nconf));
      H5Dclose(h1);
    }
    //iconf is the last one
    int iconf=nconf-1;
    sprintf(GrpName,"config%04d/coord",iconf);
    h1 = H5Dopen(h_config,GrpName);
    hid_t dataspace = H5Dget_space(h1);
    int rank = H5Sget_simple_extent_ndims(dataspace);
    int status_n = H5Sget_simple_extent_dims(dataspace, dimin, NULL);
    H5Sclose(dataspace);
    H5Dclose(h1);
    nw_pernode[mynode]+=dimin[0];
    numWalkersIn[0]->push_back(dimin[0]);
    H5Fclose(hfile);
  }
  myComm->allreduce(nw_pernode);
  NumPtcl = dimin[1];
  Dimension = dimin[2];
  OffSet.resize(2,nprocs+1);
  OffSet=0;
  int ip_in=0;
  for(int ip=0; ip<nprocs; ip++)
  {
    OffSet(0,ip+1)=OffSet(0,ip)+nw_pernode[ip];
    if(mynode ==0)
    {
      std::cout << ip << " has " << nw_pernode[ip] << std::endl;
    }
  }
  for(int ip=0; ip<=nprocs; ip++)
    OffSet(1,ip)=nc_offset[ip];
  OffSet(nprocs,1)=nc_offset[nprocs];
  MaxNumWalkers=nw_pernode[mynode];
  NumConfig=1;
}

void HDFWalkerMerger::writeHeader(hid_t gid)
{
  int cur_version=1;
  hsize_t dim=1;
  hid_t header_space= H5Screate_simple(1, &dim, NULL);
  hid_t header_set= H5Dcreate(gid, "version", H5T_NATIVE_INT, header_space, H5P_DEFAULT);
  H5Dwrite(header_set, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&cur_version);
  H5Dclose(header_set);
  header_set= H5Dcreate(gid, "ncontexts", H5T_NATIVE_INT, header_space, H5P_DEFAULT);
  H5Dwrite(header_set, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&NumCopy);
  H5Dclose(header_set);
  H5Sclose(header_space);
}

void HDFWalkerMerger::merge()
{
  init();
  char fname[128];
  char GrpName[128];
  char ConfName[128];
#if H5_VERS_RELEASE < 4
  hssize_t offset[]= {0,0,0};
#else
  hsize_t offset[]= {0,0,0};
#endif
  typedef std::vector<OHMMS_PRECISION> Container_t;
  Container_t vecLoc(MaxNumWalkers*NumPtcl*Dimension);
  hsize_t dimout[3], dimin[3];
  dimin[1]=NumPtcl;
  dimin[2]=Dimension;
  dimout[1]=NumPtcl;
  dimout[2]=Dimension;
  ////create new dataspace for each config
  //for(int iconf=0; iconf<NumConfig; iconf++) {
  //  dimout[0]=OffSet(NumCopy,iconf);
  //  sprintf(ConfName,"config%04d",iconf);
  //  hid_t cf_id = H5Gcreate(mastercf,ConfName,0);
  //  dataspace = H5Screate_simple(3,dimout,NULL);
  //  dataset = H5Dcreate(cf_id,"coord",H5T_NATIVE_DOUBLE, dataspace,H5P_DEFAULT);
  //  H5Dclose(dataset);
  //  H5Sclose(dataspace);
  //  H5Gclose(cf_id);
  //}
  int nprocs=myComm->ncontexts();
  int mynode=myComm->mycontext();
  int ndacc=0; //local data count
  for(int ip=OffSet(1,mynode),ipL=0; ip<OffSet(1,mynode+1); ip++,ipL++)
  {
    sprintf(fname,"%s.p%03d.config.h5",FileRoot.c_str(),ip);
    hid_t h0 = H5Fopen(fname,H5F_ACC_RDWR,H5P_DEFAULT);
    hid_t cfcollect = H5Gopen(h0,"config_collection");
    //get the number of walkers for (ip,iconf)
    dimin[0] = (*numWalkersIn[0])[ipL];
    sprintf(ConfName,"config%04d",0);
    std::cout << "Node " << mynode << " reads " << dimin[0] << " walkers from " << fname <<  std::endl;
    int nd=dimin[0]*NumPtcl*Dimension;
    hid_t cf_id = H5Gopen(cfcollect,ConfName);
    Container_t vecIn(nd);
    HDFAttribIO<Container_t> in(vecIn, 3, dimin);
    in.read(cf_id,"coord");
    H5Gclose(cf_id);
    copy(vecIn.begin(),vecIn.begin()+nd,vecLoc.begin()+ndacc);
    ndacc+=nd;
    H5Fclose(h0);
  }
  MPI_Info info=MPI_INFO_NULL;
  hid_t acc_tpl1=H5Pcreate(H5P_FILE_ACCESS);
  herr_t ret=H5Pset_fapl_mpio(acc_tpl1,myComm->getMPI(),info);
  sprintf(fname,"%s.config.h5",FileRoot.c_str());
  hid_t masterfile= H5Fcreate(fname,H5F_ACC_TRUNC, H5P_DEFAULT,acc_tpl1);
  ret=H5Pclose(acc_tpl1);
  dimout[0]=OffSet(0,nprocs);
  hid_t mastercf = H5Gcreate(masterfile,"config_collection",0);
  hid_t c1 = H5Gcreate(mastercf,"config0000",0);
  //create global dimensionality
  hid_t sid1=H5Screate_simple(3,dimout,NULL);
  //create dataset collectively
  hid_t dataset1=H5Dcreate(c1,"coord",H5T_NATIVE_DOUBLE,sid1,H5P_DEFAULT);
  hsize_t start[3],count[3],stride[3]= {1,1,1};
  start[0]=OffSet(0,mynode);
  start[1]=0;
  start[2]=0;
  count[0]=OffSet(0,mynode+1)-OffSet(0,mynode);
  count[1]=NumPtcl;
  count[2]=Dimension;
  //get the database
  hid_t file_dataspace=H5Dget_space(dataset1);
  ret=H5Sselect_hyperslab(file_dataspace,H5S_SELECT_SET,start,stride,count,NULL);
  //create data independently
  hid_t mem_dataspace=H5Screate_simple(3,count,NULL);
  ret=H5Dwrite(dataset1,H5T_NATIVE_DOUBLE,mem_dataspace,file_dataspace,H5P_DEFAULT, &vecLoc[0]);
  H5Sclose(file_dataspace);
  ret = H5Dclose(dataset1);
  H5Sclose(sid1);
  H5Gclose(c1);
  H5Gclose(mastercf);
  H5Fclose(masterfile);
  ////write the summaryset
  //double d=1.0/static_cast<double>(NumCopy);
  //for(int k=0; k<Summary.size(); k++) Summary[k]*=d;
  //hsize_t dim=Summary.size();
  //hid_t summary_space  = H5Screate_simple(1, &dim, NULL);
  //hid_t summary_set =  H5Dcreate(mastercf, "Summary", H5T_NATIVE_DOUBLE, summary_space, H5P_DEFAULT);
  //ret = H5Dwrite(summary_set, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,&Summary[0]);
  //H5Dclose(summary_set);
  //H5Sclose(summary_space);
  //// H5Gclose(mastercf);
  //H5Fclose(masterfile);
}

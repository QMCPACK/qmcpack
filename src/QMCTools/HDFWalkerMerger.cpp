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
using namespace qmcplusplus;


HDFWalkerMerger::HDFWalkerMerger(const std::string& aroot, int ncopy)
  : NumCopy(ncopy), FileRoot(aroot)
{
}

HDFWalkerMerger::~HDFWalkerMerger() { }

void HDFWalkerMerger::init()
{
  char fname[128];
  char GrpName[128];
  hsize_t dimin[3], summary_size=3;
  std::vector<double> summaryIn;
  //determine the minimum configuration and walkers for each config
  int min_config=10000;
  MaxNumWalkers=0;
  for(int ip=0; ip<NumCopy; ip++)
  {
    int nconf=1;
    numWalkersIn.push_back(new std::vector<hsize_t>);
    sprintf(fname,"%s.p%03d.config.h5",FileRoot.c_str(),ip);
    hid_t hfile = H5Fopen(fname,H5F_ACC_RDWR,H5P_DEFAULT);
    hid_t h_config = H5Gopen(hfile,"config_collection");
    hid_t h1=H5Dopen(h_config,"NumOfConfigurations");
    if(h1>-1)
    {
      H5Dread(h1, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&(nconf));
      H5Dclose(h1);
    }
    for(int iconf=0; iconf<nconf; iconf++)
    {
      sprintf(GrpName,"config%04d/coord",iconf);
      h1 = H5Dopen(h_config,GrpName);
      hid_t dataspace = H5Dget_space(h1);
      int rank = H5Sget_simple_extent_ndims(dataspace);
      int status_n = H5Sget_simple_extent_dims(dataspace, dimin, NULL);
      H5Sclose(dataspace);
      H5Dclose(h1);
      numWalkersIn[ip]->push_back(dimin[0]);
      MaxNumWalkers = std::max(MaxNumWalkers,dimin[0]);
    }
    h1 = H5Dopen(h_config,"Summary");
    hid_t summary_space = H5Dget_space(h1);
    H5Sget_simple_extent_ndims(summary_space);
    H5Sget_simple_extent_dims(summary_space,&summary_size, NULL);
    if(Summary.size() != summary_size)
    {
      Summary.resize(summary_size,0.0);
      summaryIn.resize(summary_size,0.0);
    }
    H5Dread(h1, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &summaryIn[0]);
    for(hsize_t k=0; k<summary_size; k++)
    {
      Summary[k] += summaryIn[k];
    }
    H5Sclose(summary_space);
    H5Dclose(h1);
    H5Fclose(hfile);
    min_config = std::min(min_config,nconf);
  }
  NumPtcl = dimin[1];
  Dimension = dimin[2];
  NumConfig = min_config;
  OffSet.resize(NumCopy+1,NumConfig);
  OffSet=0;
  for(int ip=0; ip<NumCopy; ip++)
  {
    for(int iconf=0; iconf<NumConfig; iconf++)
    {
      OffSet(ip+1,iconf) = OffSet(ip,iconf)+(*numWalkersIn[ip])[iconf];
    }
  }
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
  Container_t vecIn(MaxNumWalkers*NumPtcl*Dimension);
  sprintf(fname,"%s.config.h5",FileRoot.c_str());
  hid_t masterfile= H5Fcreate(fname,H5F_ACC_TRUNC, H5P_DEFAULT,H5P_DEFAULT);
  writeHeader(masterfile);
  hid_t mastercf = H5Gcreate(masterfile,"config_collection",0);
  hsize_t onedim=1;
  hid_t dataspace= H5Screate_simple(1, &onedim, NULL);
  hid_t dataset= H5Dcreate(mastercf, "NumOfConfigurations", H5T_NATIVE_INT, dataspace, H5P_DEFAULT);
  hid_t ret = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&NumConfig);
  H5Dclose(dataset);
  H5Sclose(dataspace);
  hsize_t dimout[3], dimin[3];
  dimin[1]=NumPtcl;
  dimin[2]=Dimension;
  dimout[1]=NumPtcl;
  dimout[2]=Dimension;
  //create new dataspace for each config
  for(int iconf=0; iconf<NumConfig; iconf++)
  {
    dimout[0]=OffSet(NumCopy,iconf);
    sprintf(ConfName,"config%04d",iconf);
    hid_t cf_id = H5Gcreate(mastercf,ConfName,0);
    dataspace = H5Screate_simple(3,dimout,NULL);
    dataset = H5Dcreate(cf_id,"coord",H5T_NATIVE_DOUBLE, dataspace,H5P_DEFAULT);
    H5Dclose(dataset);
    H5Sclose(dataspace);
    H5Gclose(cf_id);
  }
  //use random_states opposed to random_states
  for(int ip=0; ip<NumCopy; ip++)
  {
    sprintf(fname,"%s.p%03d.config.h5",FileRoot.c_str(),ip);
    hid_t h0 = H5Fopen(fname,H5F_ACC_RDWR,H5P_DEFAULT);
    hid_t cfcollect = H5Gopen(h0,"config_collection");
    const std::vector<hsize_t>& nw(*(numWalkersIn[ip]));
    for(int iconf=0; iconf<NumConfig; iconf++)
    {
      //get the number of walkers for (ip,iconf)
      dimin[0] = nw[iconf];
      sprintf(ConfName,"config%04d",iconf);
      hid_t cf_id = H5Gopen(cfcollect,ConfName);
      HDFAttribIO<Container_t> in(vecIn, 3, dimin);
      in.read(cf_id,"coord");
      H5Gclose(cf_id);
      //change the offset
      offset[0]=OffSet(ip,iconf);
      //create simple data to write
      hid_t memspace = H5Screate_simple(3, dimin, NULL);
      cf_id = H5Gopen(mastercf,ConfName);
      hid_t dataset = H5Dopen(cf_id,"coord");
      hid_t filespace = H5Dget_space(dataset);
      herr_t status = H5Sselect_hyperslab(filespace,H5S_SELECT_SET, offset,NULL,dimin,NULL);
      status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, memspace, filespace, H5P_DEFAULT, &vecIn[0]);
      H5Dclose(dataset);
      H5Sclose(filespace);
      H5Sclose(memspace);
      H5Gclose(cf_id);
    }
    H5Gclose(cfcollect);
    sprintf(GrpName,"context%04d",ip);
    hid_t mycontext = H5Gcreate(masterfile,GrpName,0);
    hid_t ranIn = H5Gopen(h0,"random_state");
    HDFAttribIO<RandomGenerator_t> r(Random);
    r.read(ranIn,"dummy");
    H5Gclose(ranIn);
    hid_t ranOut=H5Gcreate(mycontext, "random_state",0);
    r.write(ranOut,"dummy");
    H5Gclose(ranOut);
    H5Gclose(mycontext);
    H5Fclose(h0);
  }
  //write the summaryset
  double d=1.0/static_cast<double>(NumCopy);
  for(int k=0; k<Summary.size(); k++)
    Summary[k]*=d;
  hsize_t dim=Summary.size();
  hid_t summary_space  = H5Screate_simple(1, &dim, NULL);
  hid_t summary_set =  H5Dcreate(mastercf, "Summary", H5T_NATIVE_DOUBLE, summary_space, H5P_DEFAULT);
  ret = H5Dwrite(summary_set, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,&Summary[0]);
  H5Dclose(summary_set);
  H5Sclose(summary_space);
  H5Gclose(mastercf);
  H5Fclose(masterfile);
}

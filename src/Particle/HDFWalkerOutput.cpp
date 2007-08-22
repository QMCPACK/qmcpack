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
#include "Particle/MCWalkerConfiguration.h"
#include "Particle/HDFWalkerOutput.h"
#include "OhmmsPETE/OhmmsVector.h"
#include "Particle/HDFParticleAttrib.h"
#include "Numerics/HDFNumericAttrib.h"
#include "Utilities/OhmmsInfo.h"
#include "Utilities/RandomGenerator.h"
#include "Utilities/RandomGeneratorIO.h"
#include "OhmmsData/FileUtility.h"
#include "Message/Communicate.h"
//#define ENABLE_PHDFLIB
using namespace qmcplusplus;


/** specialized HDFAttribIOBase to write the positions of all the walkers.
 *
 * This class is to optimize writing walker configurations for restart and
 * for variational optimizations.
 */
struct WalkerPackedWriter: public HDFAttribIOBase {

  ///reference to the walkers
  MCWalkerConfiguration& W;
  ///if true, the content is replaced
  bool replace;

  WalkerPackedWriter(MCWalkerConfiguration& a, bool reuse=false):W(a), replace(reuse) { }

  inline void write(hid_t grp, const char* name) {
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
      //herr_t status=H5Dextend(dataset,dims);
      hid_t memspace = H5Screate_simple(rank, dims, NULL);
      hid_t ret = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, memspace, H5S_ALL, H5P_DEFAULT,tp.data());
      H5Sclose(memspace);
      H5Dclose(dataset);
    }
    else
    {
      hid_t dataspace  = H5Screate_simple(rank, dims, maxdims);
      hid_t p = H5Pcreate (H5P_DATASET_CREATE);
      H5Pset_chunk(p,rank,dims);
      hid_t dataset =  H5Dcreate(grp, name, H5T_NATIVE_DOUBLE, dataspace, p);
      hid_t memspace = H5Screate_simple(rank, dims, NULL);
      hid_t ret = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace, H5P_DEFAULT,tp.data());
      H5Sclose(memspace);
      H5Sclose(dataspace);
      H5Dclose(dataset);
      H5Pclose(p);
    }
  }

  inline void read(hid_t grp, const char* name) 
  {
  }

#if defined(HAVE_MPI)
  //collective writing: not working yet
  inline void writeAll(hid_t grp, const char* name, Communicate* myComm) 
  {
    int mynode=myComm->mycontext();
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


    //now save to grp as a named object
    const int rank = 3;
    hsize_t start[rank],count[rank],stride[]={1,1,1};
    start[0] = W.WalkerOffsets[mynode]; start[1] = 0; start[2] = 0;
    count[0] = tp.rows(); count[1] = tp.cols(); count[2] = OHMMS_DIM;


    hid_t dataset1;
    hsize_t gcount[rank];
    gcount[0] = W.getGlobalNumWalkers(); gcount[1] = tp.cols(); gcount[2] = OHMMS_DIM;

    hid_t sid1  = H5Screate_simple(rank,gcount,NULL);
    if(replace)
    {
      dataset1 =  H5Dopen(grp, name);
    }
    else
    {
      dataset1=H5Dcreate(grp,name,H5T_NATIVE_DOUBLE,sid1,H5P_DEFAULT);
    }

    hid_t file_dataspace=H5Dget_space(dataset1);
    herr_t  ret=H5Sselect_hyperslab(file_dataspace,H5S_SELECT_SET,start,stride,count,NULL);

    hid_t mem_dataspace = H5Screate_simple(rank, count, NULL);
    //hid_t xfer_plist = H5Pcreate(H5P_DATASET_XFER);
    //ret = H5Pset_dxpl_mpio(xfer_plist,H5FD_MPIO_COLLECTIVE);
    //ret = H5Dwrite(dataset1, H5T_NATIVE_DOUBLE, memspace, file_dataspace, xfer_plist, tp.data());
    ret = H5Dwrite(dataset1, H5T_NATIVE_DOUBLE, mem_dataspace, file_dataspace, H5P_DEFAULT, &(tp(0,0)[0]));

    H5Sclose(mem_dataspace);
    H5Sclose(file_dataspace);
    H5Dclose(dataset1);
    H5Sclose(sid1);

    //H5Pclose(xfer_plist);
  }
#endif
};

#if defined(ENABLE_PHDFLIB)

HDFWalkerOutput::HDFWalkerOutput(MCWalkerConfiguration& W,
    const string& aroot,Communicate* c): 
  NumOfWalkers(0), h_file(-1), h_config(-1), myComm(c)
{

  myComm = OHMMS::Controller;

  //h5FileName=aroot+".config.h5";
  h5FileName="phdftest.h5";

  MPI_Info info=MPI_INFO_NULL;
  hid_t acc_tpl1=H5Pcreate(H5P_FILE_ACCESS);
  herr_t ret=H5Pset_fapl_mpio(acc_tpl1,myComm->getMPI(),info);
  h_file = H5Fcreate(h5FileName.c_str(), H5F_ACC_TRUNC,H5P_DEFAULT,acc_tpl1);
  h_config = H5Gcreate(h_file,"config_collection",0); 
  ret=H5Pclose(acc_tpl1);

  if(myComm->master())
  {
    int nc(1); 
    HDFAttribIO<int> i(nc);
    i.write(h_config,"NumOfConfigurations");
  }

  //////create config0000 
  //char GrpName[128];
  //sprintf(GrpName,"config%04d",0);
  //hid_t group_id = H5Gcreate(h_config,GrpName,0);
  //WalkerPackedWriter wo(W);
  //wo.writeAll(group_id,"coord",myComm);
  //H5Gclose(group_id);

  //hid_t h_random = H5Gcreate(h_file,"random_state",0);
  //HDFAttribIO<RandomGenerator_t> r(Random);
  //r.write(h_random,"dummy");
  //H5Gclose(h_random);
  //H5Gclose(h_config);
  //H5Fclose(h_file); 
}

void HDFWalkerOutput::open()
{
  //MPI_Info info=MPI_INFO_NULL;
  //hid_t acc_tpl1=H5Pcreate(H5P_FILE_ACCESS);
  //herr_t ret=H5Pset_fapl_mpio(acc_tpl1,myComm->getMPI(),info);
  //h_file =  H5Fopen(h5FileName.c_str(),H5F_ACC_RDWR,acc_tpl1);
  //h_config = H5Gopen(h_file,"config_collection");
  //ret=H5Pclose(acc_tpl1);
}

/** Destructor writes the state of random numbers and close the file */
HDFWalkerOutput::~HDFWalkerOutput() {
  //if(h_file<0)
  //  h_file =  H5Fopen(h5FileName.c_str(),H5F_ACC_RDWR,H5P_DEFAULT);
  //hid_t h_random = H5Gopen(h_file,"random_state");
  //HDFAttribIO<RandomGenerator_t> r(Random);
  //r.write(h_random,"dummy");
  //H5Fclose(h_file);
  //
  if(h_config>-1) H5Gclose(h_config); 
  if(h_file>-1) H5Fclose(h_file); 
}

/** Write the set of walker configurations to the HDF5 file.  
 * @param W set of walker configurations
 */
bool HDFWalkerOutput::dump(MCWalkerConfiguration& W) {

  return true;
  const bool replace=true;
  char GrpName[128];
  sprintf(GrpName,"config%04d",0);

  hid_t group_id = H5Gopen(h_config,GrpName);
  if(NumOfWalkers != W.getGlobalNumWalkers())//need resize of the data
  {
    herr_t status=H5Gunlink(group_id,"coord");
    WalkerPackedWriter wo(W);
    wo.writeAll(group_id,"coord",myComm);
    NumOfWalkers=W.getGlobalNumWalkers();
  } 
  else
  {//simply overwrite it
    WalkerPackedWriter wo(W,replace);
    wo.writeAll(group_id,"coord",myComm);
  }
  H5Gclose(group_id);

  int c(1); 
  HDFAttribIO<int> i(c,replace);
  i.write(h_config,"NumOfConfigurations");

  return true;
}

/** Write the set of walker configurations to the HDF5 file.  
 * @param W set of walker configurations
 * @param block index for the block
 */
bool  HDFWalkerOutput::append(MCWalkerConfiguration& W, int block) {

  if(block==0) {
    return dump(W);
  }

  char GrpName[128];
  sprintf(GrpName,"config%04d",block);
  hid_t group_id = H5Gcreate(h_config,GrpName,0);
  WalkerPackedWriter wo(W);
  wo.write(group_id,"coord");
  H5Gclose(group_id);

  block++;
  
  HDFAttribIO<int> i(block,true);
  i.write(h_config,"NumOfConfigurations");
  return true;
}
#else
/** Create the HDF5 file "aroot.config.h5" for output. 
 * @param W walkers to operate on
 * @param aroot the root file name
 *
 * Constructor is responsible to create a hdf5 file and create basic groups
 * so that subsequent write can utilize existing dataspaces.
 * HDF5 contains
 * - config_collection
 *   -- NumOfConfigurations current count of the configurations
 *   -- config0000 current configuration 
 *   -- config#### configuration for each block
 * - random_state
 * Other classes can add datasets as needed.
 * open/close functions are provided so that the hdf5 file is closed during
 * the life time of this object. This is necessary so that failures do not lead
 * to unclosed hdf5.
 */
HDFWalkerOutput::HDFWalkerOutput(MCWalkerConfiguration& W,
    const string& aroot,Communicate* c): 
  NumOfWalkers(0), h_file(-1), h_config(-1), myComm(c)
{
  h5FileName=aroot+".config.h5";
  h_file = H5Fcreate(h5FileName.c_str(),H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);

  //create space in advance
  h_config = H5Gcreate(h_file,"config_collection",0); 

  int c0(0); 
  HDFAttribIO<int> i(c0);
  i.write(h_config,"NumOfConfigurations");

  //create config0000 
  char GrpName[128];
  sprintf(GrpName,"config%04d",0);
  hid_t group_id = H5Gcreate(h_config,GrpName,0);
  WalkerPackedWriter wo(W);
  wo.write(group_id,"coord");
  H5Gclose(h_config);

  hid_t h_random = H5Gcreate(h_file,"random_state",0);
  HDFAttribIO<RandomGenerator_t> r(Random);
  r.write(h_random,"dummy");
  H5Gclose(h_random);
  H5Fclose(h_file);
}

void HDFWalkerOutput::open()
{
  h_file =  H5Fopen(h5FileName.c_str(),H5F_ACC_RDWR,H5P_DEFAULT);
  h_config = H5Gopen(h_file,"config_collection");
}

/** Destructor writes the state of random numbers and close the file */
HDFWalkerOutput::~HDFWalkerOutput() {
  if(h_file<0)
    h_file =  H5Fopen(h5FileName.c_str(),H5F_ACC_RDWR,H5P_DEFAULT);
  hid_t h_random = H5Gopen(h_file,"random_state");
  HDFAttribIO<RandomGenerator_t> r(Random);
  r.write(h_random,"dummy");
  H5Fclose(h_file);
}

/** Write the set of walker configurations to the HDF5 file.  
 * @param W set of walker configurations
 */
bool HDFWalkerOutput::dump(MCWalkerConfiguration& W) {

  const bool replace=true;
  char GrpName[128];
  sprintf(GrpName,"config%04d",0);

  hid_t group_id = H5Gopen(h_config,GrpName);
  if(NumOfWalkers != W.getActiveWalkers())//need resize of the data
  {
    herr_t status=H5Gunlink(group_id,"coord");
    WalkerPackedWriter wo(W);
    wo.write(group_id,"coord");
    NumOfWalkers=W.getActiveWalkers();
  } 
  else
  {//simply overwrite it
    WalkerPackedWriter wo(W,replace);
    wo.write(group_id,"coord");
  }
  H5Gclose(group_id);

  int c(1); 
  HDFAttribIO<int> i(c,replace);
  i.write(h_config,"NumOfConfigurations");

  return true;
}

/** Write the set of walker configurations to the HDF5 file.  
 * @param W set of walker configurations
 * @param block index for the block
 */
bool  HDFWalkerOutput::append(MCWalkerConfiguration& W, int block) {

  if(block==0) {
    return dump(W);
  }

  char GrpName[128];
  sprintf(GrpName,"config%04d",block);
  hid_t group_id = H5Gcreate(h_config,GrpName,0);
  WalkerPackedWriter wo(W);
  wo.write(group_id,"coord");
  H5Gclose(group_id);

  block++;
  
  HDFAttribIO<int> i(block,true);
  i.write(h_config,"NumOfConfigurations");
  return true;
}
#endif
///** Write the set of walker configurations to the HDF5 file.  
// *@param W set of walker configurations
// */
//bool HDFWalkerOutput::get(MCWalkerConfiguration& W) {
//
//  typedef MCWalkerConfiguration::PosType PosType;
//  typedef MCWalkerConfiguration::RealType RealType;
//  typedef MCWalkerConfiguration::PropertyContainer_t PropertyContainer_t;
//
//  typedef Matrix<PosType>  PosContainer_t;
//  typedef Vector<RealType> ScalarContainer_t;
//
//  PropertyContainer_t Properties;
//  PosContainer_t tempPos(W.getActiveWalkers(),W.R.size());
//
//  //store walkers in a temporary array
//  int nw(0),item(0);
//  MCWalkerConfiguration::const_iterator it(W.begin());
//  MCWalkerConfiguration::const_iterator it_end(W.end());
//  while(it != it_end) {
//    for(int np=0; np < W.getParticleNum(); ++np) 
//      tempPos(item++) = (*it)->R(np);    
//    ++it; ++nw;
//  }
//
//  //create the group and increment counter
//  char GrpName[128];
//  sprintf(GrpName,"config%04d",Counter++);
//
//  herr_t status = H5Eset_auto(NULL, NULL);
//  status = H5Gget_objinfo (h_config, GrpName, 0, NULL);
//  if(status ==0)
//  {
//    hid_t group_id = H5Gopen(h_config,GrpName);
//    //write the dataset
//    HDFAttribIO<PosContainer_t> Pos_out(tempPos,true);
//    Pos_out.write(group_id,"coord");
//    H5Gclose(group_id);
//  }
//  else
//  {
//    hid_t group_id = H5Gcreate(h_config,GrpName,0);
//    //write the dataset
//    HDFAttribIO<PosContainer_t> Pos_out(tempPos);
//    Pos_out.write(group_id,"coord");
//    H5Gclose(group_id);
//  }
//  return true;
//}

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

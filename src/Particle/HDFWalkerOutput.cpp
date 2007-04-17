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

  inline void read(hid_t grp, const char* name) {
  }
};

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
    const string& aroot): NumOfWalkers(0), h_file(-1), h_config(-1)
{
  h5FileName=aroot+".config.h5";
  h_file = H5Fcreate(h5FileName.c_str(),H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);

  //create space in advance
  h_config = H5Gcreate(h_file,"config_collection",0); 

  int c(0); 
  HDFAttribIO<int> i(c);
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

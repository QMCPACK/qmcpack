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
#if defined(HAVE_LIBHDF5)
#include "Particle/HDFParticleAttrib.h"
#include "Numerics/HDFNumericAttrib.h"
#endif
#include "Utilities/OhmmsInfo.h"
#include "Utilities/RandomGenerator.h"
#include "OhmmsData/FileUtility.h"
using namespace qmcplusplus;

#if defined(HAVE_LIBHDF5)
/** Create the HDF5 file "aroot.config.h5" for output. 
 *@param aroot the root file name
 *@param append 
 *
 * @if append == true 
 * the configuration is appended
 * @else
 * the configuration is overwritten
 * The constructor
 * - opens a hdf5 file
 * - create the main group "/config_collection"
 */

HDFWalkerOutput::HDFWalkerOutput(const string& aroot, bool append, int count)
  : Counter(count) {

  string h5file(aroot);
  h5file.append(".config.h5");

  AppendMode=append;
  if(AppendMode)  {
    h_file =  H5Fopen(h5file.c_str(),H5F_ACC_RDWR,H5P_DEFAULT);
    h_config = H5Gopen(h_file,"config_collection");
    h_random = H5Gopen(h_file,"random_state");
    hid_t h1=H5Dopen(h_config,"NumOfConfigurations");
    if(h1>-1) {
      H5Dread(h1, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&(Counter));
      H5Dclose(h1);
    }
  }
  else {
    h_file = H5Fcreate(h5file.c_str(),H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
    h_config = H5Gcreate(h_file,"config_collection",0);
    h_random = H5Gcreate(h_file,"random_state",0);
  }
}

/** Destructor closes the HDF5 file and main group. */

HDFWalkerOutput::~HDFWalkerOutput() {

  hsize_t dim = 1;
  hid_t dataspace, dataset;

  dataspace= H5Screate_simple(1, &dim, NULL);
  if(AppendMode) {
    dataset= H5Dopen(h_config, "NumOfConfigurations");
  } else {
    dataset= H5Dcreate(h_config, "NumOfConfigurations", H5T_NATIVE_INT, dataspace, H5P_DEFAULT);
  }
  hid_t ret = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&Counter);
  H5Sclose(dataspace);
  H5Dclose(dataset);

  H5Gclose(h_config);
  Random.write(h_random,AppendMode);
  H5Gclose(h_random);
  H5Fclose(h_file);
}


/** Write the set of walker configurations to the HDF5 file.  
 *@param W set of walker configurations
 */
bool HDFWalkerOutput::get(MCWalkerConfiguration& W) {

  typedef MCWalkerConfiguration::PosType PosType;
  typedef MCWalkerConfiguration::RealType RealType;
  typedef MCWalkerConfiguration::PropertyContainer_t PropertyContainer_t;

  typedef Matrix<PosType>  PosContainer_t;
  typedef Vector<RealType> ScalarContainer_t;

  PropertyContainer_t Properties;
  PosContainer_t tempPos(W.getActiveWalkers(),W.R.size());

  //store walkers in a temporary array
  int nw(0),item(0);
  MCWalkerConfiguration::const_iterator it(W.begin());
  MCWalkerConfiguration::const_iterator it_end(W.end());
  while(it != it_end) {
    for(int np=0; np < W.getParticleNum(); ++np) 
      tempPos(item++) = (*it)->R(np);    
    ++it; ++nw;
  }

  //create the group and increment counter
  char GrpName[128];
  sprintf(GrpName,"config%04d",Counter++);
  hid_t group_id = H5Gcreate(h_config,GrpName,0);

  //write the dataset
  HDFAttribIO<PosContainer_t> Pos_out(tempPos);
  Pos_out.write(group_id,"coord");

  H5Gclose(group_id);
  
  return true;
}
#else
HDFWalkerOutput::HDFWalkerOutput(const string& aroot, bool append, int count)
  : Counter(count) {
}

/** Destructor closes the HDF5 file and main group. */

HDFWalkerOutput::~HDFWalkerOutput() {
}


/** Write the set of walker configurations to the HDF5 file.  
 *@param W set of walker configurations
 */
bool HDFWalkerOutput::get(MCWalkerConfiguration& W) {
  return false;
}
#endif

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

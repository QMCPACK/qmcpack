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
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//   Department of Physics, Ohio State University
//   Ohio Supercomputer Center
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#include "Particle/MCWalkerConfiguration.h"
#include "Particle/HDFWalkerIO.h"
#include "OhmmsPETE/OhmmsVector.h"
#include "Particle/HDFParticleAttrib.h"
#include "Numerics/HDFNumericAttrib.h"
#include "Utilities/OhmmsInfo.h"
#include "Utilities/RandomGenerator.h"
#include "OhmmsData/FileUtility.h"
using namespace ohmmsqmc;

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
  : Counter(count), AppendMode(append) {

  string h5file(aroot);
  h5file.append(".config.h5");
  if(AppendMode)  {
    h_file =  H5Fopen(h5file.c_str(),H5F_ACC_RDWR,H5P_DEFAULT);
    h_config = H5Gopen(h_file,"config_collection");
    h_random = H5Gopen(h_file,"random_state");
  }
  else {
    Counter=0;
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

/** Open the HDF5 file "aroot.config.h5" for reading. 
 *@param aroot the root file name
 *
 *@note The main group is "/config_collection"
 */
HDFWalkerInput::HDFWalkerInput(const string& aroot, int ipart, int nparts):
Counter(0), NumSets(0) {
  string h5file = aroot;

  string ext=getExtension(h5file);
  if(ext != "h5") { //if the filename does not h5 extension, add the extension
    h5file.append(".config.h5");
  }

  h_file =  H5Fopen(h5file.c_str(),H5F_ACC_RDWR,H5P_DEFAULT);
  h_config = H5Gopen(h_file,"config_collection");
  hid_t h1=H5Dopen(h_config,"NumOfConfigurations");
  if(h1>-1) {
    H5Dread(h1, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&(NumSets));
    H5Dclose(h1);
  }

  if(!NumSets) {
    H5Gget_num_objs(h_config,&NumSets);
  }

  if(NumSets) {
    int nframe=NumSets/nparts;
    FirstSet=nframe*ipart;
    LastSet=nframe*(ipart+1);
  } else {
    FirstSet=0; LastSet=0;
    ERRORMSG("File does not contain walkers!")
  }
}

/** Destructor closes the HDF5 file and main group. */
HDFWalkerInput::~HDFWalkerInput() {
  H5Gclose(h_config);
  H5Fclose(h_file);
}

/** Read a configuration from a hdf5 file.
 *@param W set of walker configurations
int HDFWalkerInput::put(MCWalkerConfiguration& W){

  //no more walker to read
  if(Counter >= NumSets) return 0;

  //current Counter 
  int ic = Counter++;

  //This configuration may be skipped
  if(ic<FirstSet || ic>=LastSet) return 2;

  put(W,ic);
  return 1;
}
 */

/**  Write the set of walker configurations to the HDF5 file.  
 *@param W set of walker configurations
 *@param ic the number of frames
 *
 * \if ic==-1
 *  use only the last frame for a restart
 * \else if ic>=0
 *  use ic frames from the file for opitimizations
 */
bool  
HDFWalkerInput::put(MCWalkerConfiguration& W, int ic){

  if(Counter<0) return false;

  int selected = ic;
  if(ic<0) {
    XMLReport("Will use the last set from " << NumSets << " of configurations.")
    selected = NumSets-1;
  }

  typedef MCWalkerConfiguration::PosType PosType;
  typedef MCWalkerConfiguration::PropertyContainer_t ProtertyContainer_t;

  typedef Matrix<PosType>  PosContainer_t;

  int nwt = 0;
  int npt = 0;
  //2D array of PosTypes (x,y,z) indexed by (walker,particle)
  PosContainer_t Pos_temp;

  //open the group
  char GrpName[128];
  sprintf(GrpName,"config%04d",selected);
  hid_t group_id = H5Gopen(h_config,GrpName);
    
  HDFAttribIO<PosContainer_t> Pos_in(Pos_temp);
  //read the dataset
  Pos_in.read(group_id,"coord");
  //close the group
  H5Gclose(group_id);

  /*check to see if the number of walkers and particles is  consistent with W */
  int nptcl = Pos_temp.cols();
  nwt = Pos_temp.rows();

  int curWalker = W.getActiveWalkers();
  if(curWalker) {
    LOGMSG("Adding " << nwt << " walkers to " << curWalker)
    W.createWalkers(nwt);
  } else {
    W.resize(nwt,nptcl); 
  }

  //assign configurations to W
  int iw=0;
  MCWalkerConfiguration::iterator it = W.begin()+curWalker; 
  MCWalkerConfiguration::iterator it_end = W.end(); 
  while(it != it_end) {
    std::copy(Pos_temp[iw],Pos_temp[iw+1], (*it)->R.begin());
    ++it;++iw;
  }

  return true;
}

bool HDFWalkerInput::append(MCWalkerConfiguration& W){

  if(Counter<0) return false;

  typedef MCWalkerConfiguration::PosType PosType;
  typedef Matrix<PosType>  PosContainer_t;
  PosContainer_t Pos_temp;
  int nw_in=0,curConfig=FirstSet;
  while(curConfig<LastSet) {
    //open the group
    char GrpName[128];
    sprintf(GrpName,"config%04d",curConfig);
    hid_t group_id = H5Gopen(h_config,GrpName);
    HDFAttribIO<PosContainer_t> Pos_in(Pos_temp);
    //read the dataset
    Pos_in.read(group_id,"coord");
    //close the group
    H5Gclose(group_id);
    /*check to see if the number of walkers and particles is  consistent with W */
    int nptcl = Pos_temp.cols();
    int nwt = Pos_temp.rows();
    int curWalker = W.getActiveWalkers();
    if(curWalker) {
      W.createWalkers(nwt);
    } else {
      W.resize(nwt,nptcl); 
    }
    MCWalkerConfiguration::iterator it = W.begin()+curWalker; 
    for(int iw=0; iw<nwt; iw++) {
      //std::copy(Post_temp[iw],Post_temp[iw+1], (*it)->R.begin());
      for(int iat=0; iat < nptcl; iat++){
        (*it)->R(iat) = Pos_temp(iw,iat);
      }
      ++it;
    }
    nw_in += nwt; 
    curConfig++;
  }

  LOGMSG("Total " << nw_in << " walkers are loaded using " << LastSet-FirstSet << " frames.")
  return true;
}

bool  
HDFWalkerInput::append(MCWalkerConfiguration& W, int nwalkers){

  if(Counter<0) return false;

  if(nwalkers<0) return put(W,-1);

  typedef MCWalkerConfiguration::PosType PosType;
  typedef Matrix<PosType>  PosContainer_t;
  PosContainer_t Pos_temp;

  int nw_in=0,curConfig=NumSets-1, numConfigIn=0;
  while(curConfig>=0 && nw_in<nwalkers){
    //open the group
    char GrpName[128];
    sprintf(GrpName,"config%04d",curConfig);
    hid_t group_id = H5Gopen(h_config,GrpName);
    HDFAttribIO<PosContainer_t> Pos_in(Pos_temp);
    //read the dataset
    Pos_in.read(group_id,"coord");
    //close the group
    H5Gclose(group_id);
    /*check to see if the number of walkers and particles is  consistent with W */
    int nptcl = Pos_temp.cols();
    int nwt = Pos_temp.rows();
    int curWalker=0;
    if(nptcl != W.getParticleNum()) {
      W.resize(nwt,nptcl); 
    } else {
      curWalker=W.getActiveWalkers();
      W.createWalkers(nwt);
    }
    MCWalkerConfiguration::iterator it = W.begin()+curWalker; 
    for(int iw=0; iw<nwt; iw++) {
      //std::copy(Post_temp[iw],Post_temp[iw+1], (*it)->R.begin());
      for(int iat=0; iat < nptcl; iat++){
        (*it)->R(iat) = Pos_temp(iw,iat);
      }
      ++it;
    }
    nw_in += nwt; 
    curConfig--;
    numConfigIn++;
  }

  LOGMSG("Total " << nw_in << " walkers are loaded using " << numConfigIn << " frames.")
  return true;
}

void  HDFWalkerInput::getRandomState(bool restart){
  if(restart) {
    hid_t h_random = H5Gopen(h_file,"random_state");
    if(h_random>-1) {
      LOGMSG("Reading the state of the random number generator from the configuration file")
      Random.read(h_random);
      H5Gclose(h_random);
    }
  }
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

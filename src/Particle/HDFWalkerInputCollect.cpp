//////////////////////////////////////////////////////////////////
// (c) Copyright 2005- by Jeongnim Kim
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
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#include "Particle/MCWalkerConfiguration.h"
#include "Particle/HDFWalkerInputCollect.h"
#include "OhmmsPETE/OhmmsVector.h"
#include "Particle/HDFParticleAttrib.h"
#include "Numerics/HDFNumericAttrib.h"
#include "Utilities/OhmmsInfo.h"
#include "Utilities/RandomGenerator.h"
#include "OhmmsData/FileUtility.h"
using namespace qmcplusplus;


/** Open the HDF5 file "aroot.config.h5" for reading. 
 *@param aroot the root file name
 *
 *@note The main group is "/config_collection"
 */
HDFWalkerInputCollect::HDFWalkerInputCollect(const string& aroot): 
 fileID(-1), prevNContexts(1), curNContexts(1) {
  open(aroot);
}

HDFWalkerInputCollect::~HDFWalkerInputCollect() {
  close();
}

bool 
HDFWalkerInputCollect::open(const string& aroot) {

  string h5file = aroot;

  string ext=getExtension(h5file);
  if(ext != "h5") { //if the filename does not h5 extension, add the extension
    h5file.append(".config.h5");
  }

  fileID =  H5Fopen(h5file.c_str(),H5F_ACC_RDWR,H5P_DEFAULT);
  if(fileID<0) {
    app_error() << "Cannot open hdf5 file " << endl;
    return false;
  }

  app_log() << "  HDFWalkerInputCollect::open " << h5file<< endl;
  curNContexts = OHMMS::Controller->ncontexts();
  OffSet.resize(curNContexts+1);

  return true;
}

bool HDFWalkerInputCollect::close() {

  if(fileID >-1) {
    prevNContexts=1;
    herr_t status = H5Eset_auto(NULL, NULL);
    hid_t h1=H5Dopen(fileID,"ncontexts");
    if(h1>-1) {
      H5Dread(h1, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(prevNContexts));
      H5Dclose(h1);
      char rname[128];
      if(prevNContexts == curNContexts) {
        app_log() << "    Restart with the random states" << endl;
        sprintf(rname,"context%04d/random_state",OHMMS::Controller->mycontext());
        hid_t ranIn = H5Gopen(fileID,rname);
        Random.read(ranIn);
        H5Gclose(ranIn);
      } else {
        app_warning() << "The number of processors has changed.\n"
          << "New random seeds are generated." << endl;
      }
    } else {
      if(curNContexts == 1) {
        hid_t h_random = H5Gopen(fileID,"random_state");
        if(h_random>-1) {
          app_log() << "Using serial random seed" << endl;
          Random.read(h_random);
          H5Gclose(h_random);
        }
      }
    }
    app_log() << "  HDFWalkerInputCollect::close " << endl;
    //read the random seeds
    H5Fclose(fileID);
    fileID=-1;
  }
  return true;
}

bool 
HDFWalkerInputCollect::put(MCWalkerConfiguration& W, int rollback) {

  int myID = OHMMS::Controller->mycontext();

  int numConfigs=1;
  char confName[128];
  char coordName[128];


  hid_t mastercf = H5Gopen(fileID,"config_collection");
  app_log() << "    The number of blocks in the file " << numConfigs << endl;
  hid_t h1=H5Dopen(mastercf,"NumOfConfigurations");

  if(h1>-1) {
    H5Dread(h1, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&(numConfigs));
    H5Dclose(h1);
  }

  int firstConf(0);
  if(rollback>0) {
    if(rollback<numConfigs) 
      firstConf=numConfigs-rollback;
    else
      firstConf = 0;
  }

  hsize_t dimTot[3], dimIn[3], offset[3];
  offset[0]=0; offset[1]=0; offset[2]=0;

  typedef MCWalkerConfiguration::PosType PosType;
  
  vector<PosType> pos;
  int nwRead=0;

  for(int iconf=firstConf; iconf<numConfigs; iconf++) { 

    sprintf(coordName,"config%04d/coord",iconf); 
    hid_t dataset = H5Dopen(mastercf,coordName);
    hid_t dataspace = H5Dget_space(dataset);
    int rank = H5Sget_simple_extent_ndims(dataspace);
    int status_n = H5Sget_simple_extent_dims(dataspace, dimTot, NULL);

    distribute(dimTot[0]);

    app_log() << "    Distibution of the walkers: ";
    std::copy(OffSet.begin(), OffSet.end(), ostream_iterator<hsize_t>(app_log(), " "));
    app_log() << endl;

    //get the input dimension
    dimIn[0]=OffSet[myID+1]-OffSet[myID];
    dimIn[1]=dimTot[1]; dimIn[2]=dimTot[2];
    offset[0]=OffSet[myID];

    vector<PosType> posIn(dimIn[0]*dimIn[1]);

    hid_t memspace = H5Screate_simple(3, dimIn, NULL);
    herr_t status = H5Sselect_hyperslab(dataspace,H5S_SELECT_SET, offset,NULL,dimIn,NULL);
    status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace, H5P_DEFAULT, &(posIn[0][0]));

    H5Sclose(memspace);
    H5Dclose(dataset);
    H5Sclose(dataspace);

    pos.insert(pos.end(), posIn.begin(), posIn.end());
    nwRead += dimIn[0];
  }

  H5Gclose(mastercf);

  int curWalker = W.getActiveWalkers();
  int nptcl=W.getTotalNum();
  if(curWalker) {
    W.createWalkers(nwRead);
  } else {
    W.resize(nwRead,nptcl);
  }

  MCWalkerConfiguration::iterator it = W.begin()+curWalker; 
  int ii=0;
  for(int iw=0; iw<nwRead; iw++) {
    //std::copy(Post_temp[iw],Post_temp[iw+1], (*it)->R.begin());
    for(int iat=0; iat < nptcl; iat++,ii++){
      (*it)->R(iat) = pos[ii];
    }
    ++it;
  }

  return true;
}

void HDFWalkerInputCollect::distribute(hsize_t nw) {
  hsize_t bat=nw/curNContexts;
  hsize_t residue = curNContexts-nw%curNContexts;
  OffSet[0] = 0;
  for(int i=0; i<curNContexts; i++) {
    if(i<residue)
      OffSet[i+1] = OffSet[i] + bat;
    else
      OffSet[i+1] = OffSet[i] + bat+1;
  }
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

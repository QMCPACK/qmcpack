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
/** @file HDFWalkerOuput.cpp
 * @breif definition  of HDFWalkerOuput  and other support class and functions
 */
#include "Particle/MCWalkerConfiguration.h"
#include "Particle/HDFWalkerOutput.h"
#include "OhmmsPETE/OhmmsVector.h"
#include "Utilities/OhmmsInfo.h"
#include "Utilities/RandomGenerator.h"
#include "Utilities/RandomGeneratorIO.h"
#include "OhmmsData/FileUtility.h"
#include "Message/Communicate.h"
#include "Particle/HDFWalkerIOEngine.h"

namespace qmcplusplus 
{
  /** utility function to add a configuration to a hdf5
   * @param h0 node to which a new block is added
   * @param W current walker configuration
   * @param block current block
   * @return block+1
   */
  int appendSingle(hid_t h0, MCWalkerConfiguration& W, int block)
  {
    hid_t h1=-1;
    bool overwrite=true;
    if(block) 
      h1 = H5Gopen(h0,hdf::config_group);
    else
    {
      overwrite=false;
      h1 = H5Gcreate(h0,hdf::config_group,0);
    }

    char newblock[16];
    sprintf(newblock,"%s%04d",hdf::append_walkers,block);
    HDFWalkerIOEngine wo(W);
    wo.write(h1,newblock);

    block++;
    HDFAttribIO<int> i(block,overwrite);
    i.write(h1,hdf::num_blocks);

    H5Gclose(h1);
    return block;
  }

  /* constructor
   * @param W walkers to operate on
   * @param aroot the root file name
   * @param c communicator
   *
   * Create the HDF5 file "aroot.config.h5" for output. 
   * Constructor is responsible to create a hdf5 file and create basic groups
   * so that subsequent write can utilize existing dataspaces.
   * HDF5 contains
   * - state_0
   *   -- walkers
   * - config_collection
   *   -- NumOfConfigurations current count of the configurations
   *   -- config0000 current configuration 
   *   -- config#### configuration for each block
   * Other classes can add datasets as needed.
   * open/close functions are provided so that the hdf5 file is closed during
   * the life time of this object. This is necessary so that failures do not lead
   * to unclosed hdf5.
   */
  HDFWalkerOutput::HDFWalkerOutput(MCWalkerConfiguration& W, const string& aroot,Communicate* c): 
    appended_blocks(0), number_of_walkers(0), 
  number_of_backups(0), max_number_of_backups(4),
  h_file(-1), h_plist(H5P_DEFAULT), xfer_plist(H5P_DEFAULT), h_state(-1), myComm(c)
  {
    FileName=myComm->getName()+hdf::config_ext;
#if defined(H5_HAVE_PARALLEL)
    if(myComm->size()>1)
    { //open file collectively
      MPI_Info info=MPI_INFO_NULL;
      h_plist = H5Pcreate(H5P_FILE_ACCESS);
      H5Pset_fapl_mpio(h_plist,myComm->getMPI(),info);
      h_file = H5Fcreate(FileName.c_str(),H5F_ACC_TRUNC,H5P_DEFAULT,h_plist);

      //use_collective
      xfer_plist = H5Pcreate(H5P_DATASET_XFER);
      H5Pset_dxpl_mpio(xfer_plist,H5FD_MPIO_COLLECTIVE);
      //H5Pset_dxpl_mpio(xfer_plist,H5FD_MPIO_INDEPENDENT);

      //write the version
      HDFVersion cur_version;
      cur_version.write(h_file,hdf::version);

      //create a state
      h_state = H5Gcreate(h_file,hdf::main_state,0);
      HDFWalkerIOEngine wo(W);
      wo.writeAll(h_state,hdf::walkers,myComm);

      number_of_walkers = W.getGlobalNumWalkers();
      ////for debugging only. Open a file per node
      //string fname=aroot+".config.h5";
      //std::pair<hid_t,hid_t> h_pair= createH5FileSingle(fname,W);
      //h_debug_file=h_pair.first;
      //H5Gclose(h_pair.second);
    }
    else
    {
      std::pair<hid_t,hid_t> h_pair= createH5FileSingle(FileName,W);
      number_of_walkers = W.getActiveWalkers();
      h_file=h_pair.first;
      h_state=h_pair.second;
    }
#else
    if(myComm->size()>1)
    {
      char fname[128];
      sprintf(fname,"%s.p%03d%s",myComm->getName().c_str(), myComm->rank(),hdf::config_ext);
      FileName=fname;
    }
    std::pair<hid_t,hid_t> h_pair= createH5FileSingle(FileName,W);
    number_of_walkers = W.getActiveWalkers();
    h_file=h_pair.first;
    h_state=h_pair.second;
#endif

    HDFAttribIO<int> nwo(number_of_walkers);
    nwo.write(h_state,hdf::num_walkers);

    H5Gclose(h_state);
    H5Fclose(h_file);
  }

  /** Write the set of walker configurations to the HDF5 file.  
   * @param W set of walker configurations
   */
  bool HDFWalkerOutput::dump(MCWalkerConfiguration& W) {
    bool success=true;
    h_file =  H5Fopen(FileName.c_str(),H5F_ACC_RDWR,h_plist);
    h_state =  H5Gopen(h_file,hdf::main_state);

    //backup problem with phdf5
    //vector<string> backups;
    //backups.push_back(hdf::main_state);
    //for(int i=1;i<=max_number_of_backups; i++)
    //{
    //  ostringstream o; o<<"state_"<<i;
    //  backups.push_back(o.str());
    //}
    //herr_t status;
    //if(number_of_backups==max_number_of_backups)
    //{
    //  status=H5Gunlink(h_file,backups[max_number_of_backups].c_str());
    //}
    //int b=number_of_backups;
    //while(b)
    //{
    //  status=H5Gmove(h_file,backups[b-1].c_str(),backups[b].c_str());
    //  b--;
    //}
    //if(number_of_backups<max_number_of_backups) number_of_backups++;

    if(h_plist== H5P_DEFAULT) // serial file
    { 
      success= dumpSingle(h_state,W);
      //HDFWalkerIOEngine wo(W);
      //wo.write(h_state,hdf::walkers);
      //number_of_walkers=W.getActiveWalkers();
      //HDFAttribIO<int> nwo(number_of_walkers);
      //nwo.write(h_state,hdf::num_walkers);
    }
    else
    {//parallel file
      if(number_of_walkers != W.getGlobalNumWalkers())
      {
        herr_t status=H5Gunlink(h_state,hdf::walkers);
        HDFWalkerIOEngine wo(W);
        wo.setTransferProperty(xfer_plist);
        wo.writeAll(h_state,hdf::walkers,myComm);

        // overwrite number of walkers
        number_of_walkers=W.getGlobalNumWalkers();
        HDFAttribIO<int> nwo(number_of_walkers);
        nwo.write(h_state,hdf::num_walkers);
      }
      else
      {
        HDFWalkerIOEngine wo(W,true);
        wo.setTransferProperty(xfer_plist);
        wo.writeAll(h_state,hdf::walkers,myComm);
      }
    }
    H5Gclose(h_state);
    H5Fclose(h_file);

    //appended_blocks=appendSingle(h_debug_file,W,appended_blocks);
    return success;
  }

  /** Write the set of walker configurations to the HDF5 file.  
   * @param W set of walker configurations
   * @param block index for the block
   */
  bool  HDFWalkerOutput::append(MCWalkerConfiguration& W) {
    hid_t h0=h_file;
    if(h0<0) h0 =  H5Fopen(FileName.c_str(),H5F_ACC_RDWR,H5P_DEFAULT);
    appended_blocks=appendSingle(h0,W,appended_blocks);
    if(h0 != h_file) H5Fclose(h0);
    return true;
  }


  std::pair<hid_t,hid_t> 
    HDFWalkerOutput::createH5FileSingle(const string& fname, MCWalkerConfiguration& W)
    {
      hid_t fid= H5Fcreate(fname.c_str(),H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);

      HDFVersion cur_version;
      cur_version.write(fid,hdf::version);

      //create space in advance 
      hid_t h1 = H5Gcreate(fid,hdf::main_state,0); 

      HDFWalkerIOEngine wo(W);
      wo.write(h1,hdf::walkers);

      return std::pair<hid_t,hid_t>(fid,h1);
    }

  /** Write the set of walker configurations to the HDF5 file.  
   * @param W set of walker configurations
   */
  bool HDFWalkerOutput::dumpSingle(hid_t gid, MCWalkerConfiguration& W) {
    int nw=W.getActiveWalkers();
    if(number_of_walkers != nw)
    {//need resize of the data: unlink and create a walkers group
      herr_t status=H5Gunlink(gid,hdf::walkers);
      HDFWalkerIOEngine wo(W);
      wo.write(gid,hdf::walkers);
      number_of_walkers=nw;
      HDFAttribIO<int> nwo(number_of_walkers);
      nwo.write(h_state,hdf::num_walkers);
    } 
    else
    {//simply overwrite it
      const bool replace=true;
      HDFWalkerIOEngine wo(W,replace);
      wo.write(h_state,hdf::walkers);
    }
    return true;
  }

  /** Destructor writes the state of random numbers and close the file */
  HDFWalkerOutput::~HDFWalkerOutput() 
  {
    if(xfer_plist!= H5P_DEFAULT) H5Pclose(xfer_plist);
    if(h_plist!= H5P_DEFAULT) H5Pclose(h_plist);
  }

  void HDFWalkerOutput::open()
  {
    h_file =  H5Fopen(FileName.c_str(),H5F_ACC_RDWR,h_plist);
    h_state = H5Gopen(h_file,hdf::main_state);
  }
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

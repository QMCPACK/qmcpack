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
#include "Particle/HDFWalkerOutput.h"
#include "Utilities/IteratorUtility.h"
#include "OhmmsData/FileUtility.h"
#include "HDFVersion.h"
#include <numeric>
#include <iostream>
#include <sstream>
#include <Message/Communicate.h>
#include <mpi/collectives.h>
#include <io/hdf_hyperslab.h>

namespace qmcplusplus
{
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
  HDFWalkerOutput::HDFWalkerOutput(MCWalkerConfiguration& W, const string& aroot,Communicate* c)
    : appended_blocks(0), number_of_walkers(0), currentConfigNumber(0)
      , number_of_backups(0), max_number_of_backups(4), myComm(c), RootName(aroot)
      , fw_out(myComm)
  {
    number_of_particles=W.getTotalNum();
    //FileName=myComm->getName()+hdf::config_ext;
    //ConfigFileName=myComm->getName()+".storeConfig.h5";
    string ConfigFileName=myComm->getName()+".storeConfig.h5";
    HDFVersion cur_version;
    int dim=OHMMS_DIM;
    fw_out.create(ConfigFileName);
    fw_out.write(cur_version.version,hdf::version);
    fw_out.write(number_of_particles,"NumberElectrons");
    fw_out.write(dim,"DIM");
  }

  /** Destructor writes the state of random numbers and close the file */
  HDFWalkerOutput::~HDFWalkerOutput()
  {
    fw_out.close();
    delete_iter(RemoteData.begin(),RemoteData.end());
    delete_iter(FWData.begin(),FWData.end());
  }

  /** Write the set of walker configurations to the HDF5 file.
   * @param W set of walker configurations
   *
   * Dump is for restart and do not preserve the file
   */
  bool HDFWalkerOutput::dump(MCWalkerConfiguration& W) 
  {
    string FileName=myComm->getName()+hdf::config_ext;
    hdf_archive dump_file(myComm,true);
    dump_file.create(FileName);
    HDFVersion cur_version;
    //version
    dump_file.write(cur_version.version,hdf::version);
    //state
    dump_file.push(hdf::main_state);
    //walkers
    write_configuration(W,dump_file);
    dump_file.close();
    return true;
  }

  void HDFWalkerOutput::write_configuration(MCWalkerConfiguration& W, hdf_archive& hout)
  {
    if (RemoteData.empty())
    {//add two buffers
      RemoteData.push_back(new BufferType);
      RemoteData.push_back(new BufferType);
    }
    const int wb=OHMMS_DIM*number_of_particles;

    //populate RemoteData[0] to dump
#if defined(H5_HAVE_PARALLEL) && defined(ENABLE_PHDF5)
    RemoteData[0]->resize(wb*W.getActiveWalkers());
    W.putConfigurations(RemoteData[0]->begin());

    //TinyVector<hsize_t,3> gcounts, counts, offset;
    hsize_t gcounts[3], counts[3], offset[3];
    gcounts[0]=W.WalkerOffsets[myComm->size()]; gcounts[1]=number_of_particles; gcounts[2]=OHMMS_DIM;
    counts[0]=W.getActiveWalkers(); counts[1]=number_of_particles; counts[2]=OHMMS_DIM;
    offset[0]=W.WalkerOffsets[myComm->rank()]; offset[1]=0; offset[2]=0;

    BufferType::value_type t;
    const hid_t etype=get_h5_datatype(t);
    hout.write(gcounts[0],hdf::num_walkers);

    hid_t gid=hout.top();
    hid_t sid1  = H5Screate_simple(3,gcounts,NULL);
    hid_t memspace=H5Screate_simple(3,counts,NULL);
    hid_t dset_id=H5Dcreate(gid,hdf::walkers,etype,sid1,H5P_DEFAULT);
    hid_t filespace=H5Dget_space(dset_id);
    herr_t ret=H5Sselect_hyperslab(filespace,H5S_SELECT_SET,offset,NULL,counts,NULL);
    ret = H5Dwrite(dset_id,etype,memspace,filespace,hout.xfer_plist,RemoteData[0]->data());
    H5Sclose(filespace);
    H5Dclose(dset_id);

#else
    if(myComm->size()==1)
    {
      number_of_walkers=W.getActiveWalkers();
      RemoteData[0]->resize(wb*W.getActiveWalkers());
      W.putConfigurations(RemoteData[0]->begin());
    }
#if defined(HAVE_MPI)
    else
    {
      RemoteData[1]->resize(wb*W.getActiveWalkers());
      W.putConfigurations(RemoteData[1]->begin());
      vector<int> displ(myComm->size()), counts(myComm->size());
      for (int i=0; i<myComm->size(); ++i)
      {
        counts[i]=wb*(W.WalkerOffsets[i+1]-W.WalkerOffsets[i]);
        displ[i]=wb*W.WalkerOffsets[i];
      }

      if(!myComm->rank()) RemoteData[0]->resize(wb*W.WalkerOffsets.back());
      mpi::gatherv(*myComm,*RemoteData[1],*RemoteData[0],counts,displ);
      //myComm->gatherv(*RemoteData[0],*RemoteData[1],counts, displ);
      number_of_walkers=W.WalkerOffsets[myComm->size()];
    }
#endif
    vector<hsize_t> inds(3);
    inds[0]=number_of_walkers;
    inds[1]=number_of_particles;
    inds[2]=OHMMS_DIM;
    hout.write(number_of_walkers,hdf::num_walkers);

    hyperslab_proxy<BufferType,3> slab(*RemoteData[0],inds);
    hout.write(slab,hdf::walkers);
#endif
    //HDFAttribIO<BufferType> po(*RemoteData[0],inds);
    //po.write(hout.top(),hdf::walkers,hout.xfer_plist);
  }

bool HDFWalkerOutput::dump(ForwardWalkingHistoryObject& FWO)
{
    string ConfigFileName=myComm->getName()+".storeConfig.h5";
    fw_out.open(ConfigFileName);

    if (myComm->size()==1)
    {
      for (int i=0; i<FWO.ForwardWalkingHistory.size(); i++ )
      {
        int fwdata_size=FWO.ForwardWalkingHistory[i]->size();
        std::stringstream sstr;
        sstr<<"Block_"<<currentConfigNumber;
        fw_out.push(sstr.str());//open the group

        vector<float> posVecs;
        //reserve enough space
        vector<long> IDs(fwdata_size,0);
        vector<long> ParentIDs(fwdata_size,0);
        vector<float>::iterator tend(posVecs.begin());
        for (int j=0;j<fwdata_size;j++) 
        {
          const ForwardWalkingData& fwdata(FWO(i,j));
          IDs[j]=fwdata.ID;
          ParentIDs[j]=fwdata.ParentID;
          fwdata.append(posVecs);
        }
        fw_out.write(posVecs,"Positions");
        fw_out.write(IDs,"WalkerID");
        fw_out.write(ParentIDs,"ParentID");
        fw_out.write(fwdata_size,hdf::num_walkers);

        fw_out.pop();//close the group
        ++currentConfigNumber;
      }
    }
#if defined(HAVE_MPI)
    else
    {
      const int n3=number_of_particles*OHMMS_DIM;
      for (int i=0; i<FWO.ForwardWalkingHistory.size(); i++ )
      {
        int fwdata_size=FWO.ForwardWalkingHistory[i]->size();
        std::stringstream sstr;
        sstr<<"Block_"<<currentConfigNumber;
        fw_out.push(sstr.str());//open the group

        vector<int> counts(myComm->size());
        mpi::all_gather(*myComm,fwdata_size,counts);

        vector<float> posVecs;
        //reserve space to minimize the allocation
        posVecs.reserve(FWO.number_of_walkers*n3);
        vector<long> myIDs(fwdata_size),pIDs(fwdata_size);
        vector<float>::iterator tend(posVecs.begin());
        for (int j=0;j<fwdata_size;j++) 
        {
          const ForwardWalkingData& fwdata(FWO(i,j));
          myIDs[j]=fwdata.ID;
          pIDs[j]=fwdata.ParentID;
          fwdata.append(posVecs);
        }
        vector<int> offsets(myComm->size()+1,0);
        for(int i=0; i<myComm->size();++i) offsets[i+1]=offsets[i]+counts[i];
        fwdata_size=offsets.back();
        fw_out.write(fwdata_size,hdf::num_walkers);

        vector<long> globalIDs;
        if(myComm->rank()==0) globalIDs.resize(fwdata_size);

        //collect WalkerID
        mpi::gatherv(*myComm, myIDs, globalIDs, counts, offsets);
        fw_out.write(globalIDs,"WalkerID");
        //collect ParentID
        mpi::gatherv(*myComm, pIDs, globalIDs, counts, offsets);
        fw_out.write(globalIDs,"ParentID");

        for(int i=0; i<counts.size();++i) counts[i]*=n3;
        for(int i=0; i<offsets.size();++i) offsets[i]*=n3;
        vector<float> gpos;
        if(myComm->rank()==0) gpos.resize(offsets.back());
        mpi::gatherv(*myComm, posVecs, gpos, counts, offsets);

        fw_out.write(gpos,"Positions");

        fw_out.pop();//close the group
        ++currentConfigNumber;
      }
    }
#endif
    fw_out.close();
    FWO.clearConfigsForForwardWalking();
    return true;
}

}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/

//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: D.C. Yang, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Cynthia Gu, zg1@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



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
#if defined(HAVE_ADIOS) && defined(ADIOS_VERIFY)
#include <adios.h>
#include "ADIOS/ADIOS_verify.h"
#endif

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
HDFWalkerOutput::HDFWalkerOutput(MCWalkerConfiguration& W, const std::string& aroot,Communicate* c)
  : appended_blocks(0), number_of_walkers(0), currentConfigNumber(0)
  , number_of_backups(0), max_number_of_backups(4), myComm(c), RootName(aroot)
//       , fw_out(myComm)
{
  number_of_particles=W.getTotalNum();
  RemoteData.reserve(4);
  RemoteData.push_back(new BufferType);
  RemoteData.push_back(new BufferType);
  block = -1;
//     //FileName=myComm->getName()+hdf::config_ext;
//     //ConfigFileName=myComm->getName()+".storeConfig.h5";
//     std::string ConfigFileName=myComm->getName()+".storeConfig.h5";
//     HDFVersion cur_version;
//     int dim=OHMMS_DIM;
//     fw_out.create(ConfigFileName);
//     fw_out.write(cur_version.version,hdf::version);
//     fw_out.write(number_of_particles,"NumberElectrons");
//     fw_out.write(dim,"DIM");
}

/** Destructor writes the state of random numbers and close the file */
HDFWalkerOutput::~HDFWalkerOutput()
{
//     fw_out.close();
  delete_iter(RemoteData.begin(),RemoteData.end());
}

#ifdef HAVE_ADIOS
uint64_t HDFWalkerOutput::get_group_size(MCWalkerConfiguration& W)
{
  int walker_num =  W.getActiveWalkers();
  int particle_num = number_of_particles;
  return sizeof(int) * 4 + sizeof(OHMMS_PRECISION) * walker_num * particle_num * OHMMS_DIM;
}

bool HDFWalkerOutput::adios_checkpoint(MCWalkerConfiguration& W, int64_t adios_handle, int nblock)
{
  //Need 4 * 3 bytes for storing integers and then storage
  //for all the walkers
  if (RemoteData.empty())
    RemoteData.push_back(new BufferType);
  int walker_num =  W.getActiveWalkers();
  int particle_num = number_of_particles;
  int walker_dim_num = OHMMS_DIM;
  if(nblock > block)
  {
    //This is just another wrapper for vector
    RemoteData[0]->resize(walker_num * particle_num * walker_dim_num);
    //Copy over all the walkers into one chunk of contigous memory
    W.putConfigurations(RemoteData[0]->begin());
    block = nblock;
  }
  void* walkers = RemoteData[0]->data();
  uint64_t adios_groupsize, adios_totalsize;
  adios_groupsize = sizeof(int) * 3 + sizeof(*(RemoteData[0]->data()));
  adios_group_size (adios_handle, adios_groupsize, &adios_totalsize);
  adios_write (adios_handle, "walker_num", &walker_num);
  adios_write (adios_handle, "particle_num", &particle_num);
  int walker_size = walker_num * particle_num * walker_dim_num;
  adios_write (adios_handle, "walker_size", &walker_size);
  adios_write (adios_handle, "walkers", walkers);
  return true;
}

#ifdef ADIOS_VERIFY

void HDFWalkerOutput::adios_checkpoint_verify(MCWalkerConfiguration& W, ADIOS_FILE*fp)
{
  if (RemoteData.empty())
  {
    APP_ABORT_TRACE(__FILE__, __LINE__, "empty RemoteData. Not possible");
  }

  //RemoteData.push_back(new BufferType);
  int walker_num =  W.getActiveWalkers();
  int particle_num = number_of_particles;
  int walker_dim_num = OHMMS_DIM;
  //RemoteData[0]->resize(walker_num * particle_num * walker_dim_num);
  //W.putConfigurations(RemoteData[0]->begin());
  void* walkers = RemoteData[0]->data();
  IO_VERIFY::adios_checkpoint_verify_variables(fp, "walker_num", &walker_num);
  IO_VERIFY::adios_checkpoint_verify_variables(fp, "particle_num", &particle_num);
  int walker_size = walker_num * particle_num * walker_dim_num;
  IO_VERIFY::adios_checkpoint_verify_variables(fp, "walker_size", &walker_size);
  IO_VERIFY::adios_checkpoint_verify_local_variables(fp, "walkers", (OHMMS_PRECISION *)walkers);
}
#endif
#endif


/** Write the set of walker configurations to the HDF5 file.
 * @param W set of walker configurations
 *
 * Dump is for restart and do not preserve the file
 * - version
 * - state_0
 *  - block (int)
 *  - number_of_walkes (int)
 *  - walker_partition (int array)
 *  - walkers (nw,np,3)
 */
bool HDFWalkerOutput::dump(MCWalkerConfiguration& W, int nblock)
{
  std::string FileName=myComm->getName()+hdf::config_ext;
  //rotate files
  //if(!myComm->rank() && currentConfigNumber)
  //{
  //  std::ostringstream o;
  //  o<<myComm->getName()<<".b"<<(currentConfigNumber%5) << hdf::config_ext;
  //  rename(prevFile.c_str(),o.str().c_str());
  //}

  //try to use collective
  hdf_archive dump_file(myComm,true);
  dump_file.create(FileName);
  HDFVersion cur_version;
  dump_file.write(cur_version.version,hdf::version);
  dump_file.push(hdf::main_state);
  dump_file.write(nblock,"block");

  write_configuration(W,dump_file, nblock);
  dump_file.close();

  currentConfigNumber++;
  prevFile=FileName;
  return true;
}

void HDFWalkerOutput::write_configuration(MCWalkerConfiguration& W, hdf_archive& hout, int nblock)
{
  const int wb=OHMMS_DIM*number_of_particles;
  if(nblock > block)
  {
    RemoteData[0]->resize(wb*W.getActiveWalkers());
    W.putConfigurations(RemoteData[0]->begin());
    block = nblock;
  }

  number_of_walkers=W.WalkerOffsets[myComm->size()];
  hout.write(number_of_walkers,hdf::num_walkers);

  TinyVector<int,3> gcounts(number_of_walkers,number_of_particles,OHMMS_DIM);

  if(hout.is_parallel())
  {
    { // write walker offset.
      // Though it is a small array, it needs to be written collectively in large scale runs.
      TinyVector<int,1> gcounts(myComm->size()+1);
      TinyVector<int,1> counts;
      TinyVector<int,1> offsets(myComm->rank());
      std::vector<int> myWalkerOffset;
      if(myComm->size()-1==myComm->rank())
      {
        counts[0]=2;
        myWalkerOffset.push_back(W.WalkerOffsets[myComm->rank()]);
        myWalkerOffset.push_back(W.WalkerOffsets[myComm->size()]);
      }
      else
      {
        counts[0]=1;
        myWalkerOffset.push_back(W.WalkerOffsets[myComm->rank()]);
      }
      hyperslab_proxy<std::vector<int>,1> slab(myWalkerOffset,gcounts,counts,offsets);
      hout.write(slab,"walker_partition");
    }
    { // write walker configuration
      TinyVector<int,3> counts(W.getActiveWalkers(),number_of_particles,OHMMS_DIM);
      TinyVector<int,3> offsets(W.WalkerOffsets[myComm->rank()],0,0);
      hyperslab_proxy<BufferType,3> slab(*RemoteData[0],gcounts,counts,offsets);
      hout.write(slab,hdf::walkers);
    }
  }
  else
  { //gaterv to the master and master writes it, could use isend/irecv
    hout.write(W.WalkerOffsets,"walker_partition");
    if(myComm->size()>1)
    {
      std::vector<int> displ(myComm->size()), counts(myComm->size());
      for (int i=0; i<myComm->size(); ++i)
      {
        counts[i]=wb*(W.WalkerOffsets[i+1]-W.WalkerOffsets[i]);
        displ[i]=wb*W.WalkerOffsets[i];
      }
      if(!myComm->rank())
        RemoteData[1]->resize(wb*W.WalkerOffsets[myComm->size()]);
      mpi::gatherv(*myComm,*RemoteData[0],*RemoteData[1],counts,displ);
    }
    int buffer_id=(myComm->size()>1)?1:0;
    hyperslab_proxy<BufferType,3> slab(*RemoteData[buffer_id],gcounts);
    hout.write(slab,hdf::walkers);
  }
}

/*
bool HDFWalkerOutput::dump(ForwardWalkingHistoryObject& FWO)
{
//     std::string ConfigFileName=myComm->getName()+".storeConfig.h5";
//     fw_out.open(ConfigFileName);
//
//     if (myComm->size()==1)
//     {
//       for (int i=0; i<FWO.ForwardWalkingHistory.size(); i++ )
//       {
//         int fwdata_size=FWO.ForwardWalkingHistory[i]->size();
//         std::stringstream sstr;
//         sstr<<"Block_"<<currentConfigNumber;
//         fw_out.push(sstr.str());//open the group
//
//         std::vector<float> posVecs;
//         //reserve enough space
//         std::vector<long> IDs(fwdata_size,0);
//         std::vector<long> ParentIDs(fwdata_size,0);
//         std::vector<float>::iterator tend(posVecs.begin());
//         for (int j=0;j<fwdata_size;j++)
//         {
//           const ForwardWalkingData& fwdata(FWO(i,j));
//           IDs[j]=fwdata.ID;
//           ParentIDs[j]=fwdata.ParentID;
//           fwdata.append(posVecs);
//         }
//         fw_out.write(posVecs,"Positions");
//         fw_out.write(IDs,"WalkerID");
//         fw_out.write(ParentIDs,"ParentID");
//         fw_out.write(fwdata_size,hdf::num_walkers);
//
//         fw_out.pop();//close the group
//         ++currentConfigNumber;
//       }
//     }
// #if defined(HAVE_MPI)
//     else
//     {
//       const int n3=number_of_particles*OHMMS_DIM;
//       for (int i=0; i<FWO.ForwardWalkingHistory.size(); i++ )
//       {
//         int fwdata_size=FWO.ForwardWalkingHistory[i]->size();
//         std::stringstream sstr;
//         sstr<<"Block_"<<currentConfigNumber;
//         fw_out.push(sstr.str());//open the group
//
//         std::vector<int> counts(myComm->size());
//         mpi::all_gather(*myComm,fwdata_size,counts);
//
//         std::vector<float> posVecs;
//         //reserve space to minimize the allocation
//         posVecs.reserve(FWO.number_of_walkers*n3);
//         std::vector<long> myIDs(fwdata_size),pIDs(fwdata_size);
//         std::vector<float>::iterator tend(posVecs.begin());
//         for (int j=0;j<fwdata_size;j++)
//         {
//           const ForwardWalkingData& fwdata(FWO(i,j));
//           myIDs[j]=fwdata.ID;
//           pIDs[j]=fwdata.ParentID;
//           fwdata.append(posVecs);
//         }
//         std::vector<int> offsets(myComm->size()+1,0);
//         for(int i=0; i<myComm->size();++i) offsets[i+1]=offsets[i]+counts[i];
//         fwdata_size=offsets.back();
//         fw_out.write(fwdata_size,hdf::num_walkers);
//
//         std::vector<long> globalIDs;
//         if(myComm->rank()==0) globalIDs.resize(fwdata_size);
//
//         //collect WalkerID
//         mpi::gatherv(*myComm, myIDs, globalIDs, counts, offsets);
//         fw_out.write(globalIDs,"WalkerID");
//         //collect ParentID
//         mpi::gatherv(*myComm, pIDs, globalIDs, counts, offsets);
//         fw_out.write(globalIDs,"ParentID");
//
//         for(int i=0; i<counts.size();++i) counts[i]*=n3;
//         for(int i=0; i<offsets.size();++i) offsets[i]*=n3;
//         std::vector<float> gpos;
//         if(myComm->rank()==0) gpos.resize(offsets.back());
//         mpi::gatherv(*myComm, posVecs, gpos, counts, offsets);
//
//         fw_out.write(gpos,"Positions");
//
//         fw_out.pop();//close the group
//         ++currentConfigNumber;
//       }
//     }
// #endif
//     fw_out.close();
//     FWO.clearConfigsForForwardWalking();
//
    return true;
}*/

}

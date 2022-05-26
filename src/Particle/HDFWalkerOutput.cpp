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


/** @file
 * @brief definition  of HDFWalkerOuput  and other support class and functions
 */
#include "HDFWalkerOutput.h"
#include "Utilities/IteratorUtility.h"
#include "OhmmsData/FileUtility.h"
#include "hdf/HDFVersion.h"
#include <numeric>
#include <iostream>
#include <sstream>
#include "Message/Communicate.h"
#include "mpi/collectives.h"
#include "hdf/hdf_hyperslab.h"

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
HDFWalkerOutput::HDFWalkerOutput(size_t num_ptcls, const std::string& aroot, Communicate* c)
    : appended_blocks(0),
      number_of_walkers_(0),
      number_of_particles_(num_ptcls),
      myComm(c),
      currentConfigNumber(0),
      RootName(aroot)
//       , fw_out(myComm)
{
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
  //     fw_out.write(number_of_particles_,"NumberElectrons");
  //     fw_out.write(dim,"DIM");
}

/** Destructor writes the state of random numbers and close the file */
HDFWalkerOutput::~HDFWalkerOutput()
{
  //     fw_out.close();
  delete_iter(RemoteData.begin(), RemoteData.end());
}

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
bool HDFWalkerOutput::dump(const WalkerConfigurations& W, int nblock)
{
  std::string FileName = myComm->getName() + hdf::config_ext;
  //rotate files
  //if(!myComm->rank() && currentConfigNumber)
  //{
  //  std::ostringstream o;
  //  o<<myComm->getName()<<".b"<<(currentConfigNumber%5) << hdf::config_ext;
  //  rename(prevFile.c_str(),o.str().c_str());
  //}

  //try to use collective
  hdf_archive dump_file(myComm, true);
  dump_file.create(FileName);
  HDFVersion cur_version;
  dump_file.write(cur_version.version, hdf::version);
  dump_file.push(hdf::main_state);
  dump_file.write(nblock, "block");

  write_configuration(W, dump_file, nblock);
  dump_file.close();

  currentConfigNumber++;
  prevFile = FileName;
  return true;
}

void HDFWalkerOutput::write_configuration(const WalkerConfigurations& W, hdf_archive& hout, int nblock)
{
  const int wb = OHMMS_DIM * number_of_particles_;
  if (nblock > block)
  {
    RemoteData[0]->resize(wb * W.getActiveWalkers());
    W.putConfigurations(RemoteData[0]->data());
    block = nblock;
  }
  std::cout << W.WalkerOffsets.size() << " " << myComm->size() << std::endl;
  //exit(1);
  number_of_walkers_ = W.WalkerOffsets[myComm->size()];
  hout.write(number_of_walkers_, hdf::num_walkers);

  std::array<size_t, 3> gcounts{number_of_walkers_, number_of_particles_, OHMMS_DIM};

  if (hout.is_parallel())
  {
    { // write walker offset.
      // Though it is a small array, it needs to be written collectively in large scale runs.
      std::array<size_t, 1> gcounts{static_cast<size_t>(myComm->size()) + 1};
      std::array<size_t, 1> counts{0};
      std::array<size_t, 1> offsets{static_cast<size_t>(myComm->rank())};
      std::vector<int> myWalkerOffset;
      if (myComm->size() - 1 == myComm->rank())
      {
        counts[0] = 2;
        myWalkerOffset.push_back(W.WalkerOffsets[myComm->rank()]);
        myWalkerOffset.push_back(W.WalkerOffsets[myComm->size()]);
      }
      else
      {
        counts[0] = 1;
        myWalkerOffset.push_back(W.WalkerOffsets[myComm->rank()]);
      }
      hyperslab_proxy<std::vector<int>, 1> slab(myWalkerOffset, gcounts, counts, offsets);
      hout.write(slab, "walker_partition");
    }
    { // write walker configuration
      std::array<size_t, 3> counts{W.getActiveWalkers(), number_of_particles_, OHMMS_DIM};
      std::array<size_t, 3> offsets{static_cast<size_t>(W.WalkerOffsets[myComm->rank()]), 0, 0};
      hyperslab_proxy<BufferType, 3> slab(*RemoteData[0], gcounts, counts, offsets);
      hout.write(slab, hdf::walkers);
    }
  }
  else
  { //gaterv to the master and master writes it, could use isend/irecv
    hout.write(W.WalkerOffsets, "walker_partition");
    if (myComm->size() > 1)
    {
      std::vector<int> displ(myComm->size()), counts(myComm->size());
      for (int i = 0; i < myComm->size(); ++i)
      {
        counts[i] = wb * (W.WalkerOffsets[i + 1] - W.WalkerOffsets[i]);
        displ[i]  = wb * W.WalkerOffsets[i];
      }
      if (!myComm->rank())
        RemoteData[1]->resize(wb * W.WalkerOffsets[myComm->size()]);
      mpi::gatherv(*myComm, *RemoteData[0], *RemoteData[1], counts, displ);
    }
    int buffer_id = (myComm->size() > 1) ? 1 : 0;
    hout.writeSlabReshaped(*RemoteData[buffer_id], gcounts, hdf::walkers);
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
//       const int n3=number_of_particles_*OHMMS_DIM;
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
//         posVecs.reserve(FWO.number_of_walkers_*n3);
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

} // namespace qmcplusplus

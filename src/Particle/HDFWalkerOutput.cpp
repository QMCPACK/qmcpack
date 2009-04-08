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
#include "Utilities/IteratorUtility.h"
#include "OhmmsData/FileUtility.h"
#include "Message/Communicate.h"
#include "Message/CommOperators.h"
#include "Message/CommUtilities.h"
#include "Particle/HDFWalkerIOEngine.h"
#include "Numerics/HDFSTLAttrib.h"

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
    hid_t h1;
    bool overwrite=true;

    if (block)
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
        appended_blocks(0), number_of_walkers(0), currentConfigNumber(0),
        number_of_backups(0), max_number_of_backups(4),
        h_file(-1), h_plist(H5P_DEFAULT), xfer_plist(H5P_DEFAULT), h_state(-1), myComm(c), c_file(-1), c_state(-1)
{
    FileName=myComm->getName()+hdf::config_ext;
    ConfigFileName=myComm->getName()+".storeConfig.h5";

    if (myComm->rank()==0)
    {
        c_file = H5Fcreate(ConfigFileName.c_str(),H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);

        HDFVersion cur_version;
        cur_version.write(c_file,hdf::version);
        int nel = W[0]->R.size();
        HDFAttribIO<int> nwo(nel);
        nwo.write(c_file,"NumberElectrons");
        int dim=OHMMS_DIM;
        HDFAttribIO<int> ndo(dim);
        ndo.write(c_file,"DIM");
        H5Fclose(c_file);
    }

    if (myComm->size()==1)
    {
        std::pair<hid_t,hid_t> h_pair= createH5FileSingle(FileName,W);
        h_file=h_pair.first;
        h_state=h_pair.second;
        number_of_walkers = W.getActiveWalkers();
        HDFAttribIO<int> nwo(number_of_walkers);
        nwo.write(h_state,hdf::num_walkers);
        H5Gclose(h_state);
        H5Fclose(h_file);

//      c_file = H5Fcreate(ConfigFileName.c_str(),H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);

//      HDFVersion cur_version;
//      cur_version.write(c_file,hdf::version);
//      c_state = H5Gcreate(c_file,hdf::main_state,0);
//      H5Gclose(c_state);
//      H5Fclose(c_file);

    }
    else
    {
#if defined(H5_HAVE_PARALLEL) && defined(ENABLE_PHDF5)
        //open file collectively
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

        HDFAttribIO<int> nwo(number_of_walkers);
        nwo.write(h_state,hdf::num_walkers);
        H5Gclose(h_state);
        H5Fclose(h_file);
#else
        //add two buffers: local and collected
        if (RemoteData.empty())
        {//add two buffers
            RemoteData.push_back(new BufferType);
            RemoteData.push_back(new BufferType);
            FWData.push_back(new  FWBufferType);
            FWData.push_back(new  FWBufferType);
            FWCountData.resize(myComm->size());
        }
        //myRequest.resize(myComm->size()-1);
        //if(myComm->rank())
        //{
        //  RemoteData.push_back(new BufferType);
        //}
        //else
        //{
        //  RemoteData.resize(myRequest.size(),0);
        //  for(int p=0; p < myRequest.size(); ++p) RemoteData[p] = new BufferType;
        //}
        dumpCollect(W);

        ////for debugging only. Open a file per node
        //string fname=aroot+".config.h5";
        //std::pair<hid_t,hid_t> h_pair= createH5FileSingle(fname,W);
        //h_debug_file=h_pair.first;
        //H5Gclose(h_pair.second);
#endif
    }

}

/** Destructor writes the state of random numbers and close the file */
HDFWalkerOutput::~HDFWalkerOutput()
{
    if (xfer_plist!= H5P_DEFAULT) H5Pclose(xfer_plist);
    if (h_plist!= H5P_DEFAULT) H5Pclose(h_plist);
    delete_iter(RemoteData.begin(),RemoteData.end());
}


/** Write the set of walker configurations to the HDF5 file.
 * @param W set of walker configurations
 */
bool HDFWalkerOutput::dump(MCWalkerConfiguration& W) {
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
    bool success=true;
    const bool overwrite=true;
    if (myComm->size()==1)
    {//dump it
        std::pair<hid_t,hid_t> h_pair= createH5FileSingle(FileName,W);
        h_file=h_pair.first;
        h_state=h_pair.second;
        number_of_walkers = W.getActiveWalkers();
        HDFAttribIO<int> nwo(number_of_walkers);
        nwo.write(h_state,hdf::num_walkers);
        H5Gclose(h_state);
        H5Fclose(h_file);
    }
    else
    {
#if defined(H5_HAVE_PARALLEL) && defined(ENABLE_PHDF5)
        if (number_of_walkers != W.getGlobalNumWalkers())
        {
            h_file = H5Fcreate(FileName.c_str(),H5F_ACC_TRUNC,H5P_DEFAULT,h_plist);
            HDFVersion cur_version;
            cur_version.write(h_file,hdf::version);

            //create a state
            h_state = H5Gcreate(h_file,hdf::main_state,0);
            HDFWalkerIOEngine wo(W);
            wo.writeAll(h_state,hdf::walkers,myComm);
            number_of_walkers = W.getGlobalNumWalkers();

            HDFAttribIO<int> nwo(number_of_walkers);
            nwo.write(h_state,hdf::num_walkers);
        }
        else
        {
            h_file =  H5Fopen(FileName.c_str(),H5F_ACC_RDWR,h_plist);
            h_state =  H5Gopen(h_file,hdf::main_state);
            //if(number_of_walkers != W.getGlobalNumWalkers())
            //{
            //  herr_t status=H5Gunlink(h_state,hdf::walkers);
            //  HDFWalkerIOEngine wo(W);
            //  wo.setTransferProperty(xfer_plist);
            //  wo.writeAll(h_state,hdf::walkers,myComm);

            //  // overwrite number of walkers
            //  number_of_walkers=W.getGlobalNumWalkers();
            //  HDFAttribIO<int> nwo(number_of_walkers,true);
            //  nwo.write(h_state,hdf::num_walkers);
            //}
            //else
            //{
            HDFWalkerIOEngine wo(W,overwrite);
            wo.setTransferProperty(xfer_plist);
            wo.writeAll(h_state,hdf::walkers,myComm);
            //}
        }
        H5Gclose(h_state);
        H5Fclose(h_file);
#else
        success=dumpCollect(W);
#endif
        //char fname[128];
        //sprintf(fname,"%s.p%03d%s",myComm->getName().c_str(), myComm->rank(),hdf::config_ext);
        //std::pair<hid_t,hid_t> h_pair= createH5FileSingle(fname,W);
        //int nw = W.getActiveWalkers();
        //HDFAttribIO<int> nwo(nw);
        //nwo.write(h_pair.first,hdf::num_walkers);
        //H5Gclose(h_pair.second);
        //H5Fclose(h_pair.first);
    }

    //app_log() << "  Dumped " << number_of_walkers << " walkers for restart." << endl;
    //appended_blocks=appendSingle(h_debug_file,W,appended_blocks);
    return success;
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
    if (number_of_walkers != nw)
    {//need resize of the data: unlink and create a walkers group
        number_of_walkers=nw;
        HDFAttribIO<int> nwo(number_of_walkers,true);
        nwo.write(h_state,hdf::num_walkers);

        herr_t status=H5Gunlink(gid,hdf::walkers);
        HDFWalkerIOEngine wo(W);
        wo.write(gid,hdf::walkers);
    }
    else
    {//simply overwrite it
        const bool replace=true;
        HDFWalkerIOEngine wo(W,replace);
        wo.write(h_state,hdf::walkers);
    }
    return true;
}

/** Write the set of walker configurations to the HDF5 file.
 * @param W set of walker configurations
 * @param block index for the block
 */
bool  HDFWalkerOutput::append(MCWalkerConfiguration& W) {
//    hid_t h0=h_file;
//    if(h0<0) h0 =  H5Fopen(FileName.c_str(),H5F_ACC_RDWR,H5P_DEFAULT);
    h_file = H5Fopen(FileName.c_str(),H5F_ACC_RDWR,H5P_DEFAULT);
    appended_blocks=appendSingle(h_file,W,appended_blocks);
//    if(h0 != h_file) H5Fclose(h0);
    if (H5Fclose(h_file) > -1) h_file = -1;
    return true;
}

bool  HDFWalkerOutput::dumpCollect(MCWalkerConfiguration& W)
{
#if defined(HAVE_MPI)
    int wb=OHMMS_DIM*W.getTotalNum();
    RemoteData[0]->resize(wb*W.getActiveWalkers());
    W.putConfigurations(RemoteData[0]->begin());
    vector<int> displ(myComm->size()), counts(myComm->size());
    for (int i=0; i<myComm->size(); ++i)
    {
        counts[i]=wb*(W.WalkerOffsets[i+1]-W.WalkerOffsets[i]);
        displ[i]=wb*W.WalkerOffsets[i];
    }
    RemoteData[1]->resize(wb*W.WalkerOffsets[myComm->size()]);
    myComm->gatherv(*RemoteData[0],*RemoteData[1],counts, displ);
    ////could use gatherv but use send/irecv for now
    ////int nw_max=W.WalkerOffSets[myComm->size()]-W.WalkerOffsets[myComm->size()-1];
    //if(myComm->rank())
    //{
    //  //RemoteData[0]->resize(wb*nw_max);
    //  RemoteData[0]->resize(wb*W.getActiveWalkers());
    //  W.putConfigurations(RemoteData[0]->begin());
    //  myComm->send(0,2000,*RemoteData[0]);
    //}
    //else
    //{
    //  for(int p=0; p<myRequest.size(); ++p)
    //  {
    //    RemoteData[p]->resize(wb*(W.WalkerOffsets[p+2]-W.WalkerOffsets[p+1]));
    //    myRequest[p]=myComm->irecv(p+1,2000,*RemoteData[p]);
    //    //RemoteData[p]->resize(wb*nw_max);
    //    //myRequest[p]=myComm->irecv(MPI_ANY_SOURCE,2000,*RemoteData[p]);
    //  }

    //  number_of_walkers=W.WalkerOffsets[myComm->size()];
    //  BufferType b(wb*number_of_walkers);
    //  W.putConfigurations(b.begin());

    //  vector<Communicate::status> myStatus(myRequest.size());
    //  MPI_Waitall(myRequest.size(),&myRequest[0],&myStatus[0]);
    //  for(int p=0; p<myRequest.size(); ++p)
    //  {
    //    std::copy(RemoteData[p]->begin(),RemoteData[p]->end(),b.begin()+wb*W.WalkerOffsets[p+1]);
    //    //int s=myStatus[p].MPI_SOURCE;
    //    //int nw_recv=W.WalkerOffsets[s+1]-W.WalkerOffsets[s];
    //    //std::copy(RemoteData[p]->begin(),RemoteData[p]->begin()+nw_recv*wb,b.begin()+wb*W.WalkerOffsets[s]);
    //  }
    number_of_walkers=W.WalkerOffsets[myComm->size()];
    if (myComm->rank()==0)
    {
        hid_t d1= H5Fcreate(FileName.c_str(),H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
        HDFVersion cur_version;
        cur_version.write(d1,hdf::version);
        hid_t d2 = H5Gcreate(d1,hdf::main_state,0);
        vector<int> inds(3);
        inds[0]=number_of_walkers;
        inds[1]=W.getTotalNum();
        inds[2]=OHMMS_DIM;
        //HDFAttribIO<BufferType> po(b,inds);
        HDFAttribIO<BufferType> po(*RemoteData[1],inds);
        po.write(d2,hdf::walkers);
        HDFAttribIO<int> nwo(number_of_walkers);
        nwo.write(d2,hdf::num_walkers);
        H5Gclose(d2);
        if (H5Fclose(d1) > -1) d1 = -1;
    }
#endif

    return true;
}

void HDFWalkerOutput::open()
{
    h_file =  H5Fopen(FileName.c_str(),H5F_ACC_RDWR,h_plist);
    h_state = H5Gopen(h_file,hdf::main_state);
}

bool HDFWalkerOutput::dump(ForwardWalkingHistoryObject& FWO)
{
    if (myComm->size()>1)
    {
#if defined(HAVE_MPI)
        int NN(myComm->size());
        int RN(myComm->rank());
        vector<int> myWalkerLayout;
        FWO.layoutOfConfigsForForwardWalking(myWalkerLayout);

        int Nblocks = myWalkerLayout.size();
        vector<int> walkerCounts(NN,Nblocks);
        vector<int> walkerOffsets(NN,0);

        for (int i=0;i<NN;i++) walkerOffsets[i]=Nblocks*i;
        vector<int> TFWCountData(Nblocks*NN,0);
        myComm->allgatherv( myWalkerLayout, TFWCountData , walkerCounts , walkerOffsets);


        vector<vector<int> > OFFSETS(NN,vector<int>(Nblocks,0)), COUNTS(NN,vector<int>(Nblocks,0)), POSCOUNTS(NN,vector<int>(Nblocks,0)), POSOFFSETS(NN,vector<int>(Nblocks,0));
        int Nstored = FWO.ForwardWalkingHistory.size();
        int nelecs = FWO.ForwardWalkingHistory[0][0].Pos.size();
        int nFloatsPerConfig = OHMMS_DIM*nelecs;


//     if(myComm->rank()==0)
//     {
        for (int i=0;i<NN;i++) for (int j=0;j<Nblocks;j++) COUNTS[j][i]=TFWCountData[i*Nblocks+j];
        for (int i=0;i<Nblocks;i++) for (int j=0;j<(NN-1);j++) OFFSETS[i][j+1] = COUNTS[i][j] + OFFSETS[i][j];
        for (int i=0;i<NN;i++) for (int j=0;j<Nblocks;j++) POSOFFSETS[j][i] = nFloatsPerConfig*OFFSETS[j][i];
        for (int i=0;i<NN;i++) for (int j=0;j<Nblocks;j++) POSCOUNTS[j][i] = nFloatsPerConfig*COUNTS[j][i];
//     }
//     for (int i=0;i<Nblocks;i++)
//     {
//       myComm->bcast(OFFSETS[i]);  /*cerr<<"OFFSET"<<endl;*/
//       myComm->bcast(POSOFFSETS[i]); /*cerr<<"POSOFFSETS"<<endl;*/
//       myComm->bcast(COUNTS[i]);  /*cerr<<"COUNTS"<<endl;*/
//       myComm->bcast(POSCOUNTS[i]); /*cerr<<"POSCOUNTS"<<endl;*/
//     }
//     if(myComm->rank()==0) for(int i=0;i<TFWCountData.size();i++) cout<<TFWCountData[i]<<" ";
//     cout<<endl;
//     if(myComm->rank()==0) for(int i=0;i<myWalkerLayout.size();i++) cout<<myWalkerLayout[i]<<" ";
//     cout<<endl;
//     if(myComm->rank()==0)
//     {
//       for(int i=0;i<NN;i++) for (int j=1;j<Nblocks;j++) cerr<<COUNTS[j][i]<<" ";
//       cerr<<endl;
//       for(int i=0;i<NN;i++) for (int j=1;j<Nblocks;j++) cerr<<OFFSETS[j][i]<<" ";
//       cerr<<endl;
//       for(int i=0;i<NN;i++) for (int j=1;j<Nblocks;j++) cerr<<POSCOUNTS[j][i]<<" ";
//       cerr<<endl;
//       for(int i=0;i<NN;i++) for (int j=1;j<Nblocks;j++) cerr<<POSOFFSETS[j][i]<<" ";
//       cerr<<endl;
//     }
//     cerr<<endl;


//     cerr<<Nstored<<" "<<nelecs<<endl;
        for (int i=0;i<Nstored;i++)
        {
            int globalNPositions = POSCOUNTS[i][POSCOUNTS[i].size()-1] + POSOFFSETS[i][POSOFFSETS[i].size()-1];
            int globalNIDs = COUNTS[i][COUNTS[i].size()-1] + OFFSETS[i][OFFSETS[i].size()-1];
            vector<float> myPositions(POSCOUNTS[i][RN],0),globalPositions(globalNPositions,0);
            vector<long> myIDs(COUNTS[i][RN],0),globalIDs(globalNIDs,0);
            vector<long> myParentIDs(COUNTS[i][RN],0),globalParentIDs(globalNIDs,0);
            myComm->gatherv( myIDs, globalIDs , COUNTS[i] , OFFSETS[i]);
            myComm->gatherv( myParentIDs, globalParentIDs , COUNTS[i] , OFFSETS[i]);
            myComm->gatherv( myPositions, globalPositions , POSCOUNTS[i] , POSOFFSETS[i]);
            //Now I have gathered all the positions and IDs to the head node. I just need to write them.
            if (myComm->rank()==0)
            {
                ++currentConfigNumber;
                c_file = H5Fopen(ConfigFileName.c_str(),H5F_ACC_RDWR,H5P_DEFAULT);
                std::stringstream sstr("");
                sstr<<"Block_"<<currentConfigNumber;
                hid_t d_file = H5Gcreate(c_file,sstr.str().c_str(),0);
                sstr.str("WalkerID");
                string groupName = sstr.str();
                const int IDrank = 1;
                hsize_t IDdims[IDrank], IDmaxdims[IDrank];
                IDdims[0] = globalNIDs;
                IDmaxdims[0] = H5S_UNLIMITED;

                hid_t dataspace  = H5Screate_simple(IDrank, IDdims, IDmaxdims);
                hid_t p = H5Pcreate (H5P_DATASET_CREATE);
                H5Pset_chunk(p,IDrank,IDdims);
                hid_t dataset =  H5Dcreate(d_file, groupName.c_str(), H5T_NATIVE_LONG, dataspace, p);
                hid_t memspace = H5Screate_simple(IDrank, IDdims, NULL);
                herr_t ret = H5Dwrite(dataset, H5T_NATIVE_LONG, memspace, dataspace, H5P_DEFAULT,&(globalIDs[0]));
//           H5Sclose(memspace);
//           H5Sclose(dataspace);
//           H5Dclose(dataset);

                sstr.str("Positions");
                groupName = sstr.str();

                const int rank = 1;
                hsize_t dims[rank], maxdims[rank];
                dims[0] = globalNPositions; /*dims[1] = nelecs; dims[2] = OHMMS_DIM;*/
                maxdims[0] = H5S_UNLIMITED; /*maxdims[1] = globalNIDs; maxdims[2] = nelecs;*/

                dataspace  = H5Screate_simple(rank, dims, maxdims);
                p = H5Pcreate (H5P_DATASET_CREATE);
                H5Pset_chunk(p,rank,dims);
                dataset =  H5Dcreate(d_file, groupName.c_str(), H5T_NATIVE_FLOAT, dataspace, p);
                memspace = H5Screate_simple(rank, dims, NULL);
                ret = H5Dwrite(dataset, H5T_NATIVE_FLOAT, memspace, dataspace, H5P_DEFAULT,&(globalPositions[0]));



                sstr.str("ParentID");
                groupName = sstr.str();

                dataspace  = H5Screate_simple(IDrank, IDdims, IDmaxdims);
                p = H5Pcreate (H5P_DATASET_CREATE);
                H5Pset_chunk(p,IDrank,IDdims);
                dataset =  H5Dcreate(d_file, groupName.c_str(), H5T_NATIVE_LONG, dataspace, p);
                memspace = H5Screate_simple(IDrank, IDdims, NULL);
                ret = H5Dwrite(dataset, H5T_NATIVE_LONG, memspace, dataspace, H5P_DEFAULT,&(globalParentIDs[0]));

                HDFAttribIO<int> nwo(globalNIDs);
                nwo.write(d_file,hdf::num_walkers);

                H5Sclose(memspace);
                H5Sclose(dataspace);
                H5Dclose(dataset);
                H5Fclose(d_file);
                H5Fclose(c_file);
            }

        }

#endif
    }
    else if (myComm->size()==1)
    {
        int Nblocks = FWO.ForwardWalkingHistory.size(), totWalkers(0);
        ///Nblocks is the number of records we collected since the last configuration dump.
        vector<int> NWalkersInBlock(Nblocks,0);
        ///NWalkersInBlock is the number of Walkers we saved for a single record
        for (int i=0;i<Nblocks;i++) {
            totWalkers+=FWO.ForwardWalkingHistory[i].size();
            NWalkersInBlock[i] = FWO.ForwardWalkingHistory[i].size();
        }
        ///nelecs is the number of R vectors we are going to have to write
        int nelecs = FWO.ForwardWalkingHistory[0][0].Pos.size();
        c_file = H5Fopen(ConfigFileName.c_str(),H5F_ACC_RDWR,H5P_DEFAULT);

        int Gsize= (2*sizeof(long)+sizeof(float)*nelecs*OHMMS_DIM)*totWalkers;

        typedef ForwardWalkingData::StoredPosType StoredPosType;
        for (int i=0; i<Nblocks; i++, ++currentConfigNumber)
        {
            std::stringstream sstr("");
            sstr<<"Block_"<<currentConfigNumber;
            hid_t d_file = H5Gcreate(c_file,sstr.str().c_str(),Gsize);
            sstr.str("Positions");
            string groupName = sstr.str();
//           Matrix<StoredPosType> tp(NWalkersInBlock[i],nelecs);
//           for (int j=0;j<NWalkersInBlock[i];j++)
//             for (int k=0;k<nelecs;k++)
//               tp[j][k] = FWO.ForwardWalkingHistory[i][j].Pos[k];

            vector<float> posVecs;
            for (int j=0;j<NWalkersInBlock[i];j++)
            {
                vector<float> tposVec;
                (FWO.ForwardWalkingHistory[i][j]).toFloat(tposVec);
                posVecs.insert(posVecs.end(),tposVec.begin(),tposVec.end());
            }

//           const int rank = 3;
//           hsize_t dims[rank], maxdims[rank];
//           dims[0] = NWalkersInBlock[i]; dims[1] = nelecs; dims[2] = StoredPosType::Size;
//           maxdims[0] = H5S_UNLIMITED; maxdims[1] = tp.cols(); maxdims[2] = StoredPosType::Size;

            const int rank = 1;
            hsize_t dims[rank], maxdims[rank];
            dims[0] = posVecs.size(); /*dims[1] = nelecs; dims[2] = StoredPosType::Size;*/
            maxdims[0] = H5S_UNLIMITED; /*maxdims[1] = tp.cols(); maxdims[2] = StoredPosType::Size;*/

            hid_t dataspace  = H5Screate_simple(rank, dims, maxdims);
            hid_t p = H5Pcreate (H5P_DATASET_CREATE);
            H5Pset_chunk(p,rank,dims);
            hid_t dataset =  H5Dcreate(d_file, groupName.c_str(), H5T_NATIVE_FLOAT, dataspace, p);
            hid_t memspace = H5Screate_simple(rank, dims, NULL);
            herr_t ret = H5Dwrite(dataset, H5T_NATIVE_FLOAT, memspace, dataspace, H5P_DEFAULT,&(posVecs[0]));
//           H5Sclose(memspace);
//           H5Sclose(dataspace);
//           H5Dclose(dataset);

            sstr.str("WalkerID");
            groupName = sstr.str();
            vector<long> IDs(NWalkersInBlock[i]);
            for (int j=0;j<NWalkersInBlock[i];j++)  IDs[j]=FWO.ForwardWalkingHistory[i][j].ID;

            const int IDrank = 1;
            hsize_t IDdims[IDrank], IDmaxdims[IDrank];
            IDdims[0] = NWalkersInBlock[i];
            IDmaxdims[0] = H5S_UNLIMITED;

            dataspace  = H5Screate_simple(IDrank, IDdims, IDmaxdims);
            p = H5Pcreate (H5P_DATASET_CREATE);
            H5Pset_chunk(p,IDrank,IDdims);
            dataset =  H5Dcreate(d_file, groupName.c_str(), H5T_NATIVE_LONG, dataspace, p);
            memspace = H5Screate_simple(IDrank, IDdims, NULL);
            ret = H5Dwrite(dataset, H5T_NATIVE_LONG, memspace, dataspace, H5P_DEFAULT,&(IDs[0]));
//           H5Sclose(memspace);
//           H5Sclose(dataspace);
//           H5Dclose(dataset);

            sstr.str("ParentID");
            groupName = sstr.str();
            vector<long> ParentIDs(NWalkersInBlock[i]);
            for (int j=0;j<NWalkersInBlock[i];j++)  ParentIDs[j]=FWO.ForwardWalkingHistory[i][j].ParentID;

            dataspace  = H5Screate_simple(IDrank, IDdims, IDmaxdims);
            p = H5Pcreate (H5P_DATASET_CREATE);
            H5Pset_chunk(p,IDrank,IDdims);
            dataset =  H5Dcreate(d_file, groupName.c_str(), H5T_NATIVE_LONG, dataspace, p);
            memspace = H5Screate_simple(IDrank, IDdims, NULL);
            ret = H5Dwrite(dataset, H5T_NATIVE_LONG, memspace, dataspace, H5P_DEFAULT,&(ParentIDs[0]));
            H5Sclose(memspace);
            H5Sclose(dataspace);
            H5Dclose(dataset);
            H5Fclose(d_file);
        }
        H5Fclose(c_file);
    }
    FWO.clearConfigsForForwardWalking();
    return true;
}

}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/

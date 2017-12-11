//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



#include "Particle/MCWalkerConfiguration.h"
#include "Particle/HDFWalkerInputCollect.h"
#include "OhmmsPETE/OhmmsVector.h"
#if defined(HAVE_LIBHDF5)
#include "Particle/HDFParticleAttrib.h"
#include "Numerics/HDFNumericAttrib.h"
#endif
#include "Utilities/RandomGenerator.h"
#include "OhmmsData/FileUtility.h"
#include "Utilities/RandomGeneratorIO.h"

using namespace qmcplusplus;

/** Open the HDF5 file "aroot.config.h5" for reading.
 *@param aroot the root file name
 *
 *@note The main group is "/config_collection"
 */
HDFWalkerInputCollect::HDFWalkerInputCollect(const std::string& aroot):
  fileID(-1), prevNContexts(1), curNContexts(1), CollectMode(true), RewindMode(false)
{
  open(aroot);
}

HDFWalkerInputCollect::~HDFWalkerInputCollect()
{
  close();
}

#if defined(HAVE_LIBHDF5)
bool
HDFWalkerInputCollect::open(const std::string& aroot)
{
  std::string h5file = aroot;
  std::string ext=getExtension(h5file);
  if(ext != "h5")
    //if the filename does not h5 extension, add the extension
  {
    h5file.append(".config.h5");
  }
  fileID =  H5Fopen(h5file.c_str(),H5F_ACC_RDWR,H5P_DEFAULT);
  if(fileID<0)
  {
    app_error() << "Cannot open hdf5 file " << std::endl;
    return false;
  }
  app_log() << "  HDFWalkerInputCollect::open " << h5file<< std::endl;
  curNContexts = OHMMS::Controller->size();
  OffSet.resize(curNContexts+1);
  return true;
}

bool HDFWalkerInputCollect::close()
{
  if(fileID >-1)
  {
    if(!RewindMode)
      readRandomState();
    app_log() << "  HDFWalkerInputCollect::close " << std::endl;
    //read the random seeds
    H5Fclose(fileID);
    fileID=-1;
  }
  return true;
}

void HDFWalkerInputCollect::readRandomState()
{
  prevNContexts=1;
  herr_t status = H5Eset_auto(NULL, NULL);
  hid_t h1=H5Dopen(fileID,"ncontexts");
  if(h1>-1)
  {
    H5Dread(h1, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(prevNContexts));
    H5Dclose(h1);
    char rname[128];
    if(prevNContexts == curNContexts)
    {
      app_log() << "    Restart with the random states" << std::endl;
      sprintf(rname,"context%04d/random_state",OHMMS::Controller->rank());
      hid_t h_random = H5Gopen(fileID,rname);
      HDFAttribIO<RandomGenerator_t> r(Random);
      r.read(h_random,"dummy");
      //Random.read(h_random);
      //H5Gclose(h_random);
    }
    else
    {
      app_warning() << "The number of processors has changed.\n"
                    << "New random seeds are generated." << std::endl;
    }
  }
  else
  {
    if(curNContexts == 1)
    {
      hid_t h_random = H5Gopen(fileID,"random_state");
      if(h_random>-1)
      {
        app_log() << "Using serial random seed" << std::endl;
        HDFAttribIO<RandomGenerator_t> r(Random);
        r.read(h_random,"dummy");
        //Random.read(h_random);
        //H5Gclose(h_random);
      }
    }
  }
}

bool
HDFWalkerInputCollect::put(MCWalkerConfiguration& W, int rollback)
{
  int numConfigs=1;
  hid_t mastercf = H5Gopen(fileID,"config_collection");
  hid_t h1=H5Dopen(mastercf,"NumOfConfigurations");
  if(h1>-1)
  {
    H5Dread(h1, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&(numConfigs));
    H5Dclose(h1);
  }
  app_log() << "    The number of blocks in the file " << numConfigs << std::endl;
  int firstConf(0);
  if(rollback>0)
  {
    if(rollback<numConfigs)
      firstConf=numConfigs-rollback;
    else
      firstConf = 0;
  }
  H5Gclose(mastercf);
  RewindMode=false;
  return read(W,firstConf,numConfigs);
}

bool HDFWalkerInputCollect::putSingle(MCWalkerConfiguration& W)
{
  CollectMode=false;
  return put(W,-1);
}

bool
HDFWalkerInputCollect::rewind(MCWalkerConfiguration& W, int rollback)
{
  int numConfigs=1;
  hid_t mastercf = H5Gopen(fileID,"config_collection");
  herr_t status = H5Eset_auto(NULL, NULL);
  hid_t h1=H5Dopen(mastercf,"NumOfConfigurations");
  if(h1>-1)
  {
    H5Dread(h1, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&(numConfigs));
    H5Dclose(h1);
  }
  H5Gclose(mastercf);
  app_log() << "    The number of blocks in the file " << numConfigs << std::endl;
  if(numConfigs ==1)
  {
    app_warning() << "HDFWalkerInputCollect::rewind fails: not enough configurations" << std::endl;
    app_warning() << "The number of walkers is unchanged." << std::endl;
    return false;
  }
  RewindMode=true;
  int firstConf=std::max(0,numConfigs-rollback);
  //check if the file has ncontexts
  prevNContexts=1;
  h1=H5Dopen(fileID,"ncontexts");
  if(h1>-1)
  {
    H5Dread(h1, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(prevNContexts));
    H5Dclose(h1);
  }
  if(prevNContexts>1)
  {
    app_log() << "  rewind of a shared file" << std::endl;
    CollectMode=true;
  }
  else
  {
    app_log() << "  rewind of a independent file" << std::endl;
    CollectMode=false;
  }
  return read(W,firstConf,numConfigs-1);
}

bool
HDFWalkerInputCollect::read(MCWalkerConfiguration& W, int firstConf, int lastConf)
{
  int myID = OHMMS::Controller->rank();
  hid_t mastercf = H5Gopen(fileID,"config_collection");
  char confName[128];
  char coordName[128];
#if H5_VERS_RELEASE < 4
  hssize_t offset[3];
#else
  hsize_t offset[3];
#endif
  hsize_t dimIn[3],dimTot[3];
  offset[0]=0;
  offset[1]=0;
  offset[2]=0;
  typedef MCWalkerConfiguration::PosType PosType;
  std::vector<PosType> pos;
  int nwRead=0;
  for(int iconf=firstConf; iconf<lastConf; iconf++)
  {
    sprintf(coordName,"config%04d/coord",iconf);
    hid_t dataset = H5Dopen(mastercf,coordName);
    hid_t dataspace = H5Dget_space(dataset);
    int rank = H5Sget_simple_extent_ndims(dataspace);
    int status_n = H5Sget_simple_extent_dims(dataspace, dimTot, NULL);
    if(CollectMode)
    {
      distribute(dimTot[0]);
    }
    else
    {
      OffSet[myID]=0;
      OffSet[myID+1]=dimTot[0];
    }
    //get the input dimension
    dimIn[0]=OffSet[myID+1]-OffSet[myID];
    dimIn[1]=dimTot[1];
    dimIn[2]=dimTot[2];
    offset[0]=OffSet[myID];
    std::vector<PosType> posIn(dimIn[0]*dimIn[1]);
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
  if(curWalker)
  {
    W.createWalkers(nwRead);
  }
  else
  {
    W.resize(nwRead,nptcl);
  }
  MCWalkerConfiguration::iterator it = W.begin()+curWalker;
  int ii=0;
  for(int iw=0; iw<nwRead; iw++)
  {
    //std::copy(Post_temp[iw],Post_temp[iw+1], (*it)->R.begin());
    for(int iat=0; iat < nptcl; iat++,ii++)
    {
      (*it)->R(iat) = pos[ii];
    }
    ++it;
  }
  return true;
}


void HDFWalkerInputCollect::distribute(hsize_t nw)
{
  hsize_t bat=nw/curNContexts;
  hsize_t residue = curNContexts-nw%curNContexts;
  OffSet[0] = 0;
  for(int i=0; i<curNContexts; i++)
  {
    if(i<residue)
      OffSet[i+1] = OffSet[i] + bat;
    else
      OffSet[i+1] = OffSet[i] + bat+1;
  }
}
#else
bool
HDFWalkerInputCollect::open(const std::string& aroot)
{
  return false;
}

bool HDFWalkerInputCollect::close()
{
  return false;
}

void HDFWalkerInputCollect::readRandomState()
{
}

bool
HDFWalkerInputCollect::put(MCWalkerConfiguration& W, int rollback)
{
  return false;
}

bool HDFWalkerInputCollect::putSingle(MCWalkerConfiguration& W)
{
  return false;
}

bool
HDFWalkerInputCollect::rewind(MCWalkerConfiguration& W, int rollback)
{
  return false;
}

bool
HDFWalkerInputCollect::read(MCWalkerConfiguration& W, int firstConf, int lastConf)
{
  return false;
}

void HDFWalkerInputCollect::distribute(hsize_t nw)
{
}
#endif

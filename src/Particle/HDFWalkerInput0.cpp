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
#include "Particle/HDFWalkerInput0.h"
#include "OhmmsPETE/OhmmsVector.h"
#if defined(HAVE_LIBHDF5)
#include "Particle/HDFParticleAttrib.h"
#include "Numerics/HDFNumericAttrib.h"
#endif
#include "Utilities/RandomGenerator.h"
#include "Utilities/RandomGeneratorIO.h"
#include "OhmmsData/FileUtility.h"
using namespace qmcplusplus;


#if defined(HAVE_LIBHDF5)
/** Open the HDF5 file "aroot.config.h5" for reading.
 *@param aroot the root file name
 *
 *@note The main group is "/config_collection"
 */
HDFWalkerInput0::HDFWalkerInput0(const std::string& aroot, int ipart, int nparts):
  Counter(0), NumSets(0)
{
  std::string h5file = aroot;
  std::string ext=getExtension(h5file);
  if(ext != "h5")
    //if the filename does not h5 extension, add the extension
  {
    h5file.append(".config.h5");
  }
  h_file =  H5Fopen(h5file.c_str(),H5F_ACC_RDWR,H5P_DEFAULT);
  h_config = H5Gopen(h_file,"config_collection");
  hid_t h1=H5Dopen(h_config,"NumOfConfigurations");
  if(h1>-1)
  {
    H5Dread(h1, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&(NumSets));
    H5Dclose(h1);
  }
  if(!NumSets)
  {
    //resolve the integer and long problem with 64bit
    hsize_t nset;
    H5Gget_num_objs(h_config,&nset);
    NumSets=nset;
  }
  if(NumSets)
  {
    int nframe=NumSets/nparts;
    FirstSet=nframe*ipart;
    LastSet=std::max(nframe*(ipart+1),NumSets);
  }
  else
  {
    FirstSet=0;
    LastSet=0;
    ERRORMSG("File does not contain walkers!")
  }
}

/** Destructor closes the HDF5 file and main group. */
HDFWalkerInput0::~HDFWalkerInput0()
{
  H5Gclose(h_config);
  H5Fclose(h_file);
}

/** Read a configuration from a hdf5 file.
 *@param W set of walker configurations
int HDFWalkerInput0::put(MCWalkerConfiguration& W){

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
HDFWalkerInput0::put(MCWalkerConfiguration& W, int ic)
{
  if(Counter<0)
    return false;
  int selected = ic;
  if(ic<0)
  {
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
  if(curWalker)
  {
    LOGMSG("Adding " << nwt << " walkers to " << curWalker)
    W.createWalkers(nwt);
  }
  else
  {
    W.resize(nwt,nptcl);
  }
  //assign configurations to W
  int iw=0;
  MCWalkerConfiguration::iterator it = W.begin()+curWalker;
  MCWalkerConfiguration::iterator it_end = W.end();
  while(it != it_end)
  {
    copy(Pos_temp[iw],Pos_temp[iw+1], (*it)->R.begin());
    ++it;
    ++iw;
  }
  return true;
}

//bool HDFWalkerInput0::append(MCWalkerConfiguration& W){
//
//  if(Counter<0) return false;
//
//  typedef MCWalkerConfiguration::PosType PosType;
//  typedef Matrix<PosType>  PosContainer_t;
//  PosContainer_t Pos_temp;
//  int nw_in=0,curConfig=FirstSet;
//  while(curConfig<LastSet) {
//    //open the group
//    char GrpName[128];
//    sprintf(GrpName,"config%04d",curConfig);
//    hid_t group_id = H5Gopen(h_config,GrpName);
//    HDFAttribIO<PosContainer_t> Pos_in(Pos_temp);
//    //read the dataset
//    Pos_in.read(group_id,"coord");
//    //close the group
//    H5Gclose(group_id);
//    /*check to see if the number of walkers and particles is  consistent with W */
//    int nptcl = Pos_temp.cols();
//    int nwt = Pos_temp.rows();
//    int curWalker = W.getActiveWalkers();
//    if(curWalker) {
//      W.createWalkers(nwt);
//    } else {
//      W.resize(nwt,nptcl);
//    }
//    MCWalkerConfiguration::iterator it = W.begin()+curWalker;
//    for(int iw=0; iw<nwt; iw++) {
//      //std::copy(Post_temp[iw],Post_temp[iw+1], (*it)->R.begin());
//      for(int iat=0; iat < nptcl; iat++){
//        (*it)->R(iat) = Pos_temp(iw,iat);
//      }
//      ++it;
//    }
//    nw_in += nwt;
//    curConfig++;
//  }
//
//  LOGMSG("Total " << nw_in << " walkers are loaded using " << LastSet-FirstSet << " frames.")
//  return true;
//}

bool
HDFWalkerInput0::append(MCWalkerConfiguration& W, int blocks)
{
  if(Counter<0)
    return false;
  //if(nwalkers<0) return put(W,-1);
  typedef MCWalkerConfiguration::PosType PosType;
  typedef Matrix<PosType>  PosContainer_t;
  PosContainer_t Pos_temp;
  int nw_in=0;
  int firstConf=std::max(0,NumSets-blocks);
  if(blocks<0)
    firstConf=0;
  for(int iconf=firstConf; iconf<NumSets; iconf++)
  {
    //open the group
    char GrpName[128];
    sprintf(GrpName,"config%04d",iconf);
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
    if(nptcl != W.getParticleNum())
    {
      W.resize(nwt,nptcl);
    }
    else
    {
      curWalker=W.getActiveWalkers();
      W.createWalkers(nwt);
    }
    MCWalkerConfiguration::iterator it = W.begin()+curWalker;
    for(int iw=0; iw<nwt; iw++)
    {
      //std::copy(Post_temp[iw],Post_temp[iw+1], (*it)->R.begin());
      for(int iat=0; iat < nptcl; iat++)
      {
        (*it)->R(iat) = Pos_temp(iw,iat);
      }
      ++it;
    }
    nw_in += nwt;
  }
  LOGMSG("Total " << nw_in << " walkers are loaded using " << NumSets-firstConf << " blocks.")
  return true;
}

void  HDFWalkerInput0::getRandomState(bool restart)
{
  if(restart)
  {
    hid_t h_random = H5Gopen(h_file,"random_state");
    if(h_random>-1)
    {
      LOGMSG("Reading the state of the random number generator from the configuration file")
      HDFAttribIO<RandomGenerator_t> r(Random);
      r.read(h_random,"dummy");
      //Random.read(h_random);
      //H5Gclose(h_random);
    }
  }
}
#else
HDFWalkerInput0::HDFWalkerInput0(const std::string& aroot, int ipart, int nparts):
  Counter(0), NumSets(0)
{
  ERRORMSG("HDF5 is disabled")
}

/** Destructor closes the HDF5 file and main group. */
HDFWalkerInput0::~HDFWalkerInput0()
{
  ERRORMSG("HDF5 is disabled")
}

bool
HDFWalkerInput0::put(MCWalkerConfiguration& W, int ic)
{
  ERRORMSG("HDF5 is disabled")
  return false;
}

bool
HDFWalkerInput0::append(MCWalkerConfiguration& W, int blocks)
{
  ERRORMSG("HDF5 is disabled")
  return false;
}

void  HDFWalkerInput0::getRandomState(bool restart)
{
  ERRORMSG("HDF5 is disabled")
}
#endif

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
    
    



#include "Particle/HDFWalkerInput_0_0.h"
#include "Particle/MCWalkerConfiguration.h"
#include "OhmmsPETE/OhmmsVector.h"
#include "Particle/HDFParticleAttrib.h"
#include "Numerics/HDFNumericAttrib.h"
#include "Utilities/RandomGenerator.h"
#include "OhmmsData/FileUtility.h"
namespace qmcplusplus
{

HDFWalkerInput_0_0::HDFWalkerInput_0_0(MCWalkerConfiguration& W, const std::string &aroot):
  targetW(W), version(0,0)
{
  FileName = aroot;
  std::string ext=getExtension(aroot);
  if(ext != "h5")
    //if the filename does not h5 extension, add the extension
  {
    FileName.append(".config.h5");
  }
}

HDFWalkerInput_0_0::~HDFWalkerInput_0_0() { }

bool  HDFWalkerInput_0_0::put(xmlNodePtr cur)
{
  hid_t h_file =  H5Fopen(FileName.c_str(),H5F_ACC_RDWR,H5P_DEFAULT);
  hid_t h_config = H5Gopen(h_file,"config_collection");
  if(h_config<0)
  {
    app_error() << " HDFWalkerInput_0_0::put  config_collection is not found in " << FileName << "." << std::endl;
    H5Fclose(h_file);
    return false;
  }
  int NumSets=0;
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
  if(!NumSets)
  {
    app_error() << " HDFWalkerInput_0_0::put  " << FileName << " does not contain walkers!" << std::endl;
    H5Gclose(h_config);
    H5Gclose(h_file);
    return false;
  }
  //select the last block
  int selected=NumSets-1;
  typedef MCWalkerConfiguration::PosType PosType;
  typedef MCWalkerConfiguration::PropertyContainer_t ProtertyContainer_t;
  typedef Matrix<TinyVector<double, OHMMS_DIM> >  PosContainer_t;
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
  //close groups
  H5Gclose(group_id);
  H5Gclose(h_config);
  H5Fclose(h_file);
  /*check to see if the number of walkers and particles is  consistent with W */
  int nptcl = Pos_temp.cols();
  nwt = Pos_temp.rows();
  int curWalker = targetW.getActiveWalkers();
  if(curWalker)
  {
    app_log() << "HDFWalkerInput_0_0::put Adding " << nwt << " walkers to " << curWalker << std::endl;
    targetW.createWalkers(nwt);
  }
  else
  {
    app_log() << "HDFWalkerInput_0_0::put Creating " << nwt << " walkers." << std::endl;
    targetW.resize(nwt,nptcl);
  }
  //assign configurations to W
  int iw=0;
  MCWalkerConfiguration::iterator it = targetW.begin()+curWalker;
  MCWalkerConfiguration::iterator it_end = targetW.end();
  while(it != it_end)
  {
    std::copy(Pos_temp[iw],Pos_temp[iw+1], (*it)->R.begin());
    ++it;
    ++iw;
  }
  return true;
}

//bool
//HDFWalkerInput_0_0::append(MCWalkerConfiguration& W, int blocks){
//
//  if(Counter<0) return false;
//
//  //if(nwalkers<0) return put(W,-1);
//
//  typedef MCWalkerConfiguration::PosType PosType;
//  typedef Matrix<PosType>  PosContainer_t;
//  PosContainer_t Pos_temp;
//
//  int nw_in=0;
//  int firstConf=std::max(0,NumSets-blocks);
//  if(blocks<0) firstConf=0;
//
//  for(int iconf=firstConf; iconf<NumSets; iconf++) {
//    //open the group
//    char GrpName[128];
//    sprintf(GrpName,"config%04d",iconf);
//    hid_t group_id = H5Gopen(h_config,GrpName);
//    HDFAttribIO<PosContainer_t> Pos_in(Pos_temp);
//    //read the dataset
//    Pos_in.read(group_id,"coord");
//    //close the group
//    H5Gclose(group_id);
//    /*check to see if the number of walkers and particles is  consistent with W */
//    int nptcl = Pos_temp.cols();
//    int nwt = Pos_temp.rows();
//    int curWalker=0;
//    if(nptcl != W.getParticleNum()) {
//      W.resize(nwt,nptcl);
//    } else {
//      curWalker=W.getActiveWalkers();
//      W.createWalkers(nwt);
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
//  }
//
//  LOGMSG("Total " << nw_in << " walkers are loaded using " << NumSets-firstConf << " blocks.")
//  return true;
//}
//void  HDFWalkerInput_0_0::getRandomState(bool restart){
//  if(restart) {
//    hid_t h_random = H5Gopen(h_file,"random_state");
//    if(h_random>-1) {
//      LOGMSG("Reading the state of the random number generator from the configuration file")
//      HDFAttribIO<RandomGenerator_t> r(Random);
//      r.read(h_random,"dummy");
//      //Random.read(h_random);
//      //H5Gclose(h_random);
//    }
//  }
//}
}

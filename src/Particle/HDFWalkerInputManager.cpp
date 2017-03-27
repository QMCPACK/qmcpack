//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Cynthia Gu, zg1@ornl.gov, Oak Ridge National Laboratory
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



#include "Particle/HDFWalkerInputManager.h"
#include "OhmmsData/AttributeSet.h"
#if defined(HAVE_LIBHDF5)
#include "Particle/HDFWalkerInput_0_0.h"
#include "Particle/HDFWalkerInput_0_4.h"
#endif
#include "Message/Communicate.h"
#include "HDFVersion.h"
#ifdef HAVE_ADIOS
#include "ADIOS/ADIOS_config.h"
#endif

namespace qmcplusplus
{

HDFWalkerInputManager::HDFWalkerInputManager(MCWalkerConfiguration& w, Communicate* c):
  targetW(w),myComm(c)
{
}

HDFWalkerInputManager::~HDFWalkerInputManager()
{
}

#if defined(HAVE_LIBHDF5) || defined(HAVE_ADIOS)
bool HDFWalkerInputManager::put(xmlNodePtr cur)
{
  //reference revision number
  HDFVersion start_version(0,4);
  //current node
  int pid=myComm->rank();
  std::string froot("0"), cfile("0");
  //string  target("e"), collect("no");
  int anode=-1, nblocks=1, nprocs=1;
  HDFVersion in_version(0,4); //set to be old version
  OhmmsAttributeSet pAttrib;
  pAttrib.add(cfile,"href");
  pAttrib.add(cfile,"file");
  pAttrib.add(froot,"fileroot");
  pAttrib.add(anode,"node");
  pAttrib.add(nprocs,"nprocs");
  //pAttrib.add(collect,"collected");
  pAttrib.add(in_version,"version");
  pAttrib.put(cur);
  bool success=false;
  if(in_version>=start_version)
  {
    HDFWalkerInput_0_4 win(targetW,myComm,in_version);
    success= win.put(cur);
    cfile=win.FileName;
  }
  else
  {
    //missing version or old file
    if(froot[0] != '0')//use nprocs
    {
      anode=pid;
      if(nprocs==1)
        cfile=froot;
      else
      {
        char *h5name=new char[froot.size()+10];
        sprintf(h5name,"%s.p%03d",froot.c_str(),pid);
        cfile=h5name;
        delete [] h5name;
      }
    }
    int pid_target= (anode<0)? pid:anode;
    if(pid_target == pid && cfile[0] != '0')
    {
      HDFWalkerInput_0_0 win(targetW,cfile);
      success= win.put(cur);
    }
  }
  if(success)
    CurrentFileRoot = cfile;
  return success;
}
#else
bool HDFWalkerInputManager::put(xmlNodePtr cur)
{
  return false;
}
#endif

void HDFWalkerInputManager::rewind(const std::string& h5root, int blocks)
{
//   HDFWalkerInputCollect WO(h5root);
//   WO.rewind(targetW,blocks);
}
}

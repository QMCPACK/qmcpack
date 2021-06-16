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


#include "HDFWalkerInputManager.h"
#include "OhmmsData/AttributeSet.h"
#if defined(HAVE_LIBHDF5)
#include "Particle/HDFWalkerInput_0_4.h"
#endif
#include "Message/Communicate.h"
#include "hdf/HDFVersion.h"

namespace qmcplusplus
{
HDFWalkerInputManager::HDFWalkerInputManager(MCWalkerConfiguration& w, Communicate* c) : targetW(w), myComm(c) {}

HDFWalkerInputManager::~HDFWalkerInputManager() {}

#if defined(HAVE_LIBHDF5)
bool HDFWalkerInputManager::put(xmlNodePtr cur)
{
  //reference revision number
  HDFVersion start_version(0, 4);
  //current node
  int pid = myComm->rank();
  std::string froot("0"), cfile("0");
  //string  target("e"), collect("no");
  int anode = -1, nprocs = 1;
  HDFVersion in_version(0, 4); //set to be old version
  OhmmsAttributeSet pAttrib;
  pAttrib.add(cfile, "href");
  pAttrib.add(cfile, "file");
  pAttrib.add(froot, "fileroot");
  pAttrib.add(anode, "node");
  pAttrib.add(nprocs, "nprocs");
  //pAttrib.add(collect,"collected");
  pAttrib.add(in_version, "version");
  pAttrib.put(cur);
  bool success = false;
  if (in_version >= start_version)
  {
    HDFWalkerInput_0_4 win(targetW, myComm, in_version);
    success = win.put(cur);
    cfile   = win.FileName;
  }
  else
    myComm->barrier_and_abort("Outdated restart file!");
  if (success)
    CurrentFileRoot = cfile;
  return success;
}
#else
bool HDFWalkerInputManager::put(xmlNodePtr cur) { return false; }
#endif

void HDFWalkerInputManager::rewind(const std::string& h5root, int blocks)
{
  //   HDFWalkerInputCollect WO(h5root);
  //   WO.rewind(targetW,blocks);
}
} // namespace qmcplusplus

//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



#ifdef HAVE_ADIOS
#include <Configuration.h>
#include <adios_read.h>
#include "Particle/Walker.h"
#include "Particle/ParticleSet.h"
#include "Particle/MCWalkerConfiguration.h"

namespace qmcplusplus
{

class AdiosWalkerInput
{
public :

  std::string FileRoot;
  typedef ParticleSet::Walker_t::ParticlePos_t R_t;
  int particle_num;
  MCWalkerConfiguration targetW;
  Communicate* myComm;

  AdiosWalkerInput(MCWalkerConfiguration& w, Communicate* c);
  ~AdiosWalkerInput();

  bool put(std::vector<xmlNodePtr>& wset);

  void read(int nprocs, ADIOS_FILE* file_handle, int& walker_win, std::vector<int>& nw);
  void append_walkers(R_t& walker_buff,int, int);
  void setMCWalker(std::vector<int> nw);

  std::string getFileRoot();
};
}

#endif



















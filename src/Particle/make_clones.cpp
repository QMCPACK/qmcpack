//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



/** @file make_clones.cpp
 * @brief implement cloning functions for ParticleSet and MCWalkerConfiguration
 *
 * Unlike normal copy operation, e.g., copy constructor, the name and ObjectTag of the clones are
 * inherited from their parent.
 */
#include "Particle/ParticleSet.h"
#include "Particle/MCWalkerConfiguration.h"
#include "Particle/DistanceTableData.h"

namespace qmcplusplus
{

void ParticleSet::make_clones(int n)
{
  //if(myClones.size()>n) return;
  if(myClones.empty())
  {
    myClones.resize(n,0);
//#pragma omp parallel for
    for(int i=0; i<n; ++i)
      if(i)
        myClones[i]=new ParticleSet(*this);
  }
  else
  {
    n -= myClones.size();
    while(n>0)
    {
      myClones.push_back(new ParticleSet(*this));
      n--;
    }
  }
  this->ThreadID=0;
  for(int i=1; i<myClones.size(); ++i)
  {
    myClones[i]->ThreadID=i;
    //myClones[i]->setName(this->getName());
  }
}

void ParticleSet::reset_clones()
{
  //if(myClones.empty()) return;
  int nt=myClones.size();
  for(int t=1;t<DistTables.size(); ++t)
  {
    const ParticleSet& other=DistTables[t]->origin();
    if(nt == other.clones_size())
    {
      for(int tid=1; tid<nt; ++tid)
      {
        myClones[tid]->DistTables[t]->reset(other.get_clone(tid));
      }
    }
  }

  //for(int tid=0; tid<nt; ++tid)
  //{
  //  std::cout << "Thread ID = " << get_clone(tid)->getName() << std::endl;
  //  for(int t=0; t<DistTables.size(); ++t)
  //    std::cout << get_clone(tid)->DistTables[t]->origin().getName() << std::endl;
  //}
}

void MCWalkerConfiguration::make_clones(int n)
{
  if(myClones.empty())
  {
    myClones.resize(n,0);
//#pragma omp parallel for
    for(int i=0; i<n; ++i)
      if(i)
        myClones[i]=new MCWalkerConfiguration(*this);
  }
  else
  {
    n -= myClones.size();
    while(n>0)
    {
      myClones.push_back(new MCWalkerConfiguration(*this));
      n--;
    }
  }
  this->ThreadID=0;
  for(int i=1; i<myClones.size(); ++i)
  {
    myClones[i]->ThreadID=i;
    myClones[i]->setName(this->getName());
  }
}
}


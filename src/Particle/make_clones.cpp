//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
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
  //  cout << "Thread ID = " << get_clone(tid)->getName() << endl;
  //  for(int t=0; t<DistTables.size(); ++t)
  //    cout << get_clone(tid)->DistTables[t]->origin().getName() << endl;
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

/***************************************************************************
 * $RCSfile$   $Author: jeongnim.kim $
 * $Revision: 5325 $   $Date: 2011-08-03 16:53:31 -0400 (Wed, 03 Aug 2011) $
 * $Id: ParticleSet.cpp 5325 2011-08-03 20:53:31Z jeongnim.kim $
 ***************************************************************************/

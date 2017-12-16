//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


/** @file VirtualParticleSet.cpp
 * A proxy class to the quantum ParticleSet
 */

#include <Configuration.h>
#include <Particle/VirtualParticleSet.h>
#include <Particle/DistanceTableData.h>
#include <Particle/DistanceTable.h>

namespace qmcplusplus
{

  VirtualParticleSet::VirtualParticleSet(const ParticleSet& p, int nptcl): refPtcl(p)
  {
    setName("virtual");

    //initialize local data structure
    Lattice = p.Lattice;
    R.resize(nptcl);
    RSoA.resize(nptcl);

    //create distancetables
    if(refPtcl.DistTables.size())
    {
      DistTables.resize(refPtcl.DistTables.size());
      for(int i=0; i<DistTables.size(); ++i)
      {
        DistTables[i]=createDistanceTable(refPtcl.DistTables[i]->origin(),*this, refPtcl.DistTables[0]->DTType);
        DistTables[i]->ID=i;
      }
    }
  }

  /// move virtual particles to new postions and update distance tables
  void VirtualParticleSet::makeMoves(int iat, const ParticlePos_t& vitualPos)
  {
    refID=iat;
    R=vitualPos;
    RSoA.copyIn(R);
    for (int i=0; i<DistTables.size(); i++)
      DistTables[i]->evaluate(*this);
  }

}

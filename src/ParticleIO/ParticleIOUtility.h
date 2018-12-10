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
    
    


#ifndef OHMMS_PARTICLE_INPUTOUTPUT_UTILITY_H
#define OHMMS_PARTICLE_INPUTOUTPUT_UTILITY_H

#include "OhmmsData/OhmmsElementBase.h"
#include "Particle/ParticleSet.h"

namespace qmcplusplus
{

//class LatticeParser {
//  ParticleSet::ParticleLayout_t& ref_;
//  public:
//  LatticeParser(ParticleSet::ParticleLayout_t& lat): ref_(lat){ }
//  //PARSEDOM
//};


//class LatticeXMLWriter {

//  ParticleSet::ParticleLayout_t& ref_;
//  public:

//  LatticeXMLWriter(ParticleSet::ParticleLayout_t& lat): ref_(lat){ }
//  bool get(std::ostream& ) const;
//};

//class ParticleParser {

//  ParticleSet& ref_;

//public:

//  ParticleParser(ParticleSet& aptcl):ref_(aptcl){ }
//  ///reading from a file
//  ///bool put(const char*);
//  //PARSEDOM
//};



template<class GI>
void PartitionGrid(ParticleSet::ParticleLayout_t& lattice, GI& grid)
{
  const int ndim = ParticleSet::ParticleLayout_t::DIM;
  ////////////////////////////////////////////////////////////////
  // grid partitions are created
  ////////////////////////////////////////////////////////////////
  int mppgrid_id = lattice.makeGrid(&grid[0][0],-1);
  int ompgrid_id = lattice.makeGrid(&grid[1][0],mppgrid_id);
  for(int idim=0; idim<ndim; idim++)
  {
    if(grid[1][idim] > grid[2][idim])
      grid[1][idim] = grid[2][idim];
    grid[2][idim] /= grid[1][idim];
  }
  lattice.makeGrid(&grid[2][0],ompgrid_id);
  for(int idim=0; idim<ndim; idim++)
  {
    lattice.Grid[idim] = grid[2][idim]*grid[1][idim];
  }
}

/** expand a particleset including lattice */
void expandSuperCell(ParticleSet& ref, Tensor<int,OHMMS_DIM>& tmat);

}
#endif


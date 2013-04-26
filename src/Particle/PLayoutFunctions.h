//////////////////////////////////////////////////////////////////
// (c) Copyright 2004- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
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
#include "Configuration.h"
#include "Particle/FastParticleOperators.h"
namespace OHMMS
{

struct PLayoutFunc : public PtclOnLatticeTraits
{

  static inline void
  applyBC(const ParticleLayout_t& Lattice, const ParticlePos_t& pin, ParticlePos_t& pout)
  {
    const bool orthogonal = ParticleLayout_t::IsOrthogonal;
    const int dim = OHMMS_DIM;
    int mode = pin.getUnit()*2+pout.getUnit();
    int last = pin.size();
    switch(mode)
    {
    case(0):
      ApplyBConds<ParticlePos_t,Tensor_t,orthogonal>::Cart2Cart(pin,Lattice.G,Lattice.R,pout,0,last);
      break;
    case(1):
      ApplyBConds<ParticlePos_t,Tensor_t,orthogonal>::Cart2Unit(pin,Lattice.G,pout,0,last);
      break;
    case(2):
      ApplyBConds<ParticlePos_t,Tensor_t,orthogonal>::Unit2Cart(pin,Lattice.R,pout,0,last);
      break;
    case(3):
      ApplyBConds<ParticlePos_t,Tensor_t,orthogonal>::Unit2Unit(pin,pout,0,last);
      break;
    }
  }

  static inline void
  applyBC(const ParticleLayout_t& Lattice, ParticlePos_t& pos)
  {
    const bool orthogonal = ParticleLayout_t::IsOrthogonal;
    const int dim = OHMMS_DIM;
    int last = pos.size();
    if(pos.getUnit()==PosUnit::LatticeUnit)
    {
      ApplyBConds<ParticlePos_t,Tensor_t,orthogonal>::Unit2Unit(pos,0,last);
    }
    else
    {
      ApplyBConds<ParticlePos_t,Tensor_t,orthogonal>::Cart2Cart(pos,Lattice.G,Lattice.R,0,last);
    }
  }
};
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/


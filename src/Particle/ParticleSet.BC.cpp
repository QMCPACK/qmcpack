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
/**@file ParticleSet.BC.cpp
 * @brief definition of functions controlling Boundary Conditions
 */
#include "Particle/ParticleSet.h"
#include "Particle/FastParticleOperators.h"
#include "Utilities/OhmmsInfo.h"
#include "Message/OpenMP.h"
#include "LongRange/StructFact.h"

namespace qmcplusplus {

  /** Creating StructureFactor
   * 
   * Currently testing only 1 component for PBCs.
   */
  void ParticleSet::createSK() {
    convert2Cart(R); //make sure that R is in Cartesian coordinates
    if(Lattice.BoxBConds[0] && SK == 0){
      LOGMSG("\n  Creating Structure Factor for periodic systems.")
      Lattice.SetLRCutoffs();
      SK = new StructFact(*this,Lattice.LR_kc);

      Lattice.print(app_log());
      //This uses the copy constructor to avoid recomputing the data.
      //SKOld = new StructFact(*SK);
    }
  }

  void ParticleSet::convert(const ParticlePos_t& pin, ParticlePos_t& pout){

    if(pin.getUnit() == pout.getUnit())   {
      pout = pin;
      return;
    }
    if(pin.getUnit() == PosUnit::LatticeUnit) { //convert to CartesianUnit
      ConvertPosUnit<ParticlePos_t,Tensor_t,OHMMS_ORTHO>::apply(pin,Lattice.R,pout,0,pin.size());
    } else { //convert to LatticeUnit
      ConvertPosUnit<ParticlePos_t,Tensor_t,OHMMS_ORTHO>::apply(pin,Lattice.G,pout,0,pin.size());
    }
  }

  void ParticleSet::convert2Unit(const ParticlePos_t& pin, ParticlePos_t& pout){
    pout.setUnit(PosUnit::LatticeUnit);
    if(pin.getUnit() == PosUnit::LatticeUnit) 
      pout = pin;
    else 
      ConvertPosUnit<ParticlePos_t,Tensor_t,OHMMS_ORTHO>::apply(pin,Lattice.G,pout,0,pin.size());
  }

  void ParticleSet::convert2Cart(const ParticlePos_t& pin, ParticlePos_t& pout) {
    pout.setUnit(PosUnit::CartesianUnit);
    if(pin.getUnit() == PosUnit::CartesianUnit) 
      pout = pin;
    else 
      ConvertPosUnit<ParticlePos_t,Tensor_t,OHMMS_ORTHO>::apply(pin,Lattice.R,pout,0,pin.size());
  }

  void ParticleSet::convert2Unit(ParticlePos_t& pinout) {
    if(pinout.getUnit() == PosUnit::LatticeUnit) 
      return;
    else {
      pinout.setUnit(PosUnit::LatticeUnit);
      ConvertPosUnit<ParticlePos_t,Tensor_t,OHMMS_ORTHO>::apply(pinout,Lattice.G,0,pinout.size());
    }
  }

  void ParticleSet::convert2Cart(ParticlePos_t& pinout) {
    if(pinout.getUnit() == PosUnit::CartesianUnit) 
      return;
    else {
      pinout.setUnit(PosUnit::CartesianUnit);
      ConvertPosUnit<ParticlePos_t,Tensor_t,OHMMS_ORTHO>::apply(pinout,Lattice.R,0,pinout.size());
    }
  }

  void ParticleSet::applyBC(const ParticlePos_t& pin, ParticlePos_t& pout) {
    applyBC(pin,pout,0,pin.size());
  }

  void ParticleSet::applyBC(const ParticlePos_t& pin, ParticlePos_t& pout, int first, int last) {
    const bool orthogonal = ParticleLayout_t::IsOrthogonal;
    const int dim = OHMMS_DIM;
    int mode = pin.getUnit()*2+pout.getUnit();
    switch(mode) {
    case(0):
      ApplyBConds<ParticlePos_t,Tensor_t,orthogonal>::Cart2Cart(pin,Lattice.G,Lattice.R,pout,first,last);
      break;
    case(1):
      ApplyBConds<ParticlePos_t,Tensor_t,orthogonal>::Cart2Unit(pin,Lattice.G,pout,first,last);
      break;
    case(2):
      ApplyBConds<ParticlePos_t,Tensor_t,orthogonal>::Unit2Cart(pin,Lattice.R,pout,first,last);
      break;
    case(3):
      ApplyBConds<ParticlePos_t,Tensor_t,orthogonal>::Unit2Unit(pin,pout,first,last);
      break;
    }
  }

  void ParticleSet::applyBC(ParticlePos_t& pos) {
    const bool orthogonal = ParticleLayout_t::IsOrthogonal;
    const int dim = OHMMS_DIM;
    if(pos.getUnit()==PosUnit::LatticeUnit) {
      ApplyBConds<ParticlePos_t,Tensor_t,orthogonal>::Unit2Unit(pos,0,LocalNum);
    } else {
      ApplyBConds<ParticlePos_t,Tensor_t,orthogonal>::Cart2Cart(pos,Lattice.G,Lattice.R,0,LocalNum);
    }
  }
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/


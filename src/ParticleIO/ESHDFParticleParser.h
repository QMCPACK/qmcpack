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
    
    


#ifndef ESHDF_PARTICLE_PARSER_H
#define ESHDF_PARTICLE_PARSER_H

#include "Particle/ParticleSet.h"
#include "OhmmsData/HDFAttribIO.h"

class Communicate;

namespace qmcplusplus
{

struct ESHDFElectronsParser
{

  ESHDFElectronsParser(ParticleSet& aptcl, hid_t h=-1, Communicate* c=0);

  bool put(xmlNodePtr cur);

  ParticleSet& ref_;
  hid_t hfile_id;
  Communicate* myComm;
};

struct ESHDFIonsParser
{

  ESHDFIonsParser(ParticleSet& aptcl, hid_t h=-1, Communicate* c=0);

  bool put(xmlNodePtr cur);
  void readESHDF();
  //expand the ionic systems
  void expand(Tensor<int,3>& tilematrix);

  ParticleSet& ref_;
  hid_t hfile_id;
  Communicate* myComm;
  std::string atomic_number_tag;
  std::string charge_tag;
  std::string mass_tag;
};

}

#endif

/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 757 $   $Date: 2005-10-31 10:10:28 -0600 (Mon, 31 Oct 2005) $
 * $Id: ESHDFParticleParser.h 757 2005-10-31 16:10:28Z jnkim $
 ***************************************************************************/

//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#ifndef QMCPLUSPLUS_BSPLINE_JASTROW_BUILDER_H
#define QMCPLUSPLUS_BSPLINE_JASTROW_BUILDER_H
#include "QMCWaveFunctions/WaveFunctionComponentBuilder.h"

namespace qmcplusplus
{
//forward declarations
class ParticleSet;

struct BsplineJastrowBuilder: public WaveFunctionComponentBuilder
{
  ParticleSet *sourcePtcl;
  // One-body constructor
  BsplineJastrowBuilder(ParticleSet& target, TrialWaveFunction& psi,
                        ParticleSet& source) :
    WaveFunctionComponentBuilder(target,psi), sourcePtcl(&source)
  {
    ClassName="BsplineJastrowBuilder";
  }
  // Two-body constructor
  BsplineJastrowBuilder(ParticleSet& target, TrialWaveFunction& psi) :
    WaveFunctionComponentBuilder(target,psi), sourcePtcl(NULL)
  {
    ClassName="BsplineJastrowBuilder";
  }

  /** create one-body Jastrow
   * @tparm OBJT one-body jastrow orbital class
   */
  template<typename OBJT, typename DOBJT>
  bool createOneBodyJastrow(xmlNodePtr cur);

  bool put(xmlNodePtr cur);
};

}
#endif

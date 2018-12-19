//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    

#ifndef OPTIMIZABLE_SPO_BUILDER_H
#define OPTIMIZABLE_SPO_BUILDER_H

#include "QMCWaveFunctions/SPOSetBuilder.h"
#include "QMCWaveFunctions/OptimizableSPOSet.h"


namespace qmcplusplus
{
class OptimizableSPOBuilder : public SPOSetBuilder
{
protected:
  typedef std::map<std::string,ParticleSet*> PtclPoolType;
  typedef std::map<std::string,SPOSet*>  SPOPoolType;
  ParticleSet *targetPtcl;
public:
  OptimizableSPOBuilder(ParticleSet& p, PtclPoolType& psets,
                        Communicate *comm,
                        xmlNodePtr cur);

  /** initialize the Antisymmetric wave function for electrons
   *@param cur the current xml node
   */
  SPOSet* createSPOSetFromXML(xmlNodePtr cur);
  //    SPOSet* createSPOSetFromXML(xmlNodePtr cur, SPOPool_t& spo_pool);

};
}

#endif

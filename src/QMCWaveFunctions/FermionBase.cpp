//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#include "QMCWaveFunctions/FermionBase.h"

namespace qmcplusplus
{

///add a new SPOSet to the list of determinants
void FermionBase::addSPO(const std::string& aname, SPOSetBase* sposet)
{
  if (mySPOSet.find(aname) == mySPOSet.end())
  {
    mySPOSet[aname] = sposet;
    sposet->objectName = aname;
  }
  else
  {
    APP_ABORT(" FermionBase::addSPO(sposet,aname) cannot reuse the " + aname );
  }
}

void FermionBase::resetSPOs(ParticleSet& P)
{
  std::map<std::string, SPOSetBasePtr>::iterator sit(mySPOSet.begin());
  while (sit != mySPOSet.end())
  {
    (*sit).second->resetTargetParticleSet(P);
    ++sit;
  }
}

void FermionBase::cloneSPOs(const spo_set_type& spos, ParticleSet& tqp)
{
  spo_set_type::const_iterator sit(spos.begin());
  while (sit != spos.end())
  {
    SPOSetBasePtr spo = (*sit).second;
    addSPO((*sit).first, spo->makeClone());
    ++sit;
  }
}

void FermionBase::copySPOs(spo_set_type& spos)
{
  spo_set_type::const_iterator sit(spos.begin());
  while (sit != spos.end())
  {
    addSPO((*sit).first, (*sit).second);
    ++sit;
  }
}

SPOSetBasePtr FermionBase::getSPO(const std::string& aname)
{
  spo_set_type::iterator sit=mySPOSet.find(aname);
  if(sit == mySPOSet.end())
  {
    APP_ABORT("FermionBase::getSPO failed. Missing SPOSet with " + aname);
  }
  return (*sit).second;
}


}




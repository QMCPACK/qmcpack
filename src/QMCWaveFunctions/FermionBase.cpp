//////////////////////////////////////////////////////////////////
// (c) Copyright 2006-  by Jeongnim Kim
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
#include "QMCWaveFunctions/FermionBase.h"

namespace qmcplusplus
{

///add a new SPOSet to the list of determinants
void FermionBase::addSPO(const string& aname, SPOSetBase* sposet)
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
  map<string, SPOSetBasePtr>::iterator sit(mySPOSet.begin());
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

SPOSetBasePtr FermionBase::getSPO(const string& aname)
{
  spo_set_type::iterator sit=mySPOSet.find(aname);
  if(sit == mySPOSet.end())
  {
    APP_ABORT("FermionBase::getSPO failed. Missing SPOSet with " + aname);
  }
  return (*sit).second;
}


}


/***************************************************************************
 * $RCSfile$   $Author: jmcminis $
 * $Revision: 5443 $   $Date: 2012-03-21 20:14:56 -0500 (Wed, 21 Mar 2012) $
 * $Id: FermionBase.h 5443 2012-03-22 01:14:56Z jmcminis $
 ***************************************************************************/


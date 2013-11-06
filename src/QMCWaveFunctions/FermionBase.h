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
#ifndef QMCPLUSPLUS_FERMIONBASE_H
#define QMCPLUSPLUS_FERMIONBASE_H

#include "QMCWaveFunctions/SPOSetBase.h"
#include "QMCWaveFunctions/BasisSetBase.h"

namespace qmcplusplus
{

/** base class for Fermion
 */
struct FermionBase
{
  typedef BasisSetBuilder* BasisSetBuilderPtr;

  typedef map<string,SPOSetBasePtr> spo_set_type;
  
  spo_set_type mySPOSet;
  
  /** add a new SPOSet
   * @param aname name of a SPOset
   * @param sposet to be added to the list
   */
  void addSPO(const string& aname, SPOSetBase* sposet);

  void addSPO(SPOSetBase* sposet);

  /** reset SPOs
   * @param P target ParticleSet
   */
  void resetSPOs(ParticleSet& P);
  
  /** copy SPOsets of a
   */
  inline void copySPOs(FermionBase& a)
  {
    copySPOs(a.mySPOSet);
  }
  /** copy spos
   */
  void copySPOs(spo_set_type& spos);

  /** clone SPOs
   * @param spos a set of SPOs
   * @param tqp target quantum particleset
   */
  void cloneSPOs(const spo_set_type& spos, ParticleSet& tqp);

  /** get SPO with phi.objectName
   * @param phi whose name is used to find a SPO
   */
  inline SPOSetBasePtr getSPO(SPOSetBasePtr phi)
  {
    return getSPO(phi->objectName);
  }

  /** get SPO with aname
   */
  SPOSetBasePtr getSPO(const string& aname);

  /** return true if this has the SPOSet with aname
   */
  inline bool found(const string& aname)
  {
    return (mySPOSet.find(aname) != mySPOSet.end());
  }
};

}
#endif

/***************************************************************************
 * $RCSfile$   $Author: jmcminis $
 * $Revision: 5443 $   $Date: 2012-03-21 20:14:56 -0500 (Wed, 21 Mar 2012) $
 * $Id: FermionBase.h 5443 2012-03-22 01:14:56Z jmcminis $
 ***************************************************************************/


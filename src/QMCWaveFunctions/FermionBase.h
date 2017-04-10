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

  typedef std::map<std::string,SPOSetBasePtr> spo_set_type;
  
  spo_set_type mySPOSet;
  
  /** add a new SPOSet
   * @param aname name of a SPOset
   * @param sposet to be added to the list
   */
  void addSPO(const std::string& aname, SPOSetBase* sposet);

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
  SPOSetBasePtr getSPO(const std::string& aname);

  /** return true if this has the SPOSet with aname
   */
  inline bool found(const std::string& aname)
  {
    return (mySPOSet.find(aname) != mySPOSet.end());
  }
};

}
#endif



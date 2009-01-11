//////////////////////////////////////////////////////////////////
// (c) Copyright 2009-  by Jeongnim Kim
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
/**@file QMCHamiltonianBase.cpp
 *@brief Definition of QMCHamiltonianBase
 */
#include <QMCHamiltonians/QMCHamiltonianBase.h>
#include <QMCHamiltonians/QMCHamiltonian.h>

namespace qmcplusplus 
{
  void QMCHamiltonianBase::add2Hamiltonian(ParticleSet& qp, TrialWaveFunction& psi
      ,QMCHamiltonian& targetH)
  {
    QMCHamiltonianBase* myclone=makeClone(qp,psi);
    if(myclone)
      targetH.addOperator(myclone,myName,UpdateMode[PHYSICAL]);
  }

  void QMCHamiltonianBase::registerObservables(vector<observable_helper*>& h5list
      , hid_t gid)  const
  {
    int loc=h5list.size();
    h5list.push_back(new observable_helper(myName));
    std::vector<int> onedim(1,1);
    h5list[loc]->set_dimensions(onedim,myIndex);
    h5list[loc]->open(gid);
  }
}
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 3437 $   $Date: 2008-12-18 13:44:04 -0600 (Thu, 18 Dec 2008) $
 * $Id: QMCHamiltonianBase.cpp 3437 2008-12-18 19:44:04Z jnkim $ 
 ***************************************************************************/


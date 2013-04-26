//////////////////////////////////////////////////////////////////
// (c) Copyright 2008-  by Jeremy McMinis and Jeongnim Kim
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
//   Department of Physics, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
/** @file ForwardWalking.h
 * @brief Declarations of ForwardWalking
 */
#ifndef QMCPLUSPLUS_FORWARDWALKING_H
#define QMCPLUSPLUS_FORWARDWALKING_H
#include "QMCHamiltonians/QMCHamiltonianBase.h"

namespace qmcplusplus
{

class QMCHamiltonian;

struct ForwardWalking: public QMCHamiltonianBase
{
  vector<int> Hindices;
  vector<int> Pindices;
  vector<vector<int> > walkerLengths;
  vector<double> Values;
  vector<string> Names;
  int blockT,nObservables,nValues,FirstHamiltonian;
  double count;

  /** constructor
   */
  ForwardWalking()
  {
    UpdateMode.set(OPTIMIZABLE,1);
  }

  ///destructor
  ~ForwardWalking() { }

  void resetTargetParticleSet(ParticleSet& P) { }

  inline Return_t rejectedMove(ParticleSet& P)
  {
    for (int i=0; i<nObservables; i++)
    {
      int lastindex = tWalker->PHindex[Pindices[i]]-1;
      if (lastindex<0)
        lastindex +=walkerLengths[i][2];
      tWalker->addPropertyHistoryPoint(Pindices[i],  tWalker->PropertyHistory[Pindices[i]][lastindex]  );
    }
    calculate(P);
    return 0.0;
  }

  inline Return_t calculate(ParticleSet& P)
  {
    vector<Return_t>::iterator Vit=Values.begin();
    for(int i=0; i<nObservables; i++)
    {
      int j=0;
      int FWindex = tWalker->PHindex[Pindices[i]]-1;
      while (j<walkerLengths[i][1])
      {
        FWindex -= walkerLengths[i][0];
        if (FWindex< 0)
          FWindex += walkerLengths[i][2];
        (*Vit) = tWalker->PropertyHistory[Pindices[i]][FWindex];
        j++;
        Vit++;
      }
    }
    std::copy(Values.begin(),Values.end(),tWalker->getPropertyBase()+FirstHamiltonian+myIndex);
    return 0.0;
  }



  inline Return_t evaluate(ParticleSet& P)
  {
    for(int i=0; i<nObservables; i++)
      tWalker->addPropertyHistoryPoint(Pindices[i],  P.PropertyList[Hindices[i]]);
    calculate(P);
    return 0.0;
  }

  inline Return_t evaluate(ParticleSet& P, vector<NonLocalData>& Txy)
  {
    return evaluate(P);
  }

  bool put(xmlNodePtr cur)
  {
    return true;
  }

  ///rename it to avoid conflicts with put
  bool putSpecial(xmlNodePtr cur, QMCHamiltonian& h, ParticleSet& P);

  bool get(std::ostream& os) const
  {
    os << "ForwardWalking";
    return true;
  }

  QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi);

  void addObservables(PropertySetType& plist);

  void addObservables(PropertySetType& plist, BufferType& collectables);

  void setObservables(PropertySetType& plist)
  {
    std::copy(Values.begin(),Values.end(),plist.begin()+myIndex);
  }

  void setParticlePropertyList(PropertySetType& plist, int offset)
  {
    std::copy(Values.begin(),Values.end(),plist.begin()+myIndex+offset);
  }
};
}
#endif

/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1581 $   $Date: 2007-01-04 10:02:14 -0600 (Thu, 04 Jan 2007) $
 * $Id: ForwardWalking.h 1581 2007-01-04 16:02:14Z jnkim $
 ***************************************************************************/

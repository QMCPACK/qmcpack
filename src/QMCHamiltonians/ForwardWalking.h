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
  std::vector<int> Hindices;
  std::vector<int> Pindices;
  std::vector<std::vector<int> > walkerLengths;
  std::vector<RealType> Values;
  std::vector<std::string> Names;
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
    std::vector<RealType>::iterator Vit=Values.begin();
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
    copy(Values.begin(),Values.end(),tWalker->getPropertyBase()+FirstHamiltonian+myIndex);
    return 0.0;
  }



  inline Return_t evaluate(ParticleSet& P)
  {
    for(int i=0; i<nObservables; i++)
      tWalker->addPropertyHistoryPoint(Pindices[i],  P.PropertyList[Hindices[i]]);
    calculate(P);
    return 0.0;
  }

  inline Return_t evaluate(ParticleSet& P, std::vector<NonLocalData>& Txy)
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
    copy(Values.begin(),Values.end(),plist.begin()+myIndex);
  }

  void setParticlePropertyList(PropertySetType& plist, int offset)
  {
    copy(Values.begin(),Values.end(),plist.begin()+myIndex+offset);
  }
};
}
#endif


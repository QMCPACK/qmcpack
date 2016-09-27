//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#ifndef QMCPLUSPLUS_HEEPOTENTIAL_TAIL_H
#define QMCPLUSPLUS_HEEPOTENTIAL_TAIL_H
#include "Particle/ParticleSet.h"
// #include "Particle/WalkerSetRef.h"
// #include "Particle/DistanceTableData.h"
// #include "Particle/DistanceTable.h"
// #include "QMCHamiltonians/QMCHamiltonianBase.h"

///Using Kelvin and Angstrom
namespace qmcplusplus
{
/** @ingroup hamiltonian
 *@brief HFDHE2Potential for the indentical source and target particle sets.
 */
struct HeePotential_tail: public QMCHamiltonianBase
{
  ParticleSet* PtclRef;

  HeePotential_tail(ParticleSet& P): PtclRef(&P)
  {
  }

  ~HeePotential_tail() { }

  void resetTargetParticleSet(ParticleSet& P)
  {
  }

  inline Return_t evaluate(ParticleSet& P)
  {
    return Value;
  }

  inline void set_TC(Return_t TCorr)
  {
    Value = TCorr;
  }

  inline Return_t evaluate(ParticleSet& P, std::vector<NonLocalData>& Txy)
  {
    return evaluate(P);
  }

  /** Do nothing */
  bool put(xmlNodePtr cur)
  {
    return true;
  }

  bool get(std::ostream& os) const
  {
    os << "HeePotentialTailcorr: " << PtclRef->getName();
    return true;
  }

  QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi)
  {
    return 0;
  }

  void addObservables(PropertySetType& plist)
  {
    myIndex=plist.add("Heetail");
  }

//     void setObservables(PropertySetType& plist)
//     {
//       plist[myIndex]=Value;
//     }
//     void setParticlePropertyList(PropertySetType& plist, int offset)
//     {
//       plist[myIndex+offset]=Value;
//     }
};
}
#endif

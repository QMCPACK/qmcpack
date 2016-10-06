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
    
    
#ifndef QMCPLUSPLUS_HFDHE2POTENTIAL_H
#define QMCPLUSPLUS_HFDHE2POTENTIAL_H
#include "Particle/ParticleSet.h"
#include "Particle/WalkerSetRef.h"
#include "QMCHamiltonians/QMCHamiltonianBase.h"

namespace qmcplusplus
{

/** HFDHE2Potential for the indentical source and target particle sets.
 *
 * Tail correction is evaluated by this object and is used by HFDHE2Potential_tail
 */
struct HFDHE2Potential: public QMCHamiltonianBase
{

  Return_t TCValue;
  Return_t kpre_factor;
  Return_t tailcorr,rc,A,alpha,c1,c2,c3,D;

  // epsilon = 3.42016039e-5, rm = 5.607384357
  // C6 = 1.460008056, C8 = 14.22016431, C10 = 187.2033646;
  HFDHE2Potential(ParticleSet& P);

  ~HFDHE2Potential() { }

  void resetTargetParticleSet(ParticleSet& P);

  Return_t evaluate(ParticleSet& P);

  inline Return_t evaluate(ParticleSet& P, std::vector<NonLocalData>& Txy)
  {
    return evaluate(P);
  }

  inline Return_t dampF(Return_t r)
  {
    if (r < D)
    {
      Return_t t1=(D/r - 1.0);
      return std::exp(-t1*t1);
    }
    else
      return 1.0;
  }

  inline Return_t mycorrection(Return_t r1)
  {
    Return_t r2 = (r1*r1);
    Return_t rm2 = 1.0/r2;
    Return_t rm6 = std::pow(rm2,3);
    Return_t rm8 = rm6*rm2;
    Return_t rm10 = rm8*rm2;
    return (A*std::exp(alpha*r1) - (c1*rm6+c2*rm8+c3*rm10)*dampF(r1));
  }

  /** Do nothing */
  bool put(xmlNodePtr cur)
  {
    return true;
  }

  bool get(std::ostream& os) const
  {
    os << "HFDHE2Potential (T/S): " << PtclRef->getName();
    return true;
  }

  void add2Hamiltonian(ParticleSet& qp, TrialWaveFunction& psi , QMCHamiltonian& targetH);

  void addCorrection(QMCHamiltonian& targetH);

  //this is not necessary
  QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi)
  {
    return new HFDHE2Potential(qp);
  }

  void addObservables(PropertySetType& plist)
  {
    myIndex=plist.add("HFDHE2");
  }

  void setObservables(PropertySetType& plist)
  {
    plist[myIndex]=Value;
  }
};

/** Correction to HFDHE2Potential
 */
struct HFDHE2Potential_tail: public QMCHamiltonianBase
{

  const HFDHE2Potential* phyH;
  Return_t KValue;

  // epsilon = 3.42016039e-5, rm = 5.607384357
  // C6 = 1.460008056, C8 = 14.22016431, C10 = 187.2033646;
  HFDHE2Potential_tail(const HFDHE2Potential* org)
    :phyH(org)
  {
  }

  ~HFDHE2Potential_tail() { }

  void resetTargetParticleSet(ParticleSet& P)
  {
  }

  inline Return_t evaluate(ParticleSet& P)
  {
    Value += phyH->TCValue;
    KValue = phyH->kpre_factor*(Value+P.PropertyList[LOCALENERGY]);
    return Value;
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
    os << "HFDHE2PotentialTailcorr: " << PtclRef->getName();
    return true;
  }

  QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi)
  {
    return 0;
  }

  void registerObservables(std::vector<observable_helper*>& h5list, hid_t gid) const;

  void addObservables(PropertySetType& plist);

  void setObservables(PropertySetType& plist)
  {
    plist[myIndex]=Value;
    plist[myIndex+1]=KValue;
  }

  void setParticlePropertyList(PropertySetType& plist, int offset)
  {
    plist[myIndex+offset]=Value;
    plist[myIndex+1+offset]=KValue;
  }
};
}
#endif

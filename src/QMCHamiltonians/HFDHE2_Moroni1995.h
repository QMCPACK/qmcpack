//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: D.C. Yang, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: D.C. Yang, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#ifndef QMCPLUSPLUS_HFDHE2MORONI1995_H
#define QMCPLUSPLUS_HFDHE2MORONI1995_H
#include "Particle/ParticleSet.h"
#include "Particle/WalkerSetRef.h"
#include "Particle/DistanceTableData.h"
#include "Particle/DistanceTable.h"
#include "QMCHamiltonians/QMCHamiltonianBase.h"

namespace qmcplusplus
{
/** @ingroup hamiltonian
 *@brief HFDHE2_Moroni1995_phy for the indentical source and target particle sets.
 */
struct HFDHE2_Moroni1995_phy: public QMCHamiltonianBase
{
  Return_t rc,A,alpha,c1,c2,c3,D,mirrorshift,smooth;
  // remember that the default units are Hartree and Bohrs
  DistanceTableData* d_table;
  ParticleSet* PtclRef;
//    HFDHE2_Moroni1995_aux* dep;

  // epsilon = 3.42016039e-5, rm = 5.607384357
  // C6 = 1.460008056, C8 = 14.22016431, C10 = 187.2033646;
  HFDHE2_Moroni1995_phy(ParticleSet& P);

  ~HFDHE2_Moroni1995_phy() { }

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

  /** Do nothing */
  bool put(xmlNodePtr cur)
  {
    return true;
  }

  bool get(std::ostream& os) const
  {
    os << "HFDHE2_Moroni1995_phy (T/S): " << PtclRef->getName();
    return true;
  }

  void add2Hamiltonian(ParticleSet& qp, TrialWaveFunction& psi, QMCHamiltonian& targetH);

  void addCorrection(QMCHamiltonian& targetH);

  QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi)
  {
    return new HFDHE2_Moroni1995_phy(qp);
  }

  void addObservables(PropertySetType& plist)
  {
    myIndex=plist.add(myName);
    /*
    for (int i=0; i<plist.Names.size(); i++)
    std::cout<<plist.Names[i];
    std::cout<< "!!!!!" << std::endl;
    */
  }

  void setObservables(PropertySetType& plist)
  {
    plist[myIndex]=Value;
  }
};


struct HFDHE2_Moroni1995_aux: public QMCHamiltonianBase
{
  const HFDHE2_Moroni1995_phy* phyH;
  Return_t tailcorr;
//    Return_t KValue;

  HFDHE2_Moroni1995_aux(const HFDHE2_Moroni1995_phy* orig): phyH(orig)
  {
    Return_t rho = phyH->PtclRef->G.size()/phyH->PtclRef->Lattice.Volume, N0 = phyH->PtclRef->G.size(), rc = phyH->PtclRef->Lattice.WignerSeitzRadius;
    tailcorr = 2.0*M_PI*rho*N0*(-26.7433377905*std::pow(rc,-7.0) - 2.8440930339*std::pow(rc,-5.0)-0.486669351961 *std::pow(rc,-3.0)+ std::exp(-2.381392669*rc)*(2.75969257875+6.571911675726*rc+7.82515114293*rc*rc) );
  }

  ~HFDHE2_Moroni1995_aux() { }

  void resetTargetParticleSet(ParticleSet &P)
  {
    Return_t rho = phyH->PtclRef->G.size()/phyH->PtclRef->Lattice.Volume, N0 = phyH->PtclRef->G.size(), rc = phyH->PtclRef->Lattice.WignerSeitzRadius;
    tailcorr = 2.0*M_PI*rho*N0*(-26.7433377905*std::pow(rc,-7.0) - 2.8440930339*std::pow(rc,-5.0)-0.486669351961 *std::pow(rc,-3.0)+ std::exp(-2.381392669*rc)*(2.75969257875+6.571911675726*rc+7.82515114293*rc*rc) );
  }

  inline Return_t evaluate(ParticleSet &P)
  {
//      Value += phyH->smooth;
    Value = -(phyH->smooth) + tailcorr;
    // atomic units vs. kelvins
//      KValue = 0.0;
    return Value;
  }

  inline Return_t evaluate(ParticleSet &P, std::vector<NonLocalData>& Txy)
  {
    return evaluate(P);
  }

  bool put(xmlNodePtr cur)
  {
    return true;
  }

  bool get(std::ostream& os) const
  {
    os << "HFDHE2_Moroni1995_aux: " << phyH->PtclRef->getName();
    return true;
  }

  QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi)
  {
    return 0;
  }

  /*
  void addObservables(PropertySetType& plist);

  void setObservables(PropertySetType& plist);

  void setParticlePropertyList(PropertySetType& plist, int offset);
  */
};
}
#endif

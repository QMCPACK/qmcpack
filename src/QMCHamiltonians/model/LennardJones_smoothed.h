//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: D.C. Yang, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: D.C. Yang, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#ifndef QMCPLUSPLUS_LENNARDJONES_SMOOTHED_H
#define QMCPLUSPLUS_LENNARDJONES_SMOOTHED_H
#include "Particle/ParticleSet.h"
#include "Particle/WalkerSetRef.h"
#include "Particle/DistanceTableData.h"
#include "Particle/DistanceTable.h"
#include "QMCHamiltonians/QMCHamiltonianBase.h"

namespace qmcplusplus
{
/** @ingroup hamiltonian
 *@brief LennardJones_smoothed_phy for the indentical source and target particle sets.
 */
struct LennardJones_smoothed_phy: public QMCHamiltonianBase
{
  RealType epsilon, sigma, s6, rc, mirrorshift, smooth;
  // remember that the default units are Kelvins and Bohrs
  DistanceTableData* d_table;
  ParticleSet* PtclRef;

  LennardJones_smoothed_phy(ParticleSet& P, RealType e=3.23648511e-5, RealType s=4.830139998);

  ~LennardJones_smoothed_phy() { }

  void resetTargetParticleSet(ParticleSet& P);

  Return_t evaluate(ParticleSet& P);

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
    os << "LennardJones_smoothed_phy: " << PtclRef->getName();
    return true;
  }

  void add2Hamiltonian(ParticleSet& qp, TrialWaveFunction& psi, QMCHamiltonian& targetH);

  void addCorrection(QMCHamiltonian& targetH);

  QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi)
  {
    return new LennardJones_smoothed_phy(qp, epsilon, sigma);
  }

  void addObservables(PropertySetType& plist)
  {
    myIndex = plist.add(myName);
  }
};



/** @ingroup hamiltonian
 *@brief LennardJones_smoothed_np for the indentical source and target particle sets.
 */
struct LennardJones_smoothed_aux: public QMCHamiltonianBase
{
  const LennardJones_smoothed_phy* phyH;
  RealType epsilon, sigma, tailcorr, rc;

  // epsilon = 10.22 K = 3.236485111e-5 Ha
  LennardJones_smoothed_aux(const LennardJones_smoothed_phy* orig, RealType e=3.23648511e-5, RealType s=4.830139998): phyH(orig), epsilon(e), sigma(s)
  {
    rc = phyH->PtclRef->Lattice.WignerSeitzRadius;
    tailcorr = 8.0*M_PI*epsilon*std::pow(sigma,6.0)*(std::pow(sigma,6.0)-3.0*std::pow(rc,6.0))*std::pow(phyH->PtclRef->getTotalNum(),2.0)/(9.0*std::pow(rc,9.0)*phyH->PtclRef->Lattice.Volume);
    // Note the 2 powers of N
    app_log() << "  LennardJones_smoothed_aux tail correction is " << tailcorr << std::endl;
  }

  ~LennardJones_smoothed_aux() { }

  void resetTargetParticleSet(ParticleSet& P)
  {
    rc = P.Lattice.WignerSeitzRadius;
    tailcorr = 8.0*M_PI*epsilon*std::pow(sigma,6.0)*(std::pow(sigma,6.0)-3.0*std::pow(rc,6.0))*std::pow(P.getTotalNum(),2.0)/(9.0*std::pow(rc,9.0)*P.Lattice.Volume);
    app_log() << "  LennardJones_smoothed_aux tail correction is " << tailcorr << std::endl;
  }

  inline Return_t evaluate(ParticleSet& P)
  {
    Value = -(phyH->smooth) + tailcorr;
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
    os << "LennardJones_smoothed_aux (T/S): " << phyH->PtclRef->getName();
    return true;
  }

  QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi)
  {
    return 0;
  }
};
}
#endif

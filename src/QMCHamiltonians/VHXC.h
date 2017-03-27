//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#ifndef QMCPLUSPLUS_VHXC_H
#define QMCPLUSPLUS_VHXC_H
#include <QMCHamiltonians/QMCHamiltonianBase.h>
#include <LongRange/LRCoulombSingleton.h>

#if defined(HAVE_EINSPLINE)
#include <einspline/bspline.h>
#else
class UBspline_3d_d;
#endif
namespace qmcplusplus
{

/** @ingroup hamiltonian
 *\brief Calculates the DFT Hartreee, exchange and correlation potential
 */

class VHXC: public QMCHamiltonianBase
{
private:
  UBspline_3d_d *VSpline[2];
  void init_spline();
public:
  ParticleSet* PtclRef;

  bool FirstTime;
  int NumSpecies;
  int ChargeAttribIndx;
  int MemberAttribIndx;
  int NParticles;
  std::vector<RealType> Zat,Zspec;
  std::vector<int> NofSpecies;

  VHXC(ParticleSet& ref);

  /// copy constructor
  // VHXC(const VHXC& c);

  ~VHXC();

  void resetTargetParticleSet(ParticleSet& P);

  Return_t evaluate(ParticleSet& P);

  inline Return_t evaluate(ParticleSet& P, std::vector<NonLocalData>& Txy)
  {
    return evaluate(P);
  }

  /** Do nothing */
  bool put(xmlNodePtr cur);

  bool get(std::ostream& os) const
  {
    os << "VHXC potential: " << PtclRef->getName();
    return true;
  }

  QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi);

};

}
#endif



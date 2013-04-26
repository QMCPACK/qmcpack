//////////////////////////////////////////////////////////////////
// (c) Copyright 2009-  by Jeongnim Kim and Ken Esler
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
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
  vector<RealType> Zat,Zspec;
  vector<int> NofSpecies;

  VHXC(ParticleSet& ref);

  /// copy constructor
  // VHXC(const VHXC& c);

  ~VHXC();

  void resetTargetParticleSet(ParticleSet& P);

  Return_t evaluate(ParticleSet& P);

  inline Return_t evaluate(ParticleSet& P, vector<NonLocalData>& Txy)
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

/***************************************************************************
 * $RCSfile$   $Author: esler $
 * $Revision: 1581 $   $Date: 2007-01-04 10:02:14 -0600 (Thu, 04 Jan 2007) $
 * $Id: VHXC.h 1581 2007-01-04 16:02:14Z esler $
 ***************************************************************************/


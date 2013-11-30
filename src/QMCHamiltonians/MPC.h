//////////////////////////////////////////////////////////////////
// (c) Copyright 2008-  by Jeongnim Kim and Ken Esler
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
#ifndef QMCPLUSPLUS_MPC_H
#define QMCPLUSPLUS_MPC_H
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
 *\brief Calculates the Model Periodic Coulomb potential using PBCs
 */

class MPC: public QMCHamiltonianBase
{
protected:
  UBspline_3d_d *VlongSpline, *DensitySpline;
  double Vconst;
  void compute_g_G(double &g_0_N, vector<double> &g_G_N, int N);
  void init_gvecs();
  void init_f_G();
  void init_spline();
  double Ecut;
  vector<TinyVector<int,OHMMS_DIM> > Gints;
  vector<PosType> Gvecs;
  vector<ComplexType> Rho_G;
  TinyVector<int,OHMMS_DIM> SplineDim;
  int MaxDim;
  Return_t evalSR(ParticleSet& P) const;
  Return_t evalLR(ParticleSet& P) const;

public:
  ParticleSet* PtclRef;
  // Store the average electron charge density in reciprocal space
  vector<ComplexType> RhoAvg_G;
  vector<RealType> f_G;
  // The G=0 component
  double f_0;

  bool FirstTime;
  int NumSpecies;
  int ChargeAttribIndx;
  int MemberAttribIndx;
  int NParticles;
  RealType myConst;
  RealType myRcut;
  vector<RealType> Zat,Zspec;
  vector<int> NofSpecies;

  MPC(ParticleSet& ref, double cutoff);

  /// copy constructor
  // MPC(const MPC& c);

  ~MPC();

  void resetTargetParticleSet(ParticleSet& P);

  Return_t evaluate(ParticleSet& P);

  inline Return_t evaluate(ParticleSet& P, vector<NonLocalData>& Txy)
  {
    return evaluate(P);
  }

  /** implement all-walker stuff */
  virtual void addEnergy(MCWalkerConfiguration &W, vector<RealType> &LocalEnergy);

  /** Do nothing */
  bool put(xmlNodePtr cur);

  bool get(std::ostream& os) const
  {
    os << "MPC potential: " << PtclRef->getName();
    return true;
  }

  QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi);

  void initBreakup();
};

}
#endif

/***************************************************************************
 * $RCSfile$   $Author: esler $
 * $Revision: 1581 $   $Date: 2007-01-04 10:02:14 -0600 (Thu, 04 Jan 2007) $
 * $Id: MPC.h 1581 2007-01-04 16:02:14Z esler $
 ***************************************************************************/


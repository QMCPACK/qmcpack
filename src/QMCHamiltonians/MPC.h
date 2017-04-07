//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
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
  void compute_g_G(double &g_0_N, std::vector<double> &g_G_N, int N);
  void init_gvecs();
  void init_f_G();
  void init_spline();
  double Ecut;
  std::vector<TinyVector<int,OHMMS_DIM> > Gints;
  std::vector<PosType> Gvecs;
  std::vector<ComplexType> Rho_G;
  TinyVector<int,OHMMS_DIM> SplineDim;
  int MaxDim;
  Return_t evalSR(ParticleSet& P) const;
  Return_t evalLR(ParticleSet& P) const;

public:
  ParticleSet* PtclRef;
  // Store the average electron charge density in reciprocal space
  std::vector<ComplexType> RhoAvg_G;
  std::vector<RealType> f_G;
  // The G=0 component
  double f_0;

  bool FirstTime;
  int NumSpecies;
  int ChargeAttribIndx;
  int MemberAttribIndx;
  int NParticles;
  RealType myConst;
  RealType myRcut;
  std::vector<RealType> Zat,Zspec;
  std::vector<int> NofSpecies;

  MPC(ParticleSet& ref, double cutoff);

  /// copy constructor
  // MPC(const MPC& c);

  ~MPC();

  void resetTargetParticleSet(ParticleSet& P);

  Return_t evaluate(ParticleSet& P);

  inline Return_t evaluate(ParticleSet& P, std::vector<NonLocalData>& Txy)
  {
    return evaluate(P);
  }

  /** implement all-walker stuff */
  virtual void addEnergy(MCWalkerConfiguration &W, std::vector<RealType> &LocalEnergy);

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



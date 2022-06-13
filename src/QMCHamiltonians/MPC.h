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

#include "QMCHamiltonians/OperatorBase.h"
#include "LongRange/LRCoulombSingleton.h"

#if defined(HAVE_EINSPLINE)
#include "einspline/bspline.h"
#else
class UBspline_3d_d;
#endif
namespace qmcplusplus
{
/** @ingroup hamiltonian
 *\brief Calculates the Model Periodic Coulomb potential using PBCs
 */

class MPC : public OperatorBase
{
protected:
  std::shared_ptr<UBspline_3d_d> VlongSpline;
  //std::shared_ptr<UBspline_3d_d> DensitySpline;
  double Vconst;
  void compute_g_G(double& g_0_N, std::vector<double>& g_G_N, int N);
  void init_gvecs();
  void init_f_G();
  void init_spline();
  double Ecut;
  std::vector<TinyVector<int, OHMMS_DIM>> Gints;
  std::vector<PosType> Gvecs;
  std::vector<ComplexType> Rho_G;
  std::array<size_t, OHMMS_DIM> SplineDim;
  int MaxDim;
  Return_t evalSR(ParticleSet& P) const;
  Return_t evalLR(ParticleSet& P) const;
  // AA table ID
  const int d_aa_ID;

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
  std::vector<RealType> Zat, Zspec;
  std::vector<int> NofSpecies;

  MPC(ParticleSet& ref, double cutoff);

  /// copy constructor
  // MPC(const MPC& c);

  ~MPC() override;

  void resetTargetParticleSet(ParticleSet& P) override;

  Return_t evaluate(ParticleSet& P) override;

  /** implement all-walker stuff */
  void addEnergy(MCWalkerConfiguration& W, std::vector<RealType>& LocalEnergy) override;

  /** Do nothing */
  bool put(xmlNodePtr cur) override;

  bool get(std::ostream& os) const override
  {
    os << "MPC potential: " << PtclRef->getName();
    return true;
  }

  std::unique_ptr<OperatorBase> makeClone(ParticleSet& qp, TrialWaveFunction& psi) override;

  void initBreakup();
};

} // namespace qmcplusplus
#endif

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
    
    
#ifndef QMCPLUSPLUS_LOCAL_COREPOLPOTENTIAL_H
#define QMCPLUSPLUS_LOCAL_COREPOLPOTENTIAL_H
#include "Particle/ParticleSet.h"
#include "Particle/WalkerSetRef.h"
#include "QMCHamiltonians/QMCHamiltonianBase.h"
#include "OhmmsPETE/OhmmsMatrix.h"

namespace qmcplusplus
{

/** @ingroup hamiltonian
 @brief The effective core polarization potential.

   \f[
   V_{CPP} = -\frac{1}{2}\sum_{C}\alpha_C {\bf f}_C \cdot {\bf f}_C
   \f]
   Electric field which acts on core \f$C\f$ due to the charges of
   valence electrons \f$i\f$ and the other cores \f$C'\f$
   \f{eqnarray*}
   {\bf f}_C &=& \sum_i \frac{ {\bf r}_{Ci} }{r_{Ci}^3} C(r_{Ci},\rho_C)
   -\sum_{C' \ne C} \frac{ {\bf R}_{CC'} }{ R_{CC'}^3 }Z_{C'} =
   {\bf f}_C^e + {\bf f}_C^n \\
   {\bf r}_{Ci} &=& {\bf r}_i - {\bf r}_C\\
   {\bf R}_{CC'} &=& {\bf R}_{C'} - {\bf R}_{C}
   \f}
   \f$ C(r_{Ci},\rho_C) \f$ is a cut-off function for \f$ {\bf f}_C^e \f$
   with an adjustable parameter \f$ \rho_C \f$.

   \f{eqnarray*}
   V_{CPP} &=& -\frac{1}{2} \sum_C \left\{
   \;\;\; \sum_i \frac{1}{r_{Ci}^4} C^2(r_{Ci},\rho_C) +
   \sum_{i \ne j}\frac{ {\bf r}_{Ci} \cdot {\bf r}_{Ci} }{r^3_{Ci}r^3_{Ci}}
   C(r_{Ci},\rho_C) C^2(r_{Cj},\rho_C)
   \right\} \\
   & & -2 \left\{ \sum_i \sum_{C' \ne C} \frac{{\bf r}_{Ci} \cdot
   {\bf R}_{CC'}}{r^3_{Ci}R^3_{CC'}} Z_{C'}C(r_{Ci},\rho_C)
   + \left| \sum_{C' \ne C}  \frac{ {\bf R}_{CC'} }{ R^3_{CC'} } Z_{C'}
   \right|^2 \;\;\; \right\}
   \f}
*/
struct LocalCorePolPotential: public QMCHamiltonianBase
{

  /** core-polarization parameters for each species
   */
  struct CPP_Param
  {
    RealType alpha, C,
             r_b, one_over_rr;
    inline CPP_Param(RealType a=1.0, RealType r=1.0):
      alpha(a),C(-0.5*a),r_b(r),one_over_rr(1/r/r) {}

    inline void set(RealType a, RealType r)
    {
      alpha=a;
      C=-0.5*a;
      r_b=r;
      one_over_rr=1/r/r;
    }
    inline RealType operator()(RealType r)
    {
      RealType z=1.0-std::exp(-r*r*one_over_rr);
      return z*z;
    }
    bool put(xmlNodePtr cur);
  };

  ///boolean to evaluate the constant once
  bool FirstTime;
  ///the number of ions
  int nCenters;
  ///the number of electrons
  int nParticles;
  ///the CoreCore term
  RealType eCoreCore;

  ///reference to the ionic system
  ParticleSet& IonConfig;
  ///the ion-electron DistanceTable
  DistanceTableData* d_ie;
  ///the ion-ion DistanceTable
  DistanceTableData* d_ii;

  ///input CPP_Param whose size is the number of species
  std::vector<CPP_Param*> InpCPP;
  ///CPP_Param for each ion
  std::vector<CPP_Param*> Centers;

  ///CoreCoreDipole(C,C') \f$ =  Z_{C'} {\bf R}_{CC'}/ R_{CC'}^3 \f$
  std::vector<PosType> CoreCoreDipole;
  ///ElCoreDipole(C,i) \f$ = {\bf r}_{Ci} f(\bar{r}_{bCi}) /r_{Ci}^3\f$
  Matrix<PosType> CoreElDipole;

  ///constructor
  LocalCorePolPotential(ParticleSet& ions, ParticleSet& els);

  /////copy constructor
  //LocalCorePolPotential(const LocalCorePolPotential& cpp);

  ~LocalCorePolPotential();

  void resetTargetParticleSet(ParticleSet& P);

  Return_t evaluate(ParticleSet& P);

  inline Return_t evaluate(ParticleSet& P, std::vector<NonLocalData>& Txy)
  {
    return evaluate(P);
  }

  bool put(xmlNodePtr cur);

  bool get(std::ostream& os) const
  {
    os << "LocalCorePolPotential: " << IonConfig.getName();
    return true;
  }

  QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi);
  //inline RealType fcpp(RealType z) {
  //  return pow((1.0-exp(-1.0*z*z)),2);
  //}

};
}
#endif



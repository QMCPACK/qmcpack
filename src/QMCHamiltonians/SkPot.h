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
    
    
/** @file SkPot.h
 * @brief Declare SkPot
 */
#ifndef QMCPLUSPLUS_SK_POT_H
#define QMCPLUSPLUS_SK_POT_H
#include <QMCHamiltonians/QMCHamiltonianBase.h>
#include <LongRange/StructFact.h>
namespace qmcplusplus
{

/** SkPot evaluate the structure factor of the target particleset
 *
 * <estimator name="sk" type="sk" debug="no"/>
 */
class SkPot: public QMCHamiltonianBase
{
public:

  SkPot(ParticleSet& elns);

  void resetTargetParticleSet(ParticleSet& P);

  Return_t evaluate(ParticleSet& P);

  inline Return_t evaluate(ParticleSet& P, std::vector<NonLocalData>& Txy)
  {
    return evaluate(P);
  }
  bool put(xmlNodePtr cur);
  bool get(std::ostream& os) const;
  QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi);

  inline void FillFk()
  {
    for(int ki=0; ki<NumK; ki++)
    {
      RealType k=dot(sourcePtcl->SK->KLists.kpts_cart[ki],sourcePtcl->SK->KLists.kpts_cart[ki]);
      k= std::sqrt(k) - K_0;
      Fk[ki] = OneOverN*V_0*std::exp(-k*k);
//         app_log()<<ki<<": "<<Fk[ki] << std::endl;
    }
  }


protected:
  ParticleSet *sourcePtcl;
  /** number of species */
  int NumSpecies;
  /** number of kpoints */
  int NumK;
  /** number of kshells */
  int MaxKshell;
  /** normalization factor */
  RealType OneOverN;
  /** kshell counters */
  std::vector<int> Kshell;
  /** instantaneous structure factor  */
  std::vector<RealType> Kmag;
  /** 1.0/degenracy for a ksell */
  std::vector<RealType> OneOverDnk;
  /** \f$rho_k = \sum_{\alpha} \rho_k^{\alpha} \f$ for species index \f$\alpha\f$ */
  Vector<ComplexType> RhokTot;
  Vector<RealType> Fk;
  RealType V_0, K_0;
  /** resize the internal data
   *
   * The argument list is not completed
   */
  void resize();

  bool hdf5_out;
};

}
#endif


//////////////////////////////////////////////////////////////////
// (c) Copyright 2008-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
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

  inline Return_t evaluate(ParticleSet& P, vector<NonLocalData>& Txy)
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
//         app_log()<<ki<<": "<<Fk[ki] <<endl;
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
  vector<int> Kshell;
  /** instantaneous structure factor  */
  vector<RealType> Kmag;
  /** 1.0/degenracy for a ksell */
  vector<RealType> OneOverDnk;
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

/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 2945 $   $Date: 2008-08-05 10:21:33 -0500 (Tue, 05 Aug 2008) $
 * $Id: ForceBase.h 2945 2008-08-05 15:21:33Z jnkim $
 ***************************************************************************/

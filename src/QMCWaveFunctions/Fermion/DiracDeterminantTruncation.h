//////////////////////////////////////////////////////////////////
// (c) Copyright 1998-2002,2003- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
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
/**@file DiracDeterminantTruncationBase.h
 * @brief Declaration of DiracDeterminantTruncation with a S(ingle)P(article)O(rbital)SetBase
 */
#ifndef QMCPLUSPLUS_DIRACDETERMINANTTRUNCATION_H
#define QMCPLUSPLUS_DIRACDETERMINANTTRUNCATION_H
#include "QMCWaveFunctions/Fermion/DiracDeterminantBase.h"
#include "QMCWaveFunctions/OrbitalBase.h"
#include "QMCWaveFunctions/SPOSetBase.h"

namespace qmcplusplus
{

class DiracDeterminantTruncation: public DiracDeterminantBase
{
public:

  typedef SPOSetBase::IndexVector_t IndexVector_t;
  typedef SPOSetBase::ValueVector_t ValueVector_t;
  typedef SPOSetBase::ValueMatrix_t ValueMatrix_t;
  typedef SPOSetBase::GradVector_t  GradVector_t;
  typedef SPOSetBase::GradMatrix_t  GradMatrix_t;

  /** constructor
   *@param spos the single-particle orbital set
   *@param first index of the first particle
   */
  DiracDeterminantTruncation(SPOSetBasePtr const &spos, int first=0);

  ///default destructor
  ~DiracDeterminantTruncation();

  /**copy constructor
   * @param s existing DiracDeterminantTruncation
   *
   * This constructor makes a shallow copy of Phi.
   * Other data members are allocated properly.
   */
  DiracDeterminantTruncation(const DiracDeterminantTruncation& s);

  DiracDeterminantTruncation& operator=(const DiracDeterminantTruncation& s);
  const DistanceTableData* d_table;
  void set_truncation(int first, int nel,double &temp_cutoff,double &temp_radius);
  double radius;
  DiracDeterminantBase::ValueType ratio(ParticleSet& P, int iat);
  DiracDeterminantBase::ValueType ratio(ParticleSet& P, int iat,ParticleSet::ParticleGradient_t& dG, ParticleSet::ParticleLaplacian_t& dL);

  void resize(int nel, int morb);

  //   void set_iterative(int first,int nel, double &temp_cutoff);

  ValueMatrix_t psiM_actual;
  ValueMatrix_t psiM2;
  ValueMatrix_t temp_psiM2;
  ValueVector_t psi_diff;
  void UpdatePsiM2(ValueVector_t &vec,int ptcl);

  void ChooseNearbyParticles(int ptcl,list<int> &nearbyPtcls);


  /** move was accepted, update the real container
   */
  void acceptMove(ParticleSet& P, int iat);



  DiracDeterminantBase::RealType
  evaluateLog(ParticleSet& P,
              ParticleSet::ParticleGradient_t& G,
              ParticleSet::ParticleLaplacian_t& L);

  DiracDeterminantBase::RealType
  evaluateLog(ParticleSet& P, PooledData<RealType>& buf) ;


  vector<list<pair<int,double> > > particleLists;
  list<pair<int,double> >  oldPtcl;
  double cutoff;

};
}

#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 3265 $   $Date: 2008-10-15 09:20:33 -0500 (Wed, 15 Oct 2008) $
 * $Id: DiracDeterminantTruncation.h 3265 2008-10-15 14:20:33Z jnkim $
 ***************************************************************************/

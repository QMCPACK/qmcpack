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
/**@file DiracDeterminantIterativeBase.h
 * @brief Declaration of DiracDeterminantIterative with a S(ingle)P(article)O(rbital)SetBase
 */
#ifndef QMCPLUSPLUS_DIRACDETERMINANTITERATIVE_H
#define QMCPLUSPLUS_DIRACDETERMINANTITERATIVE_H
#include "QMCWaveFunctions/Fermion/DiracDeterminantBase.h"
#include "QMCWaveFunctions/OrbitalBase.h"
#include "QMCWaveFunctions/SPOSetBase.h"

namespace qmcplusplus
{

class DiracDeterminantIterative: public DiracDeterminantBase
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
  DiracDeterminantIterative(SPOSetBasePtr const &spos, int first=0);

  ///default destructor
  ~DiracDeterminantIterative();

  /**copy constructor
   * @param s existing DiracDeterminantIterative
   *
   * This constructor makes a shallow copy of Phi.
   * Other data members are allocated properly.
   */
  DiracDeterminantIterative(const DiracDeterminantIterative& s);

  DiracDeterminantIterative& operator=(const DiracDeterminantIterative& s);


  DiracDeterminantBase::ValueType ratio(ParticleSet& P, int iat);
  DiracDeterminantBase::ValueType ratio(ParticleSet& P, int iat,ParticleSet::ParticleGradient_t& dG, ParticleSet::ParticleLaplacian_t& dL);

  void resize(int nel, int morb);
  void set(int first, int nel);
  void set_iterative(int first,int nel, double &temp_cutoff);

  void SparseToCSR(vector<int> &Arp, vector<int> &Ari,vector<double> &Arx);



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
 * $Id: DiracDeterminantIterative.h 3265 2008-10-15 14:20:33Z jnkim $
 ***************************************************************************/

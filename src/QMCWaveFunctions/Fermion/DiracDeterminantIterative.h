//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Bryan Clark, bclark@Princeton.edu, Princeton University
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Bryan Clark, bclark@Princeton.edu, Princeton University
//////////////////////////////////////////////////////////////////////////////////////
    
    
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

  void resize(int nel, int morb);
  void set(int first, int nel);
  void set_iterative(int first,int nel, double &temp_cutoff);

  void SparseToCSR(std::vector<int> &Arp, std::vector<int> &Ari,std::vector<double> &Arx);



  DiracDeterminantBase::RealType
  evaluateLog(ParticleSet& P,
              ParticleSet::ParticleGradient_t& G,
              ParticleSet::ParticleLaplacian_t& L);

  std::vector<std::list<std::pair<int,double> > > particleLists;
  std::list<std::pair<int,double> >  oldPtcl;
  double cutoff;

};
}

#endif

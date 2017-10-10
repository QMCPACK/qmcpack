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
    
    
#ifndef QMCPLUSPLUS_SINGLEPARTICLEORBITALSET_H
#define QMCPLUSPLUS_SINGLEPARTICLEORBITALSET_H
#include "QMCWaveFunctions/SPOSetBase.h"

namespace qmcplusplus
{

/**a set of single-particle orbitals.
 *
 * This class provides necessary interfaces for SlaterDeterminant<SPOSet> and SlaterDet
 * and can be used by any orbital representation that has an one-to-one
 * mapping between the function evaluations and the column indeices of a
 * Dirac determinant.
 * The template parameter is any function with a member function
 * which can provde the value, gradient and laplacian at a point.
 * value_type OT::evaluate(const point_type&, gradient_type&, value_type& )
 * Example classes can be found in Numerics/CosineFunction.h
 */
template<class OT>
struct SingleParticleOrbitalSet: public SPOSetBase
{

  ///the type of single-particle orbtials
  typedef OT                      SPOrbital_t;
  typedef std::vector<OT*>             SPOContainer_t;
  typedef typename OT::value_type value_type;

  SPOContainer_t Phi;

  ///constructor
  SingleParticleOrbitalSet(int norbs=0)
  {
    setOrbitalSetSize(norbs);
  }

  /**add a single-particle orbital */
  int add(SPOrbital_t* afunction )
  {
    Phi.push_back(afunction);
    return Phi.size()-1;
  }

  void setOrbitalSetSize(int norbs)
  {
    if(norbs == OrbitalSetSize )
      return;
    OrbitalSetSize=norbs;
    BasisSetSize=norbs;
  }

  void resetParameters(VarRegistry<RealType>& optVariables)
  {
    for(int i=0; i<Phi.size(); i++)
      Phi[i]->resetParameters(optVariables);
  }

  void resetTargetParticleSet(ParticleSet& P) { }

  void
  evaluate(const ParticleSet& P, int iat, ValueVector_t& psi)
  {
    for(int j=0; j<OrbitalSetSize; j++)
      psi[j]=Phi[j]->evaluate(P.R[iat]);
    //vector<SPOrbital_t*>::iterator it(Phi.begin()),it_end(Phi.end());
    //int j(0);
    //while(it != it_end) {
    //  psi[j++]=(*it)->evaluate(P.R[iat]);++it;
    //}
  }

  void
  evaluate(const ParticleSet& P, int iat,
           ValueVector_t& psi, GradVector_t& dpsi, ValueVector_t& d2psi)
  {
    for(int j=0; j<OrbitalSetSize; j++)
    {
      psi[j]=Phi[j]->evaluate(P.R[iat],dpsi[j],d2psi[j]);
    }
    //vector<SPOrbital_t*>::iterator it(Phi.begin()),it_end(Phi.end());
    //int j(0);
    //while(it != it_end) {
    //  psi[j]=(*it)->evaluate(P.R[iat],dpsi[j],d2psi[j]);++it;j++;
    //}
  }

  void evaluate(const ParticleSet& P, int first, int last,
                ValueMatrix_t& logdet, GradMatrix_t& dlogdet, ValueMatrix_t& d2logdet)
  {
    for(int j=0; j<OrbitalSetSize; j++)
    {
      for(int i=0,iat=first; iat<last; i++,iat++)
      {
        logdet(j,i)=Phi[j]->evaluate(P.R[iat],dlogdet(i,j),d2logdet(i,j));
      }
    }
    //int n = last-first;
    //int iat = first;
    //vector<SPOrbital_t*>::iterator it_end(Phi.end());
    //for(int i=0; i<n; i++,iat++) {
    //  std::vector<SPOrbital_t*>::iterator it(Phi.begin());
    //  int j(0);
    //  while(it != it_end) {
    //    logdet(j,i)= (*it)->evaluate(P.R[iat], dlogdet(i,j),d2logdet(i,j));
    //    ++it;++j;
    //  }
    //}
  }
};
}
#endif

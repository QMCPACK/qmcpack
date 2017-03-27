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
    
    
#ifndef QMCPLUSPLUS_MIXED_SINGLEPARTICLEBASISSET_H
#define QMCPLUSPLUS_MIXED_SINGLEPARTICLEBASISSET_H

#include "Configuration.h"
#include "Numerics/TriCubicSplineT.h"
#include "QMCWaveFunctions/DummyBasisSet.h"

namespace qmcplusplus
{

/** class to handle mixed single-particle orbitals
 *
 * For now, use MolecularOrbitals using numerical grid. It is possible to
 * generalize it.
 */
template<class LOType>
class MixedSPOSet: public QMCTraits
{

public:

  //@typedef numerical orbitals
  typedef TriCubicSplineT<ValueType,RealType>  NGOType;
  typedef typename LOType::BasisSet_t          BasisSet_t;

  ///the number of Orbitals
  int NumberOfOrbitals;

  ///pointer to the localized orbitals
  LOType* LocalizedOrbitals;
  ///a set of numerical orbitals
  std::vector<NGOType*> GridOrbitals;

  /** constructor
   */
  MixedSPOSet():NumberOfOrbitals(0), LocalizedOrbitals(0) { }

  /** destructor
   */
  ~MixedSPOSet() { }

  ///return the number of single particle orbitals
  inline int size() const
  {
    return NumberOfOrbitals;
  }

  /** Add a numerical orbital
   * @param ngorb a pointer to a numerical orbital
   */
  void add(NGOType* ngorb)
  {
    GridOrbitals.push_back(ngorb);
    NumberOfOrbitals=GridOrbitals.size();
  }

  void setLocalizedOrbitals(LOType* lo)
  {
    LocalizedOrbitals=lo;
  }
  /** reset the component
   *
   * GridOrbitals are not reset (it is hard to optimize it).
   */
  inline void reset()
  {
    if(LocalizedOrbitals)
      LocalizedOrbitals->reset();
    NumberOfOrbitals=GridOrbitals.size();
    //for(int j=0; j<NumberOfOrbitals; j++) GridOrbitals[j]->reset();
  }

  void resetTargetParticleSet(ParticleSet& P)
  {
    if(LocalizedOrbitals)
      LocalizedOrbitals->resetTargetParticleSet(P);
  }


  /** evaluate the values of the single-particle orbitals
   *
   * Add the value from LocalizedOrbitals and NumericalOrbtials.
   */
  inline ValueType
  evaluate(const ParticleSet& P, int iat, int jorb)
  {
    ValueType res=GridOrbitals[jorb]->evaluate(P.R[iat]);
    if(LocalizedOrbitals)
      res+=LocalizedOrbitals->evaluate(P,iat,jorb);
    return res;
  }

  /** evaluate the values of the single-particle orbitals
   *@param P input configuration containing N particles
   *@param iat particle index
   *@param psi values of the single-particle orbitals
   *
   * This function is introduced to function evaluations.
   */
  template<class VV>
  inline void
  evaluate(const ParticleSet& P, int iat, VV& psi)
  {
    if(LocalizedOrbitals)
    {
      LocalizedOrbitals->evaluate(P,iat,psi);
      for(int j=0; j<NumberOfOrbitals; j++)
      {
        psi[j]+=GridOrbitals[j]->evaluate(P.R[iat]);
      }
    }
    else
    {
      for(int j=0; j<NumberOfOrbitals; j++)
      {
        psi[j]=GridOrbitals[j]->evaluate(P.R[iat]);
      }
    }
  }

  /**@ingroup particlebyparticle
   *@brief evaluate the values of the single-particle orbitals for the iat-th particle
   *@param P input configuration containing N particles
   *@param iat particle index
   *@param psi values of the single-particle orbitals
   *@param dpsi gradients of the single-particle orbitals
   *@param d2psi laplacians of the single-particle orbitals
   *
   * This function completes a row of a Dirac Deterimant to perform ratio/update.
   * The particle index identifies the particle whose position is updated.
   */
  template<class VV, class GV>
  inline void
  evaluate(const ParticleSet& P, int iat, VV& psi, GV& dpsi, VV& d2psi)
  {
    if(LocalizedOrbitals)
    {
      LocalizedOrbitals->evaluate(P,iat,psi,dpsi,d2psi);
      typename GV::value_type grad;
      typename VV::value_type lap;
      for(int j=0; j<NumberOfOrbitals; j++)
      {
        psi[j]+=GridOrbitals[j]->evaluate(P.R[iat],grad,lap);
        dpsi[j]+=grad;
        d2psi[j]+=lap;
      }
    }
    else
    {
      for(int j=0; j<NumberOfOrbitals; j++)
      {
        psi[j]=GridOrbitals[j]->evaluate(P.R[iat],dpsi[j],d2psi[j]);
      }
    }
  }

  /** complete the values of the single-particle orbitals and their gradients and laplacians
   *@param P input configuration containing N particles
   *@param first index of the first particle
   *@param last index of the last particle
   *@param logdet \f$ logdet(j,i) = \sum_k C_{jk} \phi_{k}({\bf r}_i) \f$
   *@param dlogdet \f$ dlogdet(i,j) =\sum_k C_{jk} \nabla_i \phi_{k}({\bf r}_i) \f$
   *@param d2logdet \f$ d2logdet(i,j) = \sum_k C_{jk} \nabla^2_i \phi_{k}({\bf r}_i) \f$
   *
   *Each object handles [first,last) quantum particles.
   *The physical particle index maps to the orbital index of [0,last-first).
   */
  template<class VM, class GM>
  inline void
  evaluate(const ParticleSet& P, int first, int last,
           VM& logdet, GM& dlogdet, VM& d2logdet)
  {
    if(LocalizedOrbitals)
    {
      LocalizedOrbitals->evaluate(P,first,last,logdet,dlogdet,d2logdet);
      typename GM::value_type grad;
      typename VM::value_type lap;
      for(int i=0,iat=first; iat<last; iat++, i++)
      {
        for(int j=0; j<NumberOfOrbitals; j++)
        {
          logdet(j,i)+=GridOrbitals[j]->evaluate(P.R[iat],grad,lap);
          dlogdet(i,j)+=grad;
          d2logdet(i,j)+=lap;
        }
      }
    }
    else
    {
      for(int i=0,iat=first; iat<last; iat++, i++)
      {
        for(int j=0; j<NumberOfOrbitals; j++)
        {
          logdet(j,i)=GridOrbitals[j]->evaluate(P.R[iat],dlogdet(i,j),d2logdet(i,j));
        }
      }
    }
  }

  template<class VM, class GM>
  inline void
  evaluate(const WalkerSetRef& W, int first, int last,
           std::vector<VM>& logdet, std::vector<GM>& dlogdet, std::vector<VM>& d2logdet)
  {
    //this is useless
  }
};
}
#endif



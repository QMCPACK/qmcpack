//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



#ifndef QMCPLUSPLUS_RMCLOCALENERGYESTIMATOR_H
#define QMCPLUSPLUS_RMCLOCALENERGYESTIMATOR_H
#include "Estimators/ScalarEstimatorBase.h"
#include "QMCHamiltonians/QMCHamiltonian.h"
#include "Particle/Reptile.h"

namespace qmcplusplus
{

/** Class to accumulate the local energy and components
 *
 * Use Walker::Properties to accumulate Hamiltonian-related quantities.
 */
class RMCLocalEnergyEstimator: public ScalarEstimatorBase
{


  int FirstHamiltonian;
  int SizeOfHamiltonians;
  int NObs;
  int RMCSpecificTerms;
  const QMCHamiltonian& refH;

public:

  /** constructor
   * @param h QMCHamiltonian to define the components
   */
  RMCLocalEnergyEstimator(QMCHamiltonian& h, int nobs=2);

  /** accumulation per walker
   * @param awalker current walker
   * @param wgt weight
   *
   * Weight of observables should take into account the walkers weight. For Pure DMC. In branching DMC set weights to 1.
   */
  inline void accumulate(const Walker_t& awalker, RealType wgt)
  {
  }

  /*@{*/
  inline void accumulate(const MCWalkerConfiguration& W
                         , WalkerIterator first, WalkerIterator last, RealType wgt)
  {
    //WalkerIterator tail=first+W.activeBead+W.direction;
    //WalkerIterator head=first+W.activeBead;
    //if(tail>=last)
    //  tail-=last-first;
    //else if(tail<first)
    //  tail+=last-first;
    Walker_t& head = W.reptile->getHead();
    Walker_t& tail = W.reptile->getTail();
    Walker_t& center = W.reptile->getCenter();
//       mixed estimator stuff
    const RealType* restrict ePtr = head.getPropertyBase();
    const RealType* restrict lPtr = tail.getPropertyBase();
    const RealType* restrict cPtr = center.getPropertyBase();
    //   RealType wwght=  head.Weight;
    RealType wwght=0.5;
    //app_log()<<"~~~~~For head:  Energy:"<<ePtr[LOCALENERGY]<< std::endl;
    scalars[0](0.5*(ePtr[LOCALENERGY]+lPtr[LOCALENERGY]),wwght);
    scalars[1](0.5*(ePtr[LOCALENERGY]*ePtr[LOCALENERGY]+lPtr[LOCALENERGY]*lPtr[LOCALENERGY]),wwght);
    scalars[2](cPtr[LOCALENERGY],wwght);
    scalars[3](cPtr[LOCALENERGY]*cPtr[LOCALENERGY],wwght);
    scalars[4](ePtr[LOCALENERGY]*lPtr[LOCALENERGY],wwght);
    scalars[5](0.5*(ePtr[LOCALPOTENTIAL]+lPtr[LOCALPOTENTIAL]),wwght);
    scalars[6](cPtr[LOCALPOTENTIAL],wwght);

    for(int target=RMCSpecificTerms, source=FirstHamiltonian; source<FirstHamiltonian+SizeOfHamiltonians; ++target, ++source)
    {
      wwght=0.5;
      scalars[target](lPtr[source],wwght);
      scalars[target](ePtr[source],wwght);
    }
    for(int target=RMCSpecificTerms+SizeOfHamiltonians, source=FirstHamiltonian; source<FirstHamiltonian+SizeOfHamiltonians; ++target, ++source)
    {
      wwght=1;
      scalars[target](cPtr[source],wwght);
    }
    //     scalars[target](lPtr[source],wwght);
    //    int stride(0);
    //    int bds(last-first);
    //    if(bds%2>0)//odd
    //    {
    //      int odd=(bds+1)/2;
    //      stride=odd/NObs;
    //    }
    //   else
    //   {
    //     int odd=bds/2;
    //     stride=odd/NObs;
    //   }
    //
    //   int indx(4+SizeOfHamiltonians);
    //   for(int j=0; j < NObs; j++)
    //  {
    //    tail+= -W.direction*stride;
    //    head+=  W.direction*stride;
    //   if(tail>=last)
    //     tail-=last-first;
    //   else if(tail<first)
    //     tail+=last-first;
    //   if(head>=last)
    //     head-=last-first;
    //   else if(head<first)
    //    head+=last-first;
    //     const RealType* lPtr = tail.getPropertyBase();
    // const RealType* ePtr = head.getPropertyBase();
    // for(int target=0, source=FirstHamiltonian; target<SizeOfHamiltonians; ++target, ++source, indx++)
    //   {
    //    scalars[indx](lPtr[source]);
    // scalars[indx](ePtr[source]);
    //   }
//     }
    //   int maxage(0);
//     while(first!=last)
    //   {
    //     maxage= std::max(maxage,(*first)->Age);
    //      first++;
    //   }
    //  scalars[3](maxage);
//       for(; first != last; ++first) std::accumulate(**first,wgt);*/
  }

  void add2Record(RecordListType& record);
  void registerObservables(std::vector<observable_helper*>& h5dec, hid_t gid) {}
  ScalarEstimatorBase* clone();
  /*@}*/
};
}
#endif

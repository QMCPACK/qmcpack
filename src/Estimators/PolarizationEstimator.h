//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: D. Das, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: D. Das, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



#ifndef QMCPLUSPLUS_POLARIZATIONESTIMATOR_H
#define QMCPLUSPLUS_POLARIZATIONESTIMATOR_H
#include <fstream>
#include "OhmmsData/libxmldefs.h"
#include "Estimators/ScalarEstimatorBase.h"

namespace qmcplusplus
{

template<class T>
class PolarizationEstimator: public ScalarEstimatorBase<T>
{

  using ScalarEstimatorBase<T>::CollectSum;
  using ScalarEstimatorBase<T>::b_average;
  using ScalarEstimatorBase<T>::b_variance;
  ///local data
  //T z_sum, z_sq_sum;
  std::vector<T> z_sum;
  int pindex_0, pindex_1;
public:

  typedef typename ScalarEstimatorBase<T>::Walker_t Walker_t;
  typedef typename ScalarEstimatorBase<T>::WalkerIterator WalkerIterator;

  PolarizationEstimator()
  {
    z_sum.resize(2,0.0);
  }

  void add2Record(RecordNamedProperty<T>& record, BufferType& msg)
  {
    pindex_0 = record.add("Pol-z");
    pindex_1 = record.add("Pol-z-var");
  }

  inline void accumulate(const Walker_t& awalker, T wgt, BufferType& msg)
  {
    for(int i=0; i<awalker.size(); i++)
    {
      T z = awalker.R[i][2];
      z_sum[0] += wgt*z;
      z_sum[1] += wgt*z*z;
    }
  }

  void accumulate(const MCWalkerConfiguration& W
                  , WalkerIterator first, WalkerIterator last)
  {
    for(; first != last; ++first)
      std::accumulate(**first,(*first)->Weight);
  }

  ///reset all the cumulative sums to zero
  inline void reset()
  {
    z_sum[0]=T();
    z_sum[1]=T();
  }

  inline void report(RecordNamedProperty<T>& record, T wgtinv)
  {
#ifdef HAVE_MPI
    if(CollectSum)
      gsum(z_sum,0);
#endif
    b_average =  z_sum[0]*wgtinv;
    b_variance = z_sum[1]*wgtinv-b_average*b_average;
    record[pindex_0]=b_average;
    record[pindex_1]=b_variance;
    reset();
  }
};
}
#endif

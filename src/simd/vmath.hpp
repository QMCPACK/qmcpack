//////////////////////////////////////////////////////////////////
// (c) Copyright 2011- by Jeongnim Kim
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
#ifndef QMCPLUSPLUS_VECTORIZED_STDMATH_HPP
#define QMCPLUSPLUS_VECTORIZED_STDMATH_HPP

#include <cmath>

namespace qmcplusplus {

  namespace simd 
  {
    /**  mod on an array
     * out[i]=in[i]-floor(in[i])
     */
    template<typename T, typename SIZET>
      inline void remainder(const T* restrict in, T* restrict out, SIZET n)
      {
        for(SIZET i=0; i<n; ++i) out[i]=in[i]-std::floor(in[i]);
      }

    template<typename T, typename SIZET>
      inline void remainder(T* restrict inout, SIZET n)
      {
        for(SIZET i=0; i<n; ++i) inout[i]-=std::floor(inout[i]);
      }

  }
}
#endif
/***************************************************************************
 * $RCSfile$   $Author: jmcminis $
 * $Revision: 5077 $   $Date: 2010-12-09 03:14:51 -0600 (Thu, 09 Dec 2010) $
 * $Id: vmath.hpp 5077 2010-12-09 09:14:51Z jmcminis $ 
 ***************************************************************************/

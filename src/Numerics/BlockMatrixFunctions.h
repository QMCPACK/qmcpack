//////////////////////////////////////////////////////////////////
// (c) Copyright 2005- by Jeongnim Kim
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
#ifndef QMCPLUSPLUS_BLOCKMATRIXFUNCTIONS_H
#define QMCPLUSPLUS_BLOCKMATRIXFUNCTIONS_H
namespace qmcplusplus
{
template<typename T, unsigned D>
struct GEMV
{
  static inline
  void apply(const T* restrict a, const T* restrict x, T* restrict y, int n, int m)
  {
    for(int i=0; i<n; i++)
    {
      T tmp = 0;
      const T* restrict xcopy=x;
      for(int j=0; j<m; j++ )
      {
        tmp += (*a++)*(*xcopy++);
      }
      (*y++) = tmp;
    }
  }
};

/** specialization for D=3
 */
template<typename T>
struct GEMV<T,3>
{
  static inline
  void apply(const T* restrict a, const T* restrict x, T* restrict y, int n, int m)
  {
    y[0]=a[0]*x[0]+a[1]*x[1]+a[2]*x[2];
    y[1]=a[3]*x[0]+a[4]*x[1]+a[5]*x[2];
    y[2]=a[6]*x[0]+a[7]*x[1]+a[8]*x[2];
  }
};
}
#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1229 $   $Date: 2006-08-07 09:15:25 -0500 (Mon, 07 Aug 2006) $
 * $Id: MatrixOperators.h 1229 2006-08-07 14:15:25Z jnkim $
 ***************************************************************************/

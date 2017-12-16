//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



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

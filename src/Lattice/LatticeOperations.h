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
    
    


#ifndef OHMMS_LATTICEOPERATIONS_H
#define OHMMS_LATTICEOPERATIONS_H
#include "OhmmsPETE/TinyVector.h"
#include "OhmmsPETE/Tensor.h"

namespace qmcplusplus
{

/** Dot product between a vector and tensor
 *
 * This is a dummy template class to be specialized.
 * Each specialization implements a static function apply to return
 * the dot product of vector and tensor.
 */
template<class T1, class T2, bool ORTHO> struct DotProduct { };

/** Specialization of DotProduct for a general tensor in D-dimension.
 */
template<class T1, class T2, unsigned D>
struct DotProduct<TinyVector<T1,D>,Tensor<T2,D>,false>
{
  typedef typename BinaryReturn<T1,T2,OpMultiply>::Type_t Type_t;
  inline static TinyVector<Type_t,D>
  apply(const TinyVector<T1,D>& lhs, const Tensor<T2,D>& rhs)
  {
    return OTDot< TinyVector<T1,D> , Tensor<T2,D> > :: apply(lhs,rhs);
  }
};


/** Specialization of DotProduct for a general tensor in 3-dimension.
 */
template<class T1, class T2>
struct DotProduct<TinyVector<T1,3>,Tensor<T2,3>,false>
{
  typedef typename BinaryReturn<T1,T2,OpMultiply>::Type_t Type_t;

  inline static TinyVector<Type_t,3>
  apply(const TinyVector<T1,3>& lhs, const Tensor<T2,3>& rhs)
  {
    return TinyVector<Type_t,3>(lhs[0]*rhs[0]+lhs[1]*rhs[3]+lhs[2]*rhs[6],
                                lhs[0]*rhs[1]+lhs[1]*rhs[4]+lhs[2]*rhs[7],
                                lhs[0]*rhs[2]+lhs[1]*rhs[5]+lhs[2]*rhs[8]);
  }

  inline static TinyVector<Type_t,3>
  apply( const Tensor<T2,3>& lhs, const TinyVector<T1,3>& rhs)
  {
    return TinyVector<Type_t,3>(lhs[0]*rhs[0]+lhs[1]*rhs[1]+lhs[2]*rhs[2],
                                lhs[3]*rhs[0]+lhs[4]*rhs[1]+lhs[5]*rhs[2],
                                lhs[6]*rhs[0]+lhs[7]*rhs[1]+lhs[8]*rhs[2]);
  }
};

/** Specialization of DotProduct for a diagonal (ORTHO=true) in D-dimension.
 */
template<class T1, class T2, unsigned D>
struct DotProduct<TinyVector<T1,D>,Tensor<T2,D>,true>
{
  typedef typename BinaryReturn<T1,T2,OpMultiply>::Type_t Type_t;
  inline static TinyVector<Type_t,D>
  apply(const TinyVector<T1,D>& lhs, const Tensor<T2,D>& rhs)
  {
    TinyVector<Type_t,D> ret;
    for(unsigned d=0; d<D; ++d)
      ret[d]=lhs[d]*rhs(d,d);
    return ret;
  }
};

/** Specialization of DotProduct for a diagonal (ORTHO=true) in 3-dimension.
 */
template<class T1, class T2>
struct DotProduct<TinyVector<T1,3>,Tensor<T2,3>,true>
{
  typedef typename BinaryReturn<T1,T2,OpMultiply>::Type_t Type_t;
  inline static TinyVector<Type_t,3>
  apply(const TinyVector<T1,3>& lhs, const Tensor<T2,3>& rhs)
  {
    return TinyVector<Type_t,3>(lhs[0]*rhs[0],lhs[1]*rhs[4],lhs[2]*rhs[8]);
  }

  inline static TinyVector<Type_t,3>
  apply(const Tensor<T2,3>& lhs, const TinyVector<T1,3>& rhs)
  {
    return TinyVector<Type_t,3>(lhs[0]*rhs[0],lhs[4]*rhs[1],lhs[8]*rhs[2]);
  }
};

///////////////////////////////////////////////////////////////////
// Cartesian norm2: ||a-b|| for vectos a and b in the relative coordinate
///////////////////////////////////////////////////////////////////
/** Dummy template class to be specialized **/
template<class T1, class T2, bool ORTHO> struct CartesianNorm2 {};

/** Specialization for general metric tensor for the D-dimensional system
 */
template<class T1, class T2, unsigned D>
struct CartesianNorm2<TinyVector<T1,D>,Tensor<T2,D>,false>
{
  typedef typename BinaryReturn<T1,T2,OpMultiply>::Type_t Type_t;

  /** Function to evaluate a Cartesian norm2
   * @param a a relative vector
   * @param X metric tensor
   * @param b relative vector
   * @return dot(a,dot(X,b))
   */
  inline static Type_t
  apply(const TinyVector<T1,D>& a, const Tensor<T2,D>& X, const TinyVector<T1,D>& b)
  {
    return dot(a,dot(X,b));
  }
};

/** Specialization for general metric tensor for the 3-dimensional system
 */
template<class T1, class T2>
struct CartesianNorm2<TinyVector<T1,3>,Tensor<T2,3>,false>
{
  typedef typename BinaryReturn<T1,T2,OpMultiply>::Type_t Type_t;

  /** Function to evaluate a Cartesian norm2
   * @param a a relative vector
   * @param X metric tensor
   * @param b relative vector
   * @return dot(a,dot(X,b))
   */
  inline static Type_t
  apply(const TinyVector<T1,3>& a, const Tensor<T2,3>& X, const TinyVector<T1,3>& b)
  {
    return a[0]*(X[0]*b[0]+X[1]*b[1]+X[2]*b[2])
           + a[1]*(X[3]*b[0]+X[4]*b[1]+X[5]*b[2])
           + a[2]*(X[6]*b[0]+X[7]*b[1]+X[8]*b[2]);
  }
};


/** Specialization for orthogonal metric tensor for the D-dimensional system
 */
template<class T1, class T2, unsigned D>
struct CartesianNorm2<TinyVector<T1,D>,Tensor<T2,D>,true>
{
  typedef typename BinaryReturn<T1,T2,OpMultiply>::Type_t Type_t;
  inline static Type_t
  apply(const TinyVector<T1,D>& a, const Tensor<T2,D>& X, const TinyVector<T1,D>& b)
  {
    Type_t res = a[0]*X(0,0)*b[0];
    for(int d=1; d<D; d++)
      res += a[d]*X(d,d)*b[d];
    return res;
  }
};

/** Specialization for orthogonal metric tensor for the 3-dimensional system
 */
template<class T1, class T2>
struct CartesianNorm2<TinyVector<T1,3>,Tensor<T2,3>,true>
{
  typedef typename BinaryReturn<T1,T2,OpMultiply>::Type_t Type_t;
  inline static Type_t
  apply(const TinyVector<T1,3>& a, const Tensor<T2,3>& X, const TinyVector<T1,3>& b)
  {
    return a[0]*X[0]*b[0]+a[1]*X[4]*b[1]+a[2]*X[8]*b[2];
  }
};


template<class T, unsigned D>
struct MinimumImageBConds
{
  typedef TinyVector<T,D> Return_t;
  inline static void apply(const Tensor<T,D>& R, const Tensor<T,D>&  G, TinyVector<T,D>& r)
  {
    Return_t u=dot(r,G);
    for(int i=0; i<D; ++i)
      u[i]=u[i]-round(u[i]);
    r=dot(u,R);
  }

};

template<class T, unsigned D>
struct CheckBoxConds
{
  inline static bool inside(const TinyVector<T,D>& u)
  {
    bool yes=(u[0]>0.0 && u[0]<1.0);
    for(int i=1; i<D; ++i)
      yes &=(u[i]>0.0 && u[i]<1.0);
    return yes;
  }

  inline static bool inside(const TinyVector<T,D>& u, TinyVector<T,D>& ubox)
  {
    for(int i=0; i<D; ++i)
      ubox[i]=u[i]-std::floor(u[i]);
    return true;
  }
};

template<class T>
struct CheckBoxConds<T,1>
{
  inline static bool inside(const TinyVector<T,1>& u)
  {
    return (u[0]>0.0 && u[0]<1.0);
  }

  inline static bool inside(const TinyVector<T,1>& u, TinyVector<T,1>& ubox)
  {
    ubox[0]=u[0]-std::floor(u[0]);
    return true;
  }
};

template<class T>
struct CheckBoxConds<T,2>
{
  inline static bool inside(const TinyVector<T,2>& u)
  {
    return (u[0]>0.0 && u[0]<1.0) && (u[1]>0.0 && u[1]<1.0);
  }
  inline static bool inside(const TinyVector<T,2>& u, TinyVector<T,2>& ubox)
  {
    ubox[0]=u[0]-std::floor(u[0]);
    ubox[1]=u[1]-std::floor(u[1]);
    return true;
  }
};


template<class T>
struct CheckBoxConds<T,3>
{
  inline static bool inside(const TinyVector<T,3>& u)
  {
    return (u[0]>0.0 && u[0]<1.0) &&
           (u[1]>0.0 && u[1]<1.0) &&
           (u[2]>0.0 && u[2]<1.0);
  }

  inline static bool inside(const TinyVector<T,3>& u, const TinyVector<int,3>& bc)
  {
    return
      (bc[0] || (u[0]>0.0 && u[0]<1.0)) &&
      (bc[1] || (u[1]>0.0 && u[1]<1.0)) &&
      (bc[2] || (u[2]>0.0 && u[2]<1.0));
  }

  inline static bool inside(const TinyVector<T,3>& u, TinyVector<T,3>& ubox)
  {
    ubox[0]=u[0]-std::floor(u[0]);
    ubox[1]=u[1]-std::floor(u[1]);
    ubox[2]=u[2]-std::floor(u[2]);
    return true;
  }
};

}

#endif

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
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//   Department of Physics, Ohio State University
//   Ohio Supercomputer Center
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef OHMMS_LATTICEOPERATIONS_H
#define OHMMS_LATTICEOPERATIONS_H
#include "OhmmsPETE/TinyVector.h"
#include "OhmmsPETE/Tensor.h"

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
  apply(const TinyVector<T1,D>& lhs, const Tensor<T2,D>& rhs) {
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
  apply(const TinyVector<T1,3>& lhs, const Tensor<T2,3>& rhs) {
    return TinyVector<Type_t,3>(lhs[0]*rhs[0]+lhs[1]*rhs[3]+lhs[2]*rhs[6],
                                lhs[0]*rhs[1]+lhs[1]*rhs[4]+lhs[2]*rhs[7],    
                                lhs[0]*rhs[2]+lhs[1]*rhs[5]+lhs[2]*rhs[8]);
  }

  inline static TinyVector<Type_t,3>
  apply( const Tensor<T2,3>& lhs, const TinyVector<T1,3>& rhs) {
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
  apply(const TinyVector<T1,D>& lhs, const Tensor<T2,D>& rhs) {
    TinyVector<Type_t,D> ret;
    for(unsigned d=0; d<D; ++d) ret[d]=lhs[d]*rhs(d,d);
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
  apply(const TinyVector<T1,3>& lhs, const Tensor<T2,3>& rhs) {
    return TinyVector<Type_t,3>(lhs[0]*rhs[0],lhs[1]*rhs[4],lhs[2]*rhs[8]);
  }

  inline static TinyVector<Type_t,3>
  apply(const Tensor<T2,3>& lhs, const TinyVector<T1,3>& rhs) {
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
  apply(const TinyVector<T1,D>& a, const Tensor<T2,D>& X, const TinyVector<T1,D>& b) {
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
  apply(const TinyVector<T1,3>& a, const Tensor<T2,3>& X, const TinyVector<T1,3>& b) {
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
  apply(const TinyVector<T1,D>& a, const Tensor<T2,D>& X, const TinyVector<T1,D>& b) {
    Type_t res = a[0]*X(0,0)*b[0];
    for(int d=1;d<D; d++) res += a[d]*X(d,d)*b[d];
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
  apply(const TinyVector<T1,3>& a, const Tensor<T2,3>& X, const TinyVector<T1,3>& b) {
    return a[0]*X[0]*b[0]+a[1]*X[4]*b[1]+a[2]*X[8]*b[2];
  }
};

#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

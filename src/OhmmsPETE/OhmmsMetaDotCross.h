//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign   
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign 
//////////////////////////////////////////////////////////////////////////////////////


#ifndef OHMMS_TINYVECTORTENSOR_DOTCROSS_H
#define OHMMS_TINYVECTORTENSOR_DOTCROSS_H
namespace qmcplusplus
{
//////////////////////////////////////////////////////////////////////
//
// Definition of the struct OTDot.
//
//////////////////////////////////////////////////////////////////////

template<class T1, class T2> struct OTDot {};

//////////////////////////////////////////////////////////////////////
//
// Specializations for TinyVector dot TinyVector
//
//////////////////////////////////////////////////////////////////////

template<class T1, class T2, unsigned D>
struct OTDot< TinyVector<T1,D> , TinyVector<T2,D> >
{
  typedef typename BinaryReturn<T1,T2,OpMultiply>::Type_t Type_t;
  inline static Type_t
  apply(const TinyVector<T1,D>& lhs, const TinyVector<T2,D>& rhs)
  {
    Type_t dot = lhs[0]*rhs[0];
    for (unsigned d=1; d<D; ++d)
      dot += lhs[d]*rhs[d];
    return dot;
  }
};

template<class T1, class T2>
struct OTDot< TinyVector<T1,1> , TinyVector<T2,1> >
{
  typedef typename BinaryReturn<T1,T2,OpMultiply>::Type_t Type_t;
  inline static Type_t
  apply(const TinyVector<T1,1>& lhs, const TinyVector<T2,1>& rhs)
  {
    return lhs[0]*rhs[0];
  }
};

template<class T1, class T2>
struct OTDot< TinyVector<T1,2> , TinyVector<T2,2> >
{
  typedef typename BinaryReturn<T1,T2,OpMultiply>::Type_t Type_t;
  inline static Type_t
  apply(const TinyVector<T1,2>& lhs, const TinyVector<T2,2>& rhs)
  {
    return lhs[0]*rhs[0] + lhs[1]*rhs[1];
  }
};

template<class T1, class T2>
struct OTDot< TinyVector<T1,3> , TinyVector<T2,3> >
{
  typedef typename BinaryReturn<T1,T2,OpMultiply>::Type_t Type_t;
  inline static Type_t
  apply(const TinyVector<T1,3>& lhs, const TinyVector<T2,3>& rhs)
  {
    return lhs[0]*rhs[0] + lhs[1]*rhs[1] + lhs[2]*rhs[2];
  }
};


//////////////////////////////////////////////////////////////////////
//
// Definition of the struct OTCross.
//
//////////////////////////////////////////////////////////////////////

template<class T1, class T2> struct OTCross {};

//////////////////////////////////////////////////////////////////////
//
// Specializations for TinyVector cross TinyVector
//
//////////////////////////////////////////////////////////////////////

template<class T1, class T2, unsigned D>
struct OTCross< TinyVector<T1,D> , TinyVector<T2,D> >
{
  typedef typename BinaryReturn<T1,T2,OpMultiply>::Type_t Type_t;
  inline static TinyVector<Type_t,D>
  apply(const TinyVector<T1,D>& a, const TinyVector<T2,D>& b)
  {
    ERRORMSG("Cross-product *only* implemented for 3D; you're trying to"
             << " do it for " << D << "D." << std::endl);
    TinyVector<Type_t,D> bogusCross(-99999);
    return bogusCross;
  }
};

template<class T1, class T2>
struct OTCross< TinyVector<T1,3> , TinyVector<T2,3> >
{
  typedef typename BinaryReturn<T1,T2,OpMultiply>::Type_t Type_t;
  inline static TinyVector<Type_t,3>
  apply(const TinyVector<T1,3>& a, const TinyVector<T2,3>& b)
  {
    TinyVector<Type_t,3> cross;
    cross[0] = a[1]*b[2] - a[2]*b[1];
    cross[1] = a[2]*b[0] - a[0]*b[2];
    cross[2] = a[0]*b[1] - a[1]*b[0];
    return cross;
  }
};

//////////////////////////////////////////////////////////////////////
//
// Specializations for Tensor dot Tensor
//
//////////////////////////////////////////////////////////////////////

template<class T1, class T2, unsigned D>
struct OTDot< Tensor<T1,D> , Tensor<T2,D> >
{
  typedef typename BinaryReturn<T1,T2,OpMultiply>::Type_t Type_t;
  inline static Tensor<Type_t,D>
  apply(const Tensor<T1,D>& lhs, const Tensor<T2,D>& rhs)
  {
    Tensor<Type_t,D> dot = Tensor<Type_t,D>::DontInitialize();
    for (unsigned int i=0; i<D; ++i)
      for (unsigned int j=0; j<D; ++j)
      {
        Type_t sum = lhs(i,0) * rhs(0,j);
        for (unsigned int k=1; k<D; ++k)
          sum += lhs(i,k) * rhs(k,j);
        dot(i,j) = sum;
      }
    return dot;
  }
};

template<class T1, class T2>
struct OTDot< Tensor<T1,1> , Tensor<T2,1> >
{
  typedef typename BinaryReturn<T1,T2,OpMultiply>::Type_t Type_t;
  inline static Tensor<Type_t,1>
  apply(const Tensor<T1,1>& lhs, const Tensor<T2,1>& rhs)
  {
    return Tensor<Type_t,1>(lhs[0]*rhs[0]);
  }
};

template<class T1, class T2>
struct OTDot< Tensor<T1,2> , Tensor<T2,2> >
{
  typedef typename BinaryReturn<T1,T2,OpMultiply>::Type_t Type_t;
  inline static Tensor<Type_t,2>
  apply(const Tensor<T1,2>& lhs, const Tensor<T2,2>& rhs)
  {
    return Tensor<Type_t,2>(lhs(0,0)*rhs(0,0) + lhs(0,1)*rhs(1,0),
                            lhs(0,0)*rhs(0,1) + lhs(0,1)*rhs(1,1),
                            lhs(1,0)*rhs(0,0) + lhs(1,1)*rhs(1,0),
                            lhs(1,0)*rhs(0,1) + lhs(1,1)*rhs(1,1));
  }
};

template<class T1, class T2>
struct OTDot< Tensor<T1,3> , Tensor<T2,3> >
{
  typedef typename BinaryReturn<T1,T2,OpMultiply>::Type_t Type_t;
  inline static Tensor<Type_t,3>
  apply(const Tensor<T1,3>& lhs, const Tensor<T2,3>& rhs)
  {
    return Tensor<Type_t,3>( lhs(0,0)*rhs(0,0) + lhs(0,1)*rhs(1,0) + lhs(0,2)*rhs(2,0) ,
                             lhs(0,0)*rhs(0,1) + lhs(0,1)*rhs(1,1) + lhs(0,2)*rhs(2,1) ,
                             lhs(0,0)*rhs(0,2) + lhs(0,1)*rhs(1,2) + lhs(0,2)*rhs(2,2) ,
                             lhs(1,0)*rhs(0,0) + lhs(1,1)*rhs(1,0) + lhs(1,2)*rhs(2,0) ,
                             lhs(1,0)*rhs(0,1) + lhs(1,1)*rhs(1,1) + lhs(1,2)*rhs(2,1) ,
                             lhs(1,0)*rhs(0,2) + lhs(1,1)*rhs(1,2) + lhs(1,2)*rhs(2,2) ,
                             lhs(2,0)*rhs(0,0) + lhs(2,1)*rhs(1,0) + lhs(2,2)*rhs(2,0) ,
                             lhs(2,0)*rhs(0,1) + lhs(2,1)*rhs(1,1) + lhs(2,2)*rhs(2,1) ,
                             lhs(2,0)*rhs(0,2) + lhs(2,1)*rhs(1,2) + lhs(2,2)*rhs(2,2) );
  }
};

//////////////////////////////////////////////////////////////////////
//
// Specializations for SymTensor dot SymTensor
//
//////////////////////////////////////////////////////////////////////

template<class T1, class T2, unsigned D>
struct OTDot< SymTensor<T1,D> , SymTensor<T2,D> >
{
  typedef typename BinaryReturn<T1,T2,OpMultiply>::Type_t Type_t;
  inline static Tensor<Type_t,D>
  apply(const SymTensor<T1,D>& lhs, const SymTensor<T2,D>& rhs)
  {
    Tensor<Type_t,D> dot = Tensor<Type_t,D>::DontInitialize();
    for (unsigned int i=0; i<D; ++i)
      for (unsigned int j=i; j<D; ++j)
      {
        Type_t sum = lhs.HL(i,0) * rhs.HL(j,0);
        int k=1;
        for ( ; k<i; ++k )
          sum += lhs.HL(i,k) * rhs.HL(j,k);
        for ( ; k<j; ++k )
          sum += lhs.HL(k,i) * rhs.HL(j,k);
        for ( ; k<D; ++k )
          sum += lhs.HL(k,i) * rhs.HL(k,j);
        dot(i,j) = sum;
      }
    return dot;
  }
};

template<class T1, class T2>
struct OTDot< SymTensor<T1,1> , SymTensor<T2,1> >
{
  typedef typename BinaryReturn<T1,T2,OpMultiply>::Type_t Type_t;
  inline static Tensor<Type_t,1>
  apply(const SymTensor<T1,1>& lhs, const SymTensor<T2,1>& rhs)
  {
    return Tensor<Type_t,1>(lhs[0]*rhs[0]);
  }
};

template<class T1, class T2>
struct OTDot< SymTensor<T1,2> , SymTensor<T2,2> >
{
  typedef typename BinaryReturn<T1,T2,OpMultiply>::Type_t Type_t;
  inline static Tensor<Type_t,2>
  apply(const SymTensor<T1,2>& lhs, const SymTensor<T2,2>& rhs)
  {
    return Tensor<Type_t,2>(lhs(0,0)*rhs(0,0) + lhs(0,1)*rhs(1,0),
                            lhs(0,0)*rhs(0,1) + lhs(0,1)*rhs(1,1),
                            lhs(1,0)*rhs(0,0) + lhs(1,1)*rhs(1,0),
                            lhs(1,0)*rhs(0,1) + lhs(1,1)*rhs(1,1));
  }
};

template<class T1, class T2>
struct OTDot< SymTensor<T1,3> , SymTensor<T2,3> >
{
  typedef typename BinaryReturn<T1,T2,OpMultiply>::Type_t Type_t;
  inline static Tensor<Type_t,3>
  apply(const SymTensor<T1,3>& lhs, const SymTensor<T2,3>& rhs)
  {
    return Tensor<Type_t,3>( lhs(0,0)*rhs(0,0) + lhs(0,1)*rhs(1,0) + lhs(0,2)*rhs(2,0) ,
                             lhs(0,0)*rhs(0,1) + lhs(0,1)*rhs(1,1) + lhs(0,2)*rhs(2,1) ,
                             lhs(0,0)*rhs(0,2) + lhs(0,1)*rhs(1,2) + lhs(0,2)*rhs(2,2) ,
                             lhs(1,0)*rhs(0,0) + lhs(1,1)*rhs(1,0) + lhs(1,2)*rhs(2,0) ,
                             lhs(1,0)*rhs(0,1) + lhs(1,1)*rhs(1,1) + lhs(1,2)*rhs(2,1) ,
                             lhs(1,0)*rhs(0,2) + lhs(1,1)*rhs(1,2) + lhs(1,2)*rhs(2,2) ,
                             lhs(2,0)*rhs(0,0) + lhs(2,1)*rhs(1,0) + lhs(2,2)*rhs(2,0) ,
                             lhs(2,0)*rhs(0,1) + lhs(2,1)*rhs(1,1) + lhs(2,2)*rhs(2,1) ,
                             lhs(2,0)*rhs(0,2) + lhs(2,1)*rhs(1,2) + lhs(2,2)*rhs(2,2) );
  }
};

//////////////////////////////////////////////////////////////////////
//
// Specializations for Tensor dot TinyVector
//
//////////////////////////////////////////////////////////////////////

template<class T1, class T2, unsigned D>
struct OTDot< Tensor<T1,D> , TinyVector<T2,D> >
{
  typedef typename BinaryReturn<T1,T2,OpMultiply>::Type_t Type_t;
  inline static TinyVector<Type_t,D>
  apply(const Tensor<T1,D>& lhs, const TinyVector<T2,D>& rhs)
  {
    TinyVector<Type_t,D> ret= TinyVector<Type_t,D>::DontInitialize();
    for (unsigned int i=0; i<D; ++i)
    {
      Type_t sum = lhs(i,0)*rhs[0];
      for (unsigned int j=1; j<D; ++j)
        sum += lhs(i,j)*rhs[j];
      ret[i] = sum;
    }
    return ret;
  }
};


template<class T1, class T2>
struct OTDot< Tensor<T1,1> , TinyVector<T2,1> >
{
  typedef typename BinaryReturn<T1,T2,OpMultiply>::Type_t Type_t;
  inline static TinyVector<Type_t,1>
  apply(const Tensor<T1,1>& lhs, const TinyVector<T2,1>& rhs)
  {
    return TinyVector<Type_t,1>( lhs[0]*rhs[0] );
  }
};

template<class T1, class T2>
struct OTDot< Tensor<T1,2> , TinyVector<T2,2> >
{
  typedef typename BinaryReturn<T1,T2,OpMultiply>::Type_t Type_t;
  inline static TinyVector<Type_t,2>
  apply(const Tensor<T1,2>& lhs, const TinyVector<T2,2>& rhs)
  {
    return TinyVector<Type_t,2>( lhs(0,0)*rhs[0] + lhs(0,1)*rhs[1] ,
                                 lhs(1,0)*rhs[0] + lhs(1,1)*rhs[1] );
  }
};

template<class T1, class T2>
struct OTDot< Tensor<T1,3> , TinyVector<T2,3> >
{
  typedef typename BinaryReturn<T1,T2,OpMultiply>::Type_t Type_t;
  inline static TinyVector<Type_t,3>
  apply(const Tensor<T1,3>& lhs, const TinyVector<T2,3>& rhs)
  {
    return TinyVector<Type_t,3>( lhs(0,0)*rhs[0] + lhs(0,1)*rhs[1] + lhs(0,2)*rhs[2],
                                 lhs(1,0)*rhs[0] + lhs(1,1)*rhs[1] + lhs(1,2)*rhs[2],
                                 lhs(2,0)*rhs[0] + lhs(2,1)*rhs[1] + lhs(2,2)*rhs[2] );
  }
};

//////////////////////////////////////////////////////////////////////
//
// Specializations for TinyVector dot Tensor
//
//////////////////////////////////////////////////////////////////////

template<class T1, class T2, unsigned D>
struct OTDot< TinyVector<T1,D> , Tensor<T2,D> >
{
  typedef typename BinaryReturn<T1,T2,OpMultiply>::Type_t Type_t;
  inline static TinyVector<Type_t,D>
  apply(const TinyVector<T2,D>& lhs, const Tensor<T1,D>& rhs)
  {
    TinyVector<Type_t,D> ret = TinyVector<Type_t,D>::DontInitialize();
    for (unsigned int i=0; i<D; ++i)
    {
      Type_t sum = lhs[0]*rhs(0,i);
      for (unsigned int j=1; j<D; ++j)
        sum += lhs[j]*rhs(j,i);
      ret[i] = sum;
    }
    return ret;
  }
};


template<class T1, class T2>
struct OTDot< TinyVector<T1,1> , Tensor<T2,1> >
{
  typedef typename BinaryReturn<T1,T2,OpMultiply>::Type_t Type_t;
  inline static TinyVector<Type_t,1>
  apply(const TinyVector<T1,1>& lhs, const Tensor<T2,1>& rhs)
  {
    return TinyVector<Type_t,1>( lhs[0]*rhs[0] );
  }
};

template<class T1, class T2>
struct OTDot< TinyVector<T1,2> , Tensor<T2,2> >
{
  typedef typename BinaryReturn<T1,T2,OpMultiply>::Type_t Type_t;
  inline static TinyVector<Type_t,2>
  apply(const TinyVector<T1,2>& lhs, const Tensor<T2,2>& rhs)
  {
    return TinyVector<Type_t,2>( lhs[0]*rhs(0,0) + lhs[1]*rhs(1,0) ,
                                 lhs[0]*rhs(0,1) + lhs[1]*rhs(1,1) );
  }
};

template<class T1, class T2>
struct OTDot< TinyVector<T1,3> , Tensor<T2,3> >
{
  typedef typename BinaryReturn<T1,T2,OpMultiply>::Type_t Type_t;
  inline static TinyVector<Type_t,3>
  apply(const TinyVector<T1,3>& lhs, const Tensor<T2,3>& rhs)
  {
    return TinyVector<Type_t,3>( lhs[0]*rhs(0,0) + lhs[1]*rhs(1,0) + lhs[2]*rhs(2,0),
                                 lhs[0]*rhs(0,1) + lhs[1]*rhs(1,1) + lhs[2]*rhs(2,1),
                                 lhs[0]*rhs(0,2) + lhs[1]*rhs(1,2) + lhs[2]*rhs(2,2) );
  }
};

//////////////////////////////////////////////////////////////////////
//
// Specializations for SymTensor dot TinyVector
//
//////////////////////////////////////////////////////////////////////

template<class T1, class T2, unsigned D>
struct OTDot< SymTensor<T1,D> , TinyVector<T2,D> >
{
  typedef typename BinaryReturn<T1,T2,OpMultiply>::Type_t Type_t;
  inline static TinyVector<Type_t,D>
  apply(const SymTensor<T1,D>& lhs, const TinyVector<T2,D>& rhs)
  {
    TinyVector<Type_t,D> ret= TinyVector<Type_t,D>::DontInitialize();
    for (unsigned int i=0; i<D; ++i)
    {
      Type_t sum = lhs.HL(i,0)*rhs[0];
      unsigned int j=1;
      for ( ; j<i; ++j)
        sum += lhs.HL(i,j)*rhs[j];
      for ( ; j<D; ++j)
        sum += lhs.HL(j,i)*rhs[j];
      ret[i] = sum;
    }
    return ret;
  }
};


template<class T1, class T2>
struct OTDot< SymTensor<T1,1> , TinyVector<T2,1> >
{
  typedef typename BinaryReturn<T1,T2,OpMultiply>::Type_t Type_t;
  inline static TinyVector<Type_t,1>
  apply(const SymTensor<T1,1>& lhs, const TinyVector<T2,1>& rhs)
  {
    return TinyVector<Type_t,1>( lhs[0]*rhs[0] );
  }
};

template<class T1, class T2>
struct OTDot< SymTensor<T1,2> , TinyVector<T2,2> >
{
  typedef typename BinaryReturn<T1,T2,OpMultiply>::Type_t Type_t;
  inline static TinyVector<Type_t,2>
  apply(const SymTensor<T1,2>& lhs, const TinyVector<T2,2>& rhs)
  {
    return TinyVector<Type_t,2>( lhs(0,0)*rhs[0] + lhs(0,1)*rhs[1] ,
                                 lhs(1,0)*rhs[0] + lhs(1,1)*rhs[1] );
  }
};

template<class T1, class T2>
struct OTDot< SymTensor<T1,3> , TinyVector<T2,3> >
{
  typedef typename BinaryReturn<T1,T2,OpMultiply>::Type_t Type_t;
  inline static TinyVector<Type_t,3>
  apply(const SymTensor<T1,3>& lhs, const TinyVector<T2,3>& rhs)
  {
    return TinyVector<Type_t,3>( lhs(0,0)*rhs[0] + lhs(0,1)*rhs[1] + lhs(0,2)*rhs[2],
                                 lhs(1,0)*rhs[0] + lhs(1,1)*rhs[1] + lhs(1,2)*rhs[2],
                                 lhs(2,0)*rhs[0] + lhs(2,1)*rhs[1] + lhs(2,2)*rhs[2] );
  }
};

//////////////////////////////////////////////////////////////////////
//
// Specializations for TinyVector dot SymTensor
//
//////////////////////////////////////////////////////////////////////

template<class T1, class T2, unsigned D>
struct OTDot< TinyVector<T1,D> , SymTensor<T2,D> >
{
  typedef typename BinaryReturn<T1,T2,OpMultiply>::Type_t Type_t;
  inline static TinyVector<Type_t,D>
  apply(const TinyVector<T2,D>& lhs, const SymTensor<T1,D>& rhs)
  {
    TinyVector<Type_t,D> ret = TinyVector<Type_t,D>::DontInitialize();
    for (unsigned int i=0; i<D; ++i)
    {
      Type_t sum = lhs[0]*rhs[i*(i+1)/2];
      unsigned int j=1;
      for ( ; j<i; ++j)
        sum += lhs[j]*rhs[i*(i+1)/2+j];
      for ( ; j<D; ++j)
        sum += lhs[j]*rhs[j*(j+1)/2+i];
      ret[i] = sum;
    }
    return ret;
  }
};


template<class T1, class T2>
struct OTDot< TinyVector<T1,1> , SymTensor<T2,1> >
{
  typedef typename BinaryReturn<T1,T2,OpMultiply>::Type_t Type_t;
  inline static TinyVector<Type_t,1>
  apply(const TinyVector<T1,1>& lhs, const SymTensor<T2,1>& rhs)
  {
    return TinyVector<Type_t,1>( lhs[0]*rhs[0] );
  }
};

template<class T1, class T2>
struct OTDot< TinyVector<T1,2> , SymTensor<T2,2> >
{
  typedef typename BinaryReturn<T1,T2,OpMultiply>::Type_t Type_t;
  inline static TinyVector<Type_t,2>
  apply(const TinyVector<T1,2>& lhs, const SymTensor<T2,2>& rhs)
  {
    return TinyVector<Type_t,2>( lhs[0]*rhs(0,0) + lhs[1]*rhs(1,0) ,
                                 lhs[0]*rhs(0,1) + lhs[1]*rhs(1,1) );
  }
};

template<class T1, class T2>
struct OTDot< TinyVector<T1,3> , SymTensor<T2,3> >
{
  typedef typename BinaryReturn<T1,T2,OpMultiply>::Type_t Type_t;
  inline static TinyVector<Type_t,3>
  apply(const TinyVector<T1,3>& lhs, const SymTensor<T2,3>& rhs)
  {
    return TinyVector<Type_t,3>( lhs[0]*rhs(0,0) + lhs[1]*rhs(1,0) + lhs[2]*rhs(2,0),
                                 lhs[0]*rhs(0,1) + lhs[1]*rhs(1,1) + lhs[2]*rhs(2,1),
                                 lhs[0]*rhs(0,2) + lhs[1]*rhs(1,2) + lhs[2]*rhs(2,2) );
  }
};

//////////////////////////////////////////////////////////////////////
//
// Specializations for SymTensor dot Tensor
//
//////////////////////////////////////////////////////////////////////

template<class T1, class T2, unsigned D>
struct OTDot< SymTensor<T1,D> , Tensor<T2,D> >
{
  typedef typename BinaryReturn<T1,T2,OpMultiply>::Type_t Type_t;
  inline static Tensor<Type_t,D>
  apply(const SymTensor<T1,D>& lhs, const Tensor<T2,D>& rhs)
  {
    Tensor<Type_t,D> dot = Tensor<Type_t,D>::DontInitialize();
    for (unsigned int i=0; i<D; ++i)
      for (unsigned int j=0; j<D; ++j)
      {
        Type_t sum = lhs.HL(i,0) * rhs(0,j);
        unsigned int k = 1;
        for (; k<i; ++k)
          sum += lhs.HL(i,k) * rhs(k,j);
        for (; k<D; ++k)
          sum += lhs.HL(k,i) * rhs(k,j);
        dot(i,j) = sum;
      }
    return dot;
  }
};

template<class T1, class T2>
struct OTDot< SymTensor<T1,1> , Tensor<T2,1> >
{
  typedef typename BinaryReturn<T1,T2,OpMultiply>::Type_t Type_t;
  inline static Tensor<Type_t,1>
  apply(const SymTensor<T1,1>& lhs, const Tensor<T2,1>& rhs)
  {
    return Tensor<Type_t,1>(lhs[0]*rhs[0]);
  }
};

template<class T1, class T2>
struct OTDot< SymTensor<T1,2> , Tensor<T2,2> >
{
  typedef typename BinaryReturn<T1,T2,OpMultiply>::Type_t Type_t;
  inline static Tensor<Type_t,2>
  apply(const SymTensor<T1,2>& lhs, const Tensor<T2,2>& rhs)
  {
    return Tensor<Type_t,2>(lhs(0,0)*rhs(0,0) + lhs(0,1)*rhs(1,0),
                            lhs(0,0)*rhs(0,1) + lhs(0,1)*rhs(1,1),
                            lhs(1,0)*rhs(0,0) + lhs(1,1)*rhs(1,0),
                            lhs(1,0)*rhs(0,1) + lhs(1,1)*rhs(1,1));
  }
};

template<class T1, class T2>
struct OTDot< SymTensor<T1,3> , Tensor<T2,3> >
{
  typedef typename BinaryReturn<T1,T2,OpMultiply>::Type_t Type_t;
  inline static Tensor<Type_t,3>
  apply(const SymTensor<T1,3>& lhs, const Tensor<T2,3>& rhs)
  {
    return Tensor<Type_t,3>( lhs(0,0)*rhs(0,0) + lhs(0,1)*rhs(1,0) + lhs(0,2)*rhs(2,0) ,
                             lhs(0,0)*rhs(0,1) + lhs(0,1)*rhs(1,1) + lhs(0,2)*rhs(2,1) ,
                             lhs(0,0)*rhs(0,2) + lhs(0,1)*rhs(1,2) + lhs(0,2)*rhs(2,2) ,
                             lhs(1,0)*rhs(0,0) + lhs(1,1)*rhs(1,0) + lhs(1,2)*rhs(2,0) ,
                             lhs(1,0)*rhs(0,1) + lhs(1,1)*rhs(1,1) + lhs(1,2)*rhs(2,1) ,
                             lhs(1,0)*rhs(0,2) + lhs(1,1)*rhs(1,2) + lhs(1,2)*rhs(2,2) ,
                             lhs(2,0)*rhs(0,0) + lhs(2,1)*rhs(1,0) + lhs(2,2)*rhs(2,0) ,
                             lhs(2,0)*rhs(0,1) + lhs(2,1)*rhs(1,1) + lhs(2,2)*rhs(2,1) ,
                             lhs(2,0)*rhs(0,2) + lhs(2,1)*rhs(1,2) + lhs(2,2)*rhs(2,2) );
  }
};

//////////////////////////////////////////////////////////////////////
//
// Specializations for Tensor dot SymTensor
//
//////////////////////////////////////////////////////////////////////

template<class T1, class T2, unsigned D>
struct OTDot< Tensor<T1,D> , SymTensor<T2,D> >
{
  typedef typename BinaryReturn<T1,T2,OpMultiply>::Type_t Type_t;
  inline static Tensor<Type_t,D>
  apply(const Tensor<T1,D>& lhs, const SymTensor<T2,D>& rhs)
  {
    Tensor<Type_t,D> dot = Tensor<Type_t,D>::DontInitialize();
    for (unsigned int i=0; i<D; ++i)
      for (unsigned int j=0; j<D; ++j)
      {
        Type_t sum = lhs(i,0) * rhs.HL(j,0);
        unsigned int k=1;
        for (; k<j; ++k)
          sum += lhs(i,k) * rhs.HL(j,k);
        for (; k<D; ++k)
          sum += lhs(i,k) * rhs.HL(k,j);
        dot[i*D+j] = sum;
      }
    return dot;
  }
};

template<class T1, class T2>
struct OTDot< Tensor<T1,1> , SymTensor<T2,1> >
{
  typedef typename BinaryReturn<T1,T2,OpMultiply>::Type_t Type_t;
  inline static Tensor<Type_t,1>
  apply(const Tensor<T1,1>& lhs, const SymTensor<T2,1>& rhs)
  {
    return Tensor<Type_t,1>(lhs[0]*rhs[0]);
  }
};

template<class T1, class T2>
struct OTDot< Tensor<T1,2> , SymTensor<T2,2> >
{
  typedef typename BinaryReturn<T1,T2,OpMultiply>::Type_t Type_t;
  inline static Tensor<Type_t,2>
  apply(const Tensor<T1,2>& lhs, const SymTensor<T2,2>& rhs)
  {
    return Tensor<Type_t,2>(lhs(0,0)*rhs(0,0) + lhs(0,1)*rhs(1,0),
                            lhs(0,0)*rhs(0,1) + lhs(0,1)*rhs(1,1),
                            lhs(1,0)*rhs(0,0) + lhs(1,1)*rhs(1,0),
                            lhs(1,0)*rhs(0,1) + lhs(1,1)*rhs(1,1));
  }
};

template<class T1, class T2>
struct OTDot< Tensor<T1,3> , SymTensor<T2,3> >
{
  typedef typename BinaryReturn<T1,T2,OpMultiply>::Type_t Type_t;
  inline static Tensor<Type_t,3>
  apply(const Tensor<T1,3>& lhs, const SymTensor<T2,3>& rhs)
  {
    return Tensor<Type_t,3>( lhs(0,0)*rhs(0,0) + lhs(0,1)*rhs(1,0) + lhs(0,2)*rhs(2,0) ,
                             lhs(0,0)*rhs(0,1) + lhs(0,1)*rhs(1,1) + lhs(0,2)*rhs(2,1) ,
                             lhs(0,0)*rhs(0,2) + lhs(0,1)*rhs(1,2) + lhs(0,2)*rhs(2,2) ,
                             lhs(1,0)*rhs(0,0) + lhs(1,1)*rhs(1,0) + lhs(1,2)*rhs(2,0) ,
                             lhs(1,0)*rhs(0,1) + lhs(1,1)*rhs(1,1) + lhs(1,2)*rhs(2,1) ,
                             lhs(1,0)*rhs(0,2) + lhs(1,1)*rhs(1,2) + lhs(1,2)*rhs(2,2) ,
                             lhs(2,0)*rhs(0,0) + lhs(2,1)*rhs(1,0) + lhs(2,2)*rhs(2,0) ,
                             lhs(2,0)*rhs(0,1) + lhs(2,1)*rhs(1,1) + lhs(2,2)*rhs(2,1) ,
                             lhs(2,0)*rhs(0,2) + lhs(2,1)*rhs(1,2) + lhs(2,2)*rhs(2,2) );
  }
};

//////////////////////////////////////////////////////////////////////
//
// Specializations for TinyVector dot AntiSymTensor
//
//////////////////////////////////////////////////////////////////////

template<class T1, class T2, unsigned D>
struct OTDot< TinyVector<T1,D> , AntiSymTensor<T2,D> >
{
  typedef typename BinaryReturn<T1,T2,OpMultiply>::Type_t Type_t;
  inline static TinyVector<Type_t,D>
  apply(const TinyVector<T2,D>& lhs, const AntiSymTensor<T1,D>& rhs)
  {
    TinyVector<Type_t,D> ret = TinyVector<Type_t,D>::DontInitialize();
    for (unsigned int j=0; j<D; ++j)
    {
      double sum = 0;
      for ( i=0; i<j; i++)
        sum -= lhs[i]*rhs[((j-1)*j/2)+i];
      for ( i=j+1; i<D; ++i)
        sum += lhs[i]*rhs[((i-1)*i/2)+j];
      ret[j] = sum;
    }
    return ret;
  }
};

template<class T1, class T2>
struct OTDot< TinyVector<T1,2> , AntiSymTensor<T2,2> >
{
  typedef typename BinaryReturn<T1,T2,OpMultiply>::Type_t Type_t;
  inline static TinyVector<Type_t,2>
  apply(const TinyVector<T1,2>& lhs, const AntiSymTensor<T2,2>& rhs)
  {
    return TinyVector<Type_t,2>( lhs[0]*rhs(0,0) + lhs[1]*rhs(1,0) ,
                                 lhs[0]*rhs(0,1) + lhs[1]*rhs(1,1) );
  }
};

template<class T1, class T2>
struct OTDot< TinyVector<T1,3> , AntiSymTensor<T2,3> >
{
  typedef typename BinaryReturn<T1,T2,OpMultiply>::Type_t Type_t;
  inline static TinyVector<Type_t,3>
  apply(const TinyVector<T1,3>& lhs, const AntiSymTensor<T2,3>& rhs)
  {
    return TinyVector<Type_t,3>( lhs[0]*rhs(0,0) + lhs[1]*rhs(1,0) + lhs[2]*rhs(2,0),
                                 lhs[0]*rhs(0,1) + lhs[1]*rhs(1,1) + lhs[2]*rhs(2,1),
                                 lhs[0]*rhs(0,2) + lhs[1]*rhs(1,2) + lhs[2]*rhs(2,2) );
  }
};

//////////////////////////////////////////////////////////////////////
//
// Specializations for AntiSymTensor dot TinyVector
//
//////////////////////////////////////////////////////////////////////

template<class T1, class T2, unsigned D>
struct OTDot< AntiSymTensor<T1,D> , TinyVector<T2,D> >
{
  typedef typename BinaryReturn<T1,T2,OpMultiply>::Type_t Type_t;
  inline static TinyVector<Type_t,D>
  apply(const AntiSymTensor<T1,D>& lhs, const TinyVector<T2,D>& rhs)
  {
    TinyVector<Type_t,D> ret= TinyVector<Type_t,D>::DontInitialize();
    for (unsigned int i=0; i<D; ++i)
    {
      Type_t sum = 0;
      for (unsigned int j=0; j<i; ++j)
        sum += lhs[((i-1)*i/2)+j]*rhs[j];
      for (unsigned int j=i+1; j<D; ++j)
        sum -= lhs[((j-1)*j/2)+i]*rhs[j];
      ret[i] = sum;
    }
    return ret;
  }
};


template<class T1, class T2>
struct OTDot< AntiSymTensor<T1,1> , TinyVector<T2,1> >
{
  typedef typename BinaryReturn<T1,T2,OpMultiply>::Type_t Type_t;
  inline static TinyVector<Type_t,1>
  apply(const AntiSymTensor<T1,1>& lhs, const TinyVector<T2,1>& rhs)
  {
    return TinyVector<Type_t,1>( lhs[0]*rhs[0] );
  }
};

template<class T1, class T2>
struct OTDot< AntiSymTensor<T1,2> , TinyVector<T2,2> >
{
  typedef typename BinaryReturn<T1,T2,OpMultiply>::Type_t Type_t;
  inline static TinyVector<Type_t,2>
  apply(const AntiSymTensor<T1,2>& lhs, const TinyVector<T2,2>& rhs)
  {
    return TinyVector<Type_t,2>( lhs(0,0)*rhs[0] + lhs(0,1)*rhs[1] ,
                                 lhs(1,0)*rhs[0] + lhs(1,1)*rhs[1] );
  }
};

template<class T1, class T2>
struct OTDot< AntiSymTensor<T1,3> , TinyVector<T2,3> >
{
  typedef typename BinaryReturn<T1,T2,OpMultiply>::Type_t Type_t;
  inline static TinyVector<Type_t,3>
  apply(const AntiSymTensor<T1,3>& lhs, const TinyVector<T2,3>& rhs)
  {
    return TinyVector<Type_t,3>( lhs(0,0)*rhs[0] + lhs(0,1)*rhs[1] + lhs(0,2)*rhs[2],
                                 lhs(1,0)*rhs[0] + lhs(1,1)*rhs[1] + lhs(1,2)*rhs[2],
                                 lhs(2,0)*rhs[0] + lhs(2,1)*rhs[1] + lhs(2,2)*rhs[2] );
  }
};
}

#endif // OHMMS_TINYVECTOR_DOTCROSS_H 


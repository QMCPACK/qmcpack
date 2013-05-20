//////////////////////////////////////////////////////////////////
// (c) Copyright 1998-2002 by Jeongnim Kim
//
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
/**@file ParticleAttribOps.h
 *@brief Declaraton of ParticleAttrib<T>
 */

/*! \class ParticleAttribOps
 *  \brief A one-dimensional vector class based on PETE.
 *
 *  Closely related to PETE STL vector example and written to substitute
 *  Poomar1::ParticleInteractAttrib.
 *
 *  Equivalent to blitz::Array<T,1>, pooma::Array<1,T>.
 *
 *  class C is a container class. Default is vector<T>
 *
 *  \todo Implement openMP compatible container class or evaluate function.
 *  \todo Implement get/put member functions for MPI-like parallelism
 */
#ifndef OHMMS_PARTICLEATTRIB_OPS_H
#define OHMMS_PARTICLEATTRIB_OPS_H

#include "ParticleBase/ParticleUtility.h"

namespace qmcplusplus
{

template<class T1, class T2, unsigned D>
struct OTCDot
{
  typedef typename BinaryReturn<T1,T2,OpMultiply>::Type_t Type_t;
  inline static Type_t
  apply(const TinyVector<std::complex<T1>,D>& lhs, const TinyVector<std::complex<T2>,D>& rhs)
  {
    Type_t res = lhs[0].real()*rhs[0].real()+lhs[0].imag()*rhs[0].imag();
    for (unsigned d=1; d<D; ++d)
      res += lhs[d].real()*rhs[d].real()-lhs[d].imag()*rhs[d].imag();
    return res;
  }
};

template<class T1, class T2>
struct OTCDot<T1,T2,1>
{
  typedef typename BinaryReturn<T1,T2,OpMultiply>::Type_t Type_t;
  inline static Type_t
  apply(const TinyVector<std::complex<T1>,3>& lhs, const TinyVector<std::complex<T2>,3>& rhs)
  {
    return lhs[0].real()*rhs[0].real()-lhs[0].imag()*rhs[0].imag();
  }
};

template<class T1, class T2>
struct OTCDot<T1,T2,2>
{
  typedef typename BinaryReturn<T1,T2,OpMultiply>::Type_t Type_t;
  inline static Type_t
  apply(const TinyVector<std::complex<T1>,3>& lhs, const TinyVector<std::complex<T2>,3>& rhs)
  {
    return lhs[0].real()*rhs[0].real()-lhs[0].imag()*rhs[0].imag()
           + lhs[1].real()*rhs[1].real()-lhs[1].imag()*rhs[1].imag();
  }
};

template<class T1, class T2>
struct OTCDot<T1,T2,3>
{
  typedef typename BinaryReturn<T1,T2,OpMultiply>::Type_t Type_t;
  inline static Type_t
  apply(const TinyVector<std::complex<T1>,3>& lhs, const
        TinyVector<std::complex<T2>,3>& rhs)
  {
    return lhs[0].real()*rhs[0].real()-lhs[0].imag()*rhs[0].imag()
           + lhs[1].real()*rhs[1].real()-lhs[1].imag()*rhs[1].imag()
           + lhs[2].real()*rhs[2].real()-lhs[2].imag()*rhs[2].imag();
  }
};

template<unsigned D>
inline double Dot(const ParticleAttrib<TinyVector<std::complex<double>, D> >& pa,
                  const ParticleAttrib<TinyVector<std::complex<double>, D> >& pb)
{
  double sum = 0;
  for(int i=0; i<pa.size(); i++)
  {
    sum += OTCDot<double,double,D>::apply(pa[i],pb[i]);
  }
  return sum;
}

inline double Sum(const ParticleAttrib<std::complex<double> >& pa)
{
  double sum = 0;
  for(int i=0; i<pa.size(); i++)
  {
    sum += pa[i].real();
  }
  return sum;
}

template<class T, unsigned D>
inline void Copy(const ParticleAttrib<TinyVector<std::complex<T>,D> >& c,
                 ParticleAttrib<TinyVector<T,D> >& r)
{
  for(int i=0; i<c.size(); i++)
  {
    //r[i]=real(c[i]);
    for(int j=0; j<D; j++)
      r[i][j] = c[i][j].real();
  }
}

template<class T, unsigned D>
inline void Copy(const ParticleAttrib<TinyVector<T,D> >& c,
                 ParticleAttrib<TinyVector<T,D> >& r)
{
  r=c;
}


/** generic PAOps
 */
template<class T, unsigned D> struct PAOps { };

///specialization for three-dimension
template<class T>
struct PAOps<T,3>
{

  typedef T                        real_type;
  typedef complex<T>               complex_type;
  typedef TinyVector<T,3>          rpos_type;
  typedef TinyVector<complex<T>,3> cpos_type;

  static inline
  void scale(T a, const ParticleAttrib<cpos_type>& pa,
             ParticleAttrib<rpos_type>& pb)
  {
    for(int i=0; i<pa.size(); i++)
    {
      pb[i][0]=a*pa[i][0].real();
      pb[i][1]=a*pa[i][1].real();
      pb[i][2]=a*pa[i][2].real();
    }
  }

  static inline
  void scale(T a, const ParticleAttrib<rpos_type>& pa,
             ParticleAttrib<rpos_type>& pb)
  {
    pb=a*pa;
  }

  static inline
  void axpy(T a, const ParticleAttrib<cpos_type>& pa, ParticleAttrib<rpos_type>& pb)
  {
    for(int i=0; i<pa.size(); ++i)
    {
      pb[i][0]+=a*pa[i][0].real();
      pb[i][1]+=a*pa[i][1].real();
      pb[i][2]+=a*pa[i][2].real();
    }
  }

  static inline
  void axpy(T a, const ParticleAttrib<rpos_type>& pa, ParticleAttrib<rpos_type>& pb)
  {
    for(int i=0; i<pa.size(); ++i)
    {
      pb[i][0]+=a*pa[i][0];
      pb[i][1]+=a*pa[i][1];
      pb[i][2]+=a*pa[i][2];
    }
  }

  static inline
  void axpy(T a, const ParticleAttrib<cpos_type>& pa, const ParticleAttrib<rpos_type>& pb,
            ParticleAttrib<rpos_type>& py)
  {
    for(int i=0; i<pa.size(); ++i)
    {
      py[i][0]=a*pa[i][0].real()+pb[i][0];
      py[i][1]=a*pa[i][1].real()+pb[i][1];
      py[i][2]=a*pa[i][2].real()+pb[i][2];
    }
  }

  static inline
  void axpy(T a, const ParticleAttrib<rpos_type>& pa, const ParticleAttrib<rpos_type>& pb,
            ParticleAttrib<rpos_type>& py)
  {
    for(int i=0; i<pa.size(); ++i)
    {
      py[i][0]=a*pa[i][0]+pb[i][0];
      py[i][1]=a*pa[i][1]+pb[i][1];
      py[i][2]=a*pa[i][2]+pb[i][2];
    }
  }

  static inline
  void copy(const ParticleAttrib<rpos_type>& px, ParticleAttrib<rpos_type>& py)
  {
    py=px;
  }

  static inline
  void copy(const ParticleAttrib<cpos_type>& px, ParticleAttrib<rpos_type>& py)
  {
    for(int i=0; i<px.size(); ++i)
    {
      py[i][0]=px[i][0].real();
      py[i][1]=px[i][1].real();
      py[i][2]=px[i][2].real();
    }
  }
};

///specialization for 2-dimension
template<class T>
struct PAOps<T,2>
{

  typedef T                        real_type;
  typedef complex<T>               complex_type;
  typedef TinyVector<T,2>          rpos_type;
  typedef TinyVector<complex<T>,2> cpos_type;

  static inline
  void scale(T a, const ParticleAttrib<cpos_type>& pa,
             ParticleAttrib<rpos_type>& pb)
  {
    for(int i=0; i<pa.size(); i++)
    {
      pb[i][0]=a*pa[i][0].real();
      pb[i][1]=a*pa[i][1].real();
    }
  }

  static inline
  void scale(T a, const ParticleAttrib<rpos_type>& pa,
             ParticleAttrib<rpos_type>& pb)
  {
    pb=a*pa;
  }


  static inline
  void axpy(T a, const ParticleAttrib<cpos_type>& pa, ParticleAttrib<rpos_type>& pb)
  {
    for(int i=0; i<pa.size(); i++)
    {
      pb[i][0]+=a*pa[i][0].real();
      pb[i][1]+=a*pa[i][1].real();
    }
  }

  static inline
  void axpy(T a, const ParticleAttrib<rpos_type>& pa, ParticleAttrib<rpos_type>& pb)
  {
    for(int i=0; i<pa.size(); i++)
    {
      pb[i][0]+=a*pa[i][0];
      pb[i][1]+=a*pa[i][1];
    }
  }

  static inline
  void axpy(T a, const ParticleAttrib<cpos_type>& pa, const ParticleAttrib<rpos_type>& pb,
            ParticleAttrib<rpos_type>& py)
  {
    for(int i=0; i<pa.size(); i++)
    {
      py[i][0]=a*pa[i][0].real()+pb[i][0];
      py[i][1]=a*pa[i][1].real()+pb[i][1];
    }
  }

  static inline
  void axpy(T a, const ParticleAttrib<rpos_type>& pa, const ParticleAttrib<rpos_type>& pb,
            ParticleAttrib<rpos_type>& py)
  {
    for(int i=0; i<pa.size(); i++)
    {
      py[i][0]=a*pa[i][0]+pb[i][0];
      py[i][1]=a*pa[i][1]+pb[i][1];
    }
  }
};
}
#endif // OHMMS_PARTICLEATTRIB_OPS_H


/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/

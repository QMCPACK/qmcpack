//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    

#ifndef QMCPLUSPLUS_SCLAR_TRAITS_H
#define QMCPLUSPLUS_SCLAR_TRAITS_H
#include <complex>
#include <OhmmsPETE/OhmmsMatrix.h>
#include <OhmmsPETE/Tensor.h>
#include <OhmmsPETE/OhmmsVector.h>
#include <OhmmsPETE/TinyVector.h>
#include "ParticleBase/ParticleAttribOps.h"
namespace qmcplusplus
{

template<class T> struct scalar_traits
{
  enum {DIM=1};
  typedef T real_type;
  typedef T value_type;
  static inline T* get_address(T* a)
  {
    return a;
  }
};

template<typename T>
struct scalar_traits<std::complex<T> >
{
  enum {DIM=2};
  typedef T          real_type;
  typedef std::complex<T> value_type;
  static inline T* get_address(std::complex<T>* a)
  {
    return reinterpret_cast<T*>(a);
  }
};

/** generic conversion from type T1 to type T2 using implicit conversion
*/
template<typename T1, typename T2>
inline void convert(const T1& in, T2& out)
{
  out=static_cast<T2>(in);
}

/** specialization of conversion from complex to real
*/
template<typename T1>
inline void convert(const std::complex<T1>& in, double& out)
{
  out=in.real();
}

template<typename T1>
inline void convert(const std::complex<T1>& in, float& out)
{
  out=in.real();
}

/* specialization of D-dim vectors
 *
 */
template<typename T1, typename T2, unsigned D>
inline void convert(const TinyVector<T1,D>& in, TinyVector<T2,D>& out)
{
  for(int i=0; i<D; ++i)
    convert(in[i],out[i]);
}

/** specialization for 3D */
template<typename T1, typename T2>
inline void convert(const TinyVector<T1,3>& in, TinyVector<T2,3>& out)
{
  convert(in[0],out[0]);
  convert(in[1],out[1]);
  convert(in[2],out[2]);
}

/** specialization for D tensory*/
template<typename T1, typename T2, unsigned D>
inline void convert(const Tensor<T1,D>& in, Tensor<T2,D>& out)
{
  for(int i=0; i<D*D; ++i)
    convert(in[i],out[i]);
}

/** generic function to convert arrays
 * @param in starting address of type T1
 * @param out starting address of type T2
 * @param n size of in/out
 */
template<typename T1, typename T2>
inline void convert(const T1* restrict in, T2* restrict out, std::size_t n)
{
  for(int i=0; i<n; ++i)
    convert(in[i],out[i]);
}

/** specialization for a vector */
template<typename T1, typename T2>
inline void convert(const Vector<T1>& in, Vector<T2>& out)
{
  convert(in.data(),out.data(),in.size());
}

/** specialization for a vector */
template<typename T1, typename T2>
inline void convert(const Matrix<T1>& in, Matrix<T2>& out)
{
  convert(in.data(),out.data(),in.size());
}

/** specialization for a vector */
template<typename T1, typename T2>
inline void convert(const Tensor<T1,3>& in, Tensor<T2,3>& out)
{
  convert(in.data(),out.data(),in.size());
}


// Fix to allow real, imag, conj on scalar and complex types
///real part of a scalar
inline float real(const float& c)
{
  return c;
}
inline double real(const double& c)
{
  return c;
}
inline float real(const std::complex<float>& c)
{
  return c.real();
}
inline double real(const std::complex<double>& c)
{
  return c.real();
}
///imaginary part of a scalar
inline float imag(const float& c)
{
  return 0;
}
inline double imag(const double& c)
{
  return 0;
}
inline float imag(const std::complex<float>& c)
{
  return c.imag();
}
inline double imag(const std::complex<double>& c)
{
  return c.imag();
}
///complex conjugate of a scalar
inline float conj(const float& c)
{
  return c;
}
inline double conj(const double& c)
{
    return c;
}
inline std::complex<float> conj(const std::complex<float>& c)
{
    return std::conj(c);
}
inline std::complex<double> conj(const std::complex<double>& c)
{
    return std::conj(c);
}


///real part of product of scalars
template<typename T>
inline T prod_real(const T& a, const T& b)
{
  return a*b;
}

template<typename T>
inline T prod_real(const T& a, const std::complex<T>& b)
{
  return a*b.real();
}

template<typename T>
inline T prod_real(const std::complex<T>& a, const T& b)
{
  return a.real()*b;
}

template<typename T>
inline T prod_real(const std::complex<T>& a, const std::complex<T>& b)
{
  return a.real()*b.real()-a.imag()*b.imag();
}

///imaginary part of product of scalars
template<typename T>
inline T prod_imag(const T& a, const T& b)
{
  return static_cast<T>(0);
}

template<typename T>
inline T prod_imag(const T& a, const std::complex<T>& b)
{
  return a*b.imag();
}

template<typename T>
inline T prod_imag(const std::complex<T>& a, const T& b)
{
  return a.imag()*b;
}

template<typename T>
inline T prod_imag(const std::complex<T>& a, const std::complex<T>& b)
{
  return a.real()*b.imag()+a.imag()*b.real();
}


///real part of dot product
template<typename T, unsigned D>
inline T dot_real(const TinyVector<T,D>& v1, const TinyVector<T,D>& v2)
{
  return dot(v1,v2);
}

template<typename T, unsigned D>
inline T dot_real(const TinyVector<std::complex<T>,D>& v1, const TinyVector<std::complex<T>,D>& v2)
{
  return OTCDot<T,T,D>::apply(v1,v2);
}

template<typename T, unsigned D>
inline T dot_real(const TinyVector<T,D>& v1, const TinyVector<std::complex<T>,D>& v2)
{
  T res = prod_real(v1[0],v2[0]);
  for(unsigned d=1; d<D; ++d)
    res += prod_real(v1[d],v2[d]);
  return res;
}

template<typename T, unsigned D>
inline T dot_real(const TinyVector<std::complex<T>,D>& v1, const TinyVector<T,D>& v2)
{
  T res = prod_real(v1[0],v2[0]);
  for(unsigned d=1; d<D; ++d)
    res += prod_real(v1[d],v2[d]);
  return res;
}


///imaginary part of dot product
template<typename T, unsigned D>
inline T dot_imag(const TinyVector<T,D>& v1, const TinyVector<T,D>& v2)
{
  return dot(v1,v2);
}

template<typename T, unsigned D>
inline T dot_imag(const TinyVector<std::complex<T>,D>& v1, const TinyVector<std::complex<T>,D>& v2)
{
  T res = prod_imag(v1[0],v2[0]);
  for(unsigned d=1; d<D; ++d)
    res += prod_imag(v1[d],v2[d]);
  return res;
}

template<typename T, unsigned D>
inline T dot_imag(const TinyVector<T,D>& v1, const TinyVector<std::complex<T>,D>& v2)
{
  T res = prod_imag(v1[0],v2[0]);
  for(unsigned d=1; d<D; ++d)
    res += prod_imag(v1[d],v2[d]);
  return res;
}

template<typename T, unsigned D>
inline T dot_imag(const TinyVector<std::complex<T>,D>& v1, const TinyVector<T,D>& v2)
{
  T res = prod_imag(v1[0],v2[0]);
  for(unsigned d=1; d<D; ++d)
    res += prod_imag(v1[d],v2[d]);
  return res;
}




}
#endif

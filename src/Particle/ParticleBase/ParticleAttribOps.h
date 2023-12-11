//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


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
 *  class C is a container class. Default is std::vector<T>
 *
 *  \todo Implement openMP compatible container class or evaluate function.
 *  \todo Implement get/put member functions for MPI-like parallelism
 */
#ifndef OHMMS_PARTICLEATTRIB_OPS_H
#define OHMMS_PARTICLEATTRIB_OPS_H

#include "ParticleUtility.h"

namespace qmcplusplus
{
template<class T1, class T2, unsigned D>
struct OTCDot
{
  using Type_t = typename BinaryReturn<T1, T2, OpMultiply>::Type_t;
  inline static Type_t apply(const TinyVector<std::complex<T1>, D>& lhs, const TinyVector<std::complex<T2>, D>& rhs)
  {
    Type_t res = lhs[0].real() * rhs[0].real() - lhs[0].imag() * rhs[0].imag();
    for (unsigned d = 1; d < D; ++d)
      res += lhs[d].real() * rhs[d].real() - lhs[d].imag() * rhs[d].imag();
    return res;
  }

  inline static std::complex<Type_t> cplx_apply(const TinyVector<std::complex<T1>, D>& lhs,
                                                const TinyVector<std::complex<T2>, D>& rhs)
  {
    std::complex<Type_t> res = lhs[0] * rhs[0];
    for (unsigned d = 1; d < D; ++d)
      res += lhs[d] * rhs[d];
    return res;
  }
};

// Use complex conjugate of the second argument
template<class T1, class T2, unsigned D>
struct OTCDot_CC
{
  using Type_t = typename BinaryReturn<T1, T2, OpMultiply>::Type_t;
  inline static Type_t apply(const TinyVector<std::complex<T1>, D>& lhs, const TinyVector<std::complex<T2>, D>& rhs)
  {
    Type_t res = lhs[0].real() * rhs[0].real() + lhs[0].imag() * rhs[0].imag();
    for (unsigned d = 1; d < D; ++d)
      res += lhs[d].real() * rhs[d].real() + lhs[d].imag() * rhs[d].imag();
    return res;
  }
};

template<class T1, class T2>
struct OTCDot<T1, T2, 3>
{
  using Type_t = typename BinaryReturn<T1, T2, OpMultiply>::Type_t;
  inline static Type_t apply(const TinyVector<std::complex<T1>, 3>& lhs, const TinyVector<std::complex<T2>, 3>& rhs)
  {
    return lhs[0].real() * rhs[0].real() - lhs[0].imag() * rhs[0].imag() + lhs[1].real() * rhs[1].real() -
        lhs[1].imag() * rhs[1].imag() + lhs[2].real() * rhs[2].real() - lhs[2].imag() * rhs[2].imag();
  }

  inline static std::complex<Type_t> cplx_apply(const TinyVector<std::complex<T1>, 3>& lhs,
                                                const TinyVector<std::complex<T2>, 3>& rhs)
  {
    return lhs[0] * rhs[0] + lhs[1] * rhs[1] + lhs[2] * rhs[2];
  }
};

// Use complex conjugate of the second argument
template<class T1, class T2>
struct OTCDot_CC<T1, T2, 3>
{
  using Type_t = typename BinaryReturn<T1, T2, OpMultiply>::Type_t;
  inline static Type_t apply(const TinyVector<std::complex<T1>, 3>& lhs, const TinyVector<std::complex<T2>, 3>& rhs)
  {
    return lhs[0].real() * rhs[0].real() + lhs[0].imag() * rhs[0].imag() + lhs[1].real() * rhs[1].real() +
        lhs[1].imag() * rhs[1].imag() + lhs[2].real() * rhs[2].real() + lhs[2].imag() * rhs[2].imag();
  }
};

template<class T, unsigned D>
void normalize(ParticleAttrib<TinyVector<T, D>>& pa)
{
  T factor = Dot(pa, pa);
  factor   = 1.0 / sqrt(factor);
  pa *= factor;
}

/////////////////////////////////////////////////////////////////
/*\fn template<class T, unsigned D>
 *    T Dot(const ParticleAttrib<TinyVector<T, D> >& pa,
 *      const ParticleAttrib<TinyVector<T, D> >& pb)
 * \return a dot product of an array
 */
/////////////////////////////////////////////////////////////////
template<typename T, unsigned D>
inline T Dot(const ParticleAttrib<TinyVector<T, D>>& pa, const ParticleAttrib<TinyVector<T, D>>& pb)
{
  T sum = 0;
  for (int i = 0; i < pa.size(); i++)
  {
    sum += dot(pa[i], pb[i]);
  }
  return sum;
}

template<typename T, unsigned D>
inline T Dot(const ParticleAttrib<TinyVector<std::complex<T>, D>>& pa,
             const ParticleAttrib<TinyVector<std::complex<T>, D>>& pb)
{
  T sum = 0;
  for (int i = 0; i < pa.size(); i++)
  {
    sum += OTCDot<T, T, D>::apply(pa[i], pb[i]);
  }
  return sum;
}

template<typename T, unsigned D>
inline std::complex<T> CplxDot(const ParticleAttrib<TinyVector<std::complex<T>, D>>& pa,
                               const ParticleAttrib<TinyVector<std::complex<T>, D>>& pb)
{
  std::complex<T> sum(0.0, 0.0);
  for (int i = 0; i < pa.size(); i++)
  {
    sum += OTCDot<T, T, D>::cplx_apply(pa[i], pb[i]);
  }
  return sum;
}

// Use complex conjugate of the second argument
template<unsigned D>
inline double Dot_CC(const ParticleAttrib<TinyVector<std::complex<double>, D>>& pa,
                     const ParticleAttrib<TinyVector<std::complex<double>, D>>& pb)
{
  double sum = 0;
  for (int i = 0; i < pa.size(); i++)
  {
    sum += OTCDot_CC<double, double, D>::apply(pa[i], pb[i]);
  }
  return sum;
}

template<typename T>
inline T Sum(const ParticleAttrib<T>& pa)
{
  T sum = 0;
  for (int i = 0; i < pa.size(); i++)
  {
    sum += pa[i];
  }
  return sum;
}

template<typename T>
inline T Sum(const ParticleAttrib<std::complex<T>>& pa)
{
  T sum = 0;
  for (int i = 0; i < pa.size(); i++)
  {
    sum += pa[i].real();
  }
  return sum;
}

template<typename T>
inline std::complex<T> CplxSum(const ParticleAttrib<std::complex<T>>& pa)
{
  std::complex<T> sum(0.0, 0.0);
  for (int i = 0; i < pa.size(); i++)
  {
    sum += pa[i];
  }
  return sum;
}

template<class T, unsigned D>
inline void Copy(const ParticleAttrib<TinyVector<std::complex<T>, D>>& c, ParticleAttrib<TinyVector<T, D>>& r)
{
  for (int i = 0; i < c.size(); i++)
  {
    //r[i]=real(c[i]);
    for (int j = 0; j < D; j++)
      r[i][j] = c[i][j].real();
  }
}

template<class T, unsigned D>
inline void Copy(const ParticleAttrib<TinyVector<T, D>>& c, ParticleAttrib<TinyVector<T, D>>& r)
{
  r = c;
}


/** generic PAOps
 */
template<class T, unsigned D, class T1 = T>
struct PAOps
{};

///specialization for three-dimension
template<class T, class T1>
struct PAOps<T, 3, T1>
{
  using real_type    = T;
  using complex_type = std::complex<T>;
  using rpos_type    = TinyVector<T, 3>;
  using ipos_type    = TinyVector<T1, 3>;
  using cpos_type    = TinyVector<std::complex<T1>, 3>;

  static inline void scale(T a, const ParticleAttrib<cpos_type>& pa, ParticleAttrib<rpos_type>& pb)
  {
    for (int i = 0; i < pa.size(); i++)
    {
      pb[i][0] = a * pa[i][0].real();
      pb[i][1] = a * pa[i][1].real();
      pb[i][2] = a * pa[i][2].real();
    }
  }

  static inline void scale(T a, const ParticleAttrib<ipos_type>& pa, ParticleAttrib<rpos_type>& pb) { pb = a * pa; }

  static inline void axpy(T a, const ParticleAttrib<cpos_type>& pa, ParticleAttrib<rpos_type>& pb)
  {
    for (int i = 0; i < pa.size(); ++i)
    {
      pb[i][0] += a * pa[i][0].real();
      pb[i][1] += a * pa[i][1].real();
      pb[i][2] += a * pa[i][2].real();
    }
  }

  static inline void axpy(T a, const ParticleAttrib<ipos_type>& pa, ParticleAttrib<rpos_type>& pb)
  {
    for (int i = 0; i < pa.size(); ++i)
    {
      pb[i][0] += a * pa[i][0];
      pb[i][1] += a * pa[i][1];
      pb[i][2] += a * pa[i][2];
    }
  }

  static inline void axpy(T a,
                          const ParticleAttrib<cpos_type>& pa,
                          const ParticleAttrib<ipos_type>& pb,
                          ParticleAttrib<rpos_type>& py)
  {
    for (int i = 0; i < pa.size(); ++i)
    {
      py[i][0] = a * pa[i][0].real() + pb[i][0];
      py[i][1] = a * pa[i][1].real() + pb[i][1];
      py[i][2] = a * pa[i][2].real() + pb[i][2];
    }
  }

  static inline void axpy(T a,
                          const ParticleAttrib<ipos_type>& pa,
                          const ParticleAttrib<ipos_type>& pb,
                          ParticleAttrib<rpos_type>& py)
  {
    for (int i = 0; i < pa.size(); ++i)
    {
      py[i][0] = a * pa[i][0] + pb[i][0];
      py[i][1] = a * pa[i][1] + pb[i][1];
      py[i][2] = a * pa[i][2] + pb[i][2];
    }
  }

  static inline void copy(const ParticleAttrib<ipos_type>& px, ParticleAttrib<rpos_type>& py) { py = px; }

  static inline void copy(const ParticleAttrib<cpos_type>& px, ParticleAttrib<rpos_type>& py)
  {
    for (int i = 0; i < px.size(); ++i)
    {
      py[i][0] = px[i][0].real();
      py[i][1] = px[i][1].real();
      py[i][2] = px[i][2].real();
    }
  }
};

///specialization for 2-dimension
template<class T, class T1>
struct PAOps<T, 2, T1>
{
  using real_type    = T;
  using complex_type = std::complex<T>;
  using rpos_type    = TinyVector<T, 2>;
  using ipos_type    = TinyVector<T1, 2>;
  using cpos_type    = TinyVector<std::complex<T1>, 2>;

  static inline void scale(T a, const ParticleAttrib<cpos_type>& pa, ParticleAttrib<rpos_type>& pb)
  {
    for (int i = 0; i < pa.size(); i++)
    {
      pb[i][0] = a * pa[i][0].real();
      pb[i][1] = a * pa[i][1].real();
    }
  }

  static inline void scale(T a, const ParticleAttrib<ipos_type>& pa, ParticleAttrib<rpos_type>& pb) { pb = a * pa; }


  static inline void axpy(T a, const ParticleAttrib<cpos_type>& pa, ParticleAttrib<rpos_type>& pb)
  {
    for (int i = 0; i < pa.size(); i++)
    {
      pb[i][0] += a * pa[i][0].real();
      pb[i][1] += a * pa[i][1].real();
    }
  }

  static inline void axpy(T a, const ParticleAttrib<ipos_type>& pa, ParticleAttrib<rpos_type>& pb)
  {
    for (int i = 0; i < pa.size(); i++)
    {
      pb[i][0] += a * pa[i][0];
      pb[i][1] += a * pa[i][1];
    }
  }

  static inline void axpy(T a,
                          const ParticleAttrib<cpos_type>& pa,
                          const ParticleAttrib<ipos_type>& pb,
                          ParticleAttrib<rpos_type>& py)
  {
    for (int i = 0; i < pa.size(); i++)
    {
      py[i][0] = a * pa[i][0].real() + pb[i][0];
      py[i][1] = a * pa[i][1].real() + pb[i][1];
    }
  }

  static inline void axpy(T a,
                          const ParticleAttrib<ipos_type>& pa,
                          const ParticleAttrib<ipos_type>& pb,
                          ParticleAttrib<rpos_type>& py)
  {
    for (int i = 0; i < pa.size(); i++)
    {
      py[i][0] = a * pa[i][0] + pb[i][0];
      py[i][1] = a * pa[i][1] + pb[i][1];
    }
  }
};
} // namespace qmcplusplus
#endif // OHMMS_PARTICLEATTRIB_OPS_H

//////////////////////////////////////////////////////////////////
// (c) Copyright 2008-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
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
/**@file fft_help.h
 * @brief common header file for FFT tests
 */
#ifndef QMCPLUSPLUS_FFT_HELP_H
#define QMCPLUSPLUS_FFT_HELP_H
#include <OhmmsPETE/OhmmsMatrix.h>
#include <OhmmsPETE/OhmmsArray.h>
#include <fft/fft.h>

namespace qmcplusplus
{
template<typename T>
inline T abs_diff(T a, T b)
{
  return std::abs(a-b);
}

template<typename T>
inline T abs_diff(const std::complex<T>& a, const std::complex<T>& b)
{
  return std::abs(a.real()-b.real())+std::abs(a.imag()-b.imag());
}

template<typename T>
inline T abs_diff(const std::complex<T>& a, T b)
{
  return std::abs(a.real()-b);
}
template<typename T, class RG>
inline void randomize(T* restrict a, size_t n, RG& rng)
{
  for(int i=0; i<n; i++)
    a[i] = rng();
}

//template<typename T>
//  inline void init_array(Array<T,3>& in, MKL_LONG* strides_in)
//  {
//    long offset=in.size(1)*in.size(2);
//    for (int i = 0; i < in.size(0); i++)
//      init_multiple_2d_columns_z(in.data()+i*offset, in.size(1), in.size(2), strides_in);
//  }

template<typename T>
inline void init_array(Array<T,3>& in)
{
  long strides_in[3];
  typedef typename scalar_traits<T>::real_type real_type;
  real_type sqrt3=std::sqrt(3.0),step;
  for(int i=0; i<in.size(0); ++i)
  {
    for(int j=0; j<in.size(1); ++j)
      for(int k=0; k<in.size(2); ++k)
      {
        step=std::sin(static_cast<real_type>(j+1)*static_cast<real_type>(j+k+1));
        in(i,j,k)=complex<double>(step*sqrt3*0.5,step/sqrt3);
      }
  }
}

template<typename T>
inline void init_array(Matrix<T>& in)
{
  typedef typename scalar_traits<T>::real_type real_type;
  real_type sqrt3=std::sqrt(3.0),step;
  for(int j=0; j<in.rows(); ++j)
    for(int k=0; k<in.cols(); ++k)
    {
      step=std::sin(static_cast<real_type>(k+1)*static_cast<real_type>(in.cols()+1));
      in(j,k)=complex<double>(step*sqrt3*0.5,step/sqrt3);
    }
}

template<typename T>
inline void transpose(const Matrix<T>& in, Matrix<T>& out)
{
  typedef typename scalar_traits<T>::real_type real_type;
  real_type sqrt3=std::sqrt(3.0),step;
  for(int j=0; j<in.rows(); ++j)
    for(int k=0; k<in.cols(); ++k)
      out(k,j)=in(j,k);
}



//template<typename T>
//  inline void print_array(Array<T,3>& in, MKL_LONG* strides_in)
//  {
//    long offset=in.size(1)*in.size(2);
//    printf("\n INPUT  X three 2D columns\n");
//    for (int i = 0; i < in.size(0); i++)
//    {
//      printf("\n Transform number = %d\n", i);
//      print_three_2d_columns_z(in.data()+i*offset, in.size(1), strides_in);
//    }
//  }
//
template<typename T>
inline void print_array(Matrix<T>& in)
{
  for(int j=0; j<in.rows(); ++j)
  {
    for(int k=0; k<in.cols(); ++k)
      printf(" (%8.3f, %8.3f) ", in(j,k).real(), in(j,k).imag());
    printf("\n");
  }
}

template<typename T>
inline void print_array(Array<T,3>& in)
{
  for(int i=0; i<in.size(0); ++i)
  {
    printf("\n Transform number = %d\n", i);
    for(int j=0; j<in.size(1); ++j)
    {
      for(int k=0; k<in.size(2); ++k)
        printf(" (%8.3f, %8.3f) ", in(i,j,k).real(), in(i,j,k).imag());
      printf("\n");
    }
  }
}
//template<typename T, class RG>
//  inline void randomize(std::complex<T>* restrict a, size_t n, RG& rng)
//  {
//    for(int i=0; i<n;i++) a[i] = std::complex<T>(rng(),rng());
//  }

template<typename T>
inline bool check_array(const T* restrict a, const T* restrict b, size_t n
                        , typename scalar_traits<T>::real_type norm)
{
  typedef typename scalar_traits<T>::real_type real_type;
  bool no_problem=true;
  real_type eps=1000.*numeric_limits<real_type>::epsilon();
  for(int i=0; i<n; ++i)
  {
    if(abs_diff(a[i]*norm,b[i])>eps)
    {
      cerr << "WRONG " << a[i]*norm << " " << b[i] <<  " " << endl;
      no_problem=false;
    }
  }
  return no_problem;
}

#if defined(HAVE_LIBFFTW)
template<typename T1>
inline fftw_complex* mangle(T1* in)
{
  return reinterpret_cast<fftw_complex*>(in);
}
#endif

}//end-of-namespace
#endif

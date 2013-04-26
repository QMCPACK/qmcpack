//////////////////////////////////////////////////////////////////
// (c) Copyright 2010-  by Jeongnim Kim
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
#ifndef QMCPLUSPLUS_TRANSPOSE_OPS
#define QMCPLUSPLUS_TRANSPOSE_OPS

#include <fft/fft.h>
#if defined(HAVE_MKL)
#include <mkl_trans.h>
#endif
#if defined(HAVE_ESSL)
#include <essl.h>
#endif
#include <OhmmsPETE/OhmmsMatrix.h>

#if !defined(__IBMCPP__)
#define transpose_xy transpose_xy_
#define transpose_yx transpose_yx_
#define transpose_1 transpose_1_
#define trans_b trans_b_
#endif

extern "C"
{
  void transpose_1(const int& nx,const int& ny,const int& first_x, const int& last_x
                   ,const std::complex<double>* input, std::complex<double>* output);

  void transpose_xy(const int* nx, const int* ny, const int* howmany
                    , const int* first, const int* last
                    , const std::complex<double>* restrict, std::complex<double>* restrict);

  void transpose_yx(const int* nx, const int* ny, const int* howmany
                    , const int* first, const int* last
                    , const std::complex<double>* restrict, std::complex<double>* restrict);

  void trans_b(const int &nx_th, const int  &ny, const std::complex<double>* restrict in
               ,const  int &ny_th,const  int &nx, std::complex<double>* restrict out
               , int *i_start, int *j_start);
}



namespace qmcplusplus
{
enum {DUMMY_TRANSPOSER=0
                       , SWAP_TRANSPOSER
      , MKL_TRANSPOSER
      , ESSL_TRANSPOSER
     };

/** C version of transpose
 */
template<typename T>
inline void transpose_b(int nx_th, int  ny, const T* restrict in
                        , int ny_th, int nx, T* restrict out, int i_start, int j_start)
{
  for(int i=0,ii=i_start; i<ny_th; ++i,++ii)
    for(int j=0,jj=j_start; j<nx_th; ++j,++jj)
      out[i*nx+jj]=in[j*ny+ii];
}

template<typename MT>
inline void transpose_block(const vector<MT*>& in, MT& out, int ip)
{
  int npy=out.rows();
  for(int jp=0; jp<in.size(); ++jp)
  {
    int npx=in[jp]->rows();
    int i_offset=ip*npy;
    int j_offset=jp*npx;
    trans_b(npx, in[jp]->cols(), in[jp]->data()
            , npy, out.cols(), out.data()
            , &i_offset, &j_offset);
    //transpose_b(npx, in[jp]->cols(), in[jp]->data()
    //    , npy, out.cols(), out.data()
    //    , ip*npy, jp*npx);
  }
}

template<typename MT>
inline void transpose_block(const vector<MT*>& in, MT& out
                            , int ip, int np, int ivar, int howmany)
{
  int npy=out.rows();
  for(int jp=0; jp<np; ++jp)
  {
    int jj=jp*howmany+ivar;
    int npx=in[jj]->rows();
    int i_offset=ip*npy;
    int j_offset=jp*npx;
    trans_b(npx, in[jj]->cols(), in[jj]->data()
            , npy, out.cols(), out.data()
            , &i_offset, &j_offset);
    //transpose_b(npx, in[jp]->cols(), in[jp]->data()
    //    , npy, out.cols(), out.data()
    //    , ip*npy, jp*npx);
  }
}

/** dummy declaration to be specialized */
template <typename T, unsigned ENG> struct Transpose2D {};

/** basic implementation of 2D transposes
 */
template <typename T>
struct Transpose2D<T,DUMMY_TRANSPOSER>
{
  typedef T value_type;
  static inline void apply(const Matrix<T>& in, Matrix<T>& out)
  {
    int n=out.rows();
    int m=out.cols();
    T* restrict optr=out.data();
    for(int i=0; i<n; ++i)
    {
      const T* restrict sptr=in.data()+i;
      for(int j=0; j<m; ++j,sptr+=n,++optr)
        *optr=*sptr;
    }
  }
};

/** swap target and source pointer access */
template <typename T>
struct Transpose2D<T,SWAP_TRANSPOSER>
{
  typedef T value_type;
  static inline void apply(const Matrix<T>& in, Matrix<T>& out)
  {
    int n=out.rows();
    int m=out.cols();
    const T* restrict iptr=in.data();
    for(int i=0; i<m; ++i)
    {
      T* restrict optr=out.data()+i;
      for(int j=0; j<n; ++j,optr+=m)
        *optr=*iptr++;
    }
  }
};

#if defined(HAVE_MKL)
/** specialization for double with mkl transpose function */
template<>
struct Transpose2D<double,MKL_TRANSPOSER>
{
  typedef double value_type;
  static inline void apply(Matrix<value_type>& in, Matrix<value_type>& out)
  {
    mkl_domatcopy('R','T',in.rows(),in.cols(),1.0,in.data(),in.cols(),out.data(),out.cols());
  }

  static inline void apply(value_type* restrict iptr, int rows, int cols, value_type* restrict optr)
  {
    mkl_domatcopy('R','T',rows,cols,1.0,iptr,cols,optr,rows);
  }

};

template<>
struct Transpose2D<std::complex<double>,MKL_TRANSPOSER>
{
  typedef std::complex<double> value_type;
  static inline void apply(Matrix<value_type>& in, Matrix<value_type>& out)
  {
    _MKL_Complex16 one;
    one.real=1.0;
    one.imag=0.0;
    mkl_zomatcopy('R','T',in.rows(),in.cols(),one,mkl_mangle(in.data())
                  ,in.cols(),mkl_mangle(out.data()),out.cols());
  }

  static inline void apply(value_type* restrict iptr, int rows, int cols, value_type* restrict optr)
  {
    _MKL_Complex16 one;
    one.real=1.0;
    one.imag=0.0;
    mkl_zomatcopy('R','T',rows,cols,one,mkl_mangle(iptr),cols,mkl_mangle(optr),rows);
  }

};

/** specialization for float with mkl transpose function */
template<>
struct Transpose2D<float,MKL_TRANSPOSER>
{
  typedef float value_type;
  static inline void apply(Matrix<value_type>& in, Matrix<value_type>& out)
  {
    mkl_somatcopy('R','T',in.rows(),in.cols(),1.0,in.data(),in.cols(),out.data(),out.cols());
  }
};

template<>
struct Transpose2D<std::complex<float>,MKL_TRANSPOSER>
{
  typedef std::complex<float> value_type;
  static inline void apply(Matrix<value_type>& in, Matrix<value_type>& out)
  {
    _MKL_Complex8 one;
    one.real=1.0;
    one.imag=0.0;
    mkl_comatcopy('R','T',in.rows(),in.cols(),one,mkl_mangle(in.data())
                  ,in.cols(),mkl_mangle(out.data()),out.cols());
  }

  static inline void apply(value_type* restrict iptr, int rows, int cols, value_type* restrict optr)
  {
    _MKL_Complex8 one;
    one.real=1.0;
    one.imag=0.0;
    mkl_comatcopy('R','T',rows,cols,one,mkl_mangle(iptr),cols,mkl_mangle(optr),rows);
  }
};

#endif
#if defined(HAVE_ESSL)
/** specialization for double with essl transpose function */
template<>
struct Transpose2D<double,ESSL_TRANSPOSER>
{
  typedef double value_type;
  static inline void apply(const Matrix<value_type>& in, Matrix<value_type>& out)
  {
    dgetmo(in.data(),in.cols(),in.cols(),in.rows(),out.data(),out.cols());
  }
};

template<>
struct Transpose2D<std::complex<double>,ESSL_TRANSPOSER>
{
  typedef std::complex<double> value_type;
  static inline void apply(const Matrix<value_type>& in, Matrix<value_type>& out)
  {
    zgetmo(in.data(),in.cols(),in.cols(),in.rows(),out.data(),out.cols());
  }
};

template<>
struct Transpose2D<float,ESSL_TRANSPOSER>
{
  typedef float value_type;
  static inline void apply(const Matrix<value_type>& in, Matrix<value_type>& out)
  {
    sgetmo(in.data(),in.cols(),in.cols(),in.rows(),out.data(),out.cols());
  }
};

template<>
struct Transpose2D<std::complex<float>,ESSL_TRANSPOSER>
{
  typedef std::complex<float> value_type;
  static inline void apply(const Matrix<value_type>& in, Matrix<value_type>& out)
  {
    cgetmo(in.data(),in.cols(),in.cols(),in.rows(),out.data(),out.cols());
  }
};
#endif

//template<typename MT>
//inline void transpose(const vector<MT*>& in, MT& out, int ip)
//{
//  const int np=in.size();
//  for(int i=0, ii=ip*out.rows(); i<out.rows(); ++i,++ii)
//  {
//    for(int jp=0; jp<np; ++jp)
//      zcopy(in[jp]->rows()
//          ,in[jp]->data()+ii,in[jp]->cols()
//          ,out[i]+jp*in[jp]->rows(),1
//          );
//  }
//}
}
#endif
/***************************************************************************
 * $RCSfile$   $Author: jmcminis $
 * $Revision: 4631 $   $Date: 2010-02-18 14:08:30 -0600 (Thu, 18 Feb 2010) $
 * $Id: transpose.h 4631 2010-02-18 20:08:30Z jmcminis $
 ***************************************************************************/

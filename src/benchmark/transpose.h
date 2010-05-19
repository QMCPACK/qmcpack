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
#include <OhmmsPETE/OhmmsMatrix.h>

#if defined(HAVE_MKL)
#include <mkl_trans.h>
#endif
#if defined(HAVE_ESSL)
#include <essl.h>
#endif

namespace qmcplusplus
{
  enum {DUMMY_TRANSPOSER=0
    , SWAP_TRANSPOSER
    , MKL_TRANSPOSER
    , ESSL_TRANSPOSER
  };

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
        for(int j=0; j<m; ++j,sptr+=n,++optr) *optr=*sptr;
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
        for(int j=0; j<n; ++j,optr+=m) *optr=*iptr++;
      }
    }
  };

#if defined(HAVE_MKL)
  /** specialization for double with mkl transpose function */
  template<>
    struct Transpose2D<double,MKL_TRANSPOSER>
    {
      typedef double value_type;
      static inline void apply(const Matrix<value_type>& in, Matrix<value_type>& out)
      {
        mkl_domatcopy('R','T',in.rows(),in.cols(),1.0,in.data(),in.cols(),out.data(),out.cols());
      }

      typedef double value_type;
      static inline void apply(const value_type* restrict iptr, int rows, int cols, value_type* restrict optr)
      {
        mkl_domatcopy('R','T',rows,cols,1.0,iptr,cols,optr,rows);
      }

    };

  template<>
    struct Transpose2D<std::complex<double>,MKL_TRANSPOSER>
    {
      typedef std::complex<double> value_type;
      static inline void apply(const Matrix<value_type>& in, Matrix<value_type>& out)
      {
        mkl_zomatcopy('R','T',in.rows(),in.cols(),1.0,in.data(),in.cols(),out.data(),out.cols());
      }

      typedef double value_type;
      static inline void apply(const value_type* restrict iptr, int rows, int cols, value_type* restrict optr)
      {
        mkl_zomatcopy('R','T',rows,cols,1.0,iptr,cols,optr,rows);
      }

    };

  /** specialization for float with mkl transpose function */
  template<>
    struct Transpose2D<float,MKL_TRANSPOSER>
    {
      typedef float value_type;
      static inline void apply(const Matrix<value_type>& in, Matrix<value_type>& out)
      {
        mkl_somatcopy('R','T',in.rows(),in.cols(),1.0,in.data(),in.cols(),out.data(),out.cols());
      }
    };

  template<>
    struct Transpose2D<std::complex<float>,MKL_TRANSPOSER>
    {
      typedef std::complex<float> value_type;
      static inline void apply(const Matrix<value_type>& in, Matrix<value_type>& out)
      {
        mkl_zomatcopy('R','T',in.rows(),in.cols(),1.0,in.data(),in.cols(),out.data(),out.cols());
      }

      typedef float value_type;
      static inline void apply(const value_type* restrict iptr, int rows, int cols, value_type* restrict optr)
      {
        mkl_zomatcopy('R','T',rows,cols,1.0,iptr,cols,optr,rows);
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

}
#endif
/***************************************************************************
 * $RCSfile$   $Author: jmcminis $
 * $Revision: 4631 $   $Date: 2010-02-18 14:08:30 -0600 (Thu, 18 Feb 2010) $
 * $Id: transpose.h 4631 2010-02-18 20:08:30Z jmcminis $ 
 ***************************************************************************/

//////////////////////////////////////////////////////////////////
// (c) Copyright 2010-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
/** @file mkl_engine.h
 * @brief definition of fft_engine_base<fftw_real,FFTMKL_ENG> 
 */
#ifndef QMCPLUSPLUS_MKL_ENGINE_BASE_H
#define QMCPLUSPLUS_MKL_ENGINE_BASE_H
#include <mkl_dfti.h>

namespace qmcplusplus
{
  inline _MKL_Complex16* mkl_mangle(std::complex<double>* in)
  {
    return reinterpret_cast<_MKL_Complex16*>(in);
  }

  inline _MKL_Complex8* mkl_mangle(std::complex<float>* in)
  {
    return reinterpret_cast<_MKL_Complex8*>(in);
  }

  inline DFTI_CONFIG_VALUE dfti_get_precision(double a)
  {
    return DFTI_DOUBLE;
  }
  inline DFTI_CONFIG_VALUE dfti_get_precision(float a)
  {
    return DFTI_SINGLE;
  }


  /** fft_engine_base<T,FFTMKL_ENG> to be specialized for mkl
   *
   * DFT_DESCRIPTOR can be created only once.
   */
  template<typename T>
    class fft_engine_base<T,FFTMKL_ENG>
    {
      public:
        typedef T real_type;
        typedef std::complex<T> complex_type;
        typedef DFTI_DESCRIPTOR_HANDLE fft_plan_type;
        fft_plan_type my_handle;

        ///default constructor
        fft_engine_base():my_handle(0) { } 

        ///virtual destructor to clean up the plans
        virtual ~fft_engine_base() 
        { 
          if(my_handle) DftiFreeDescriptor(&my_handle); 
        }

        /** plan for inplace/outplace, complex-to-complex  transform
         */
        void create_plan(int dims, int howmany, complex_type* in, complex_type* out , int dir, unsigned uflag)
        {
          create_c2c_desc(dims,howmany,in!=out);
        }

        /** plan for inplace or outplace, real-to-complex  transform
         */
        void create_plan(int dims, int howmany, real_type* in, real_type* out , int dir, unsigned uflag)
        {
          create_r2c_desc(dims,howmany,in!=out);
        }

        /** plan for outplace, real-to-complex */
        void create_plan(int dims, int howmany, real_type* in, complex_type* out, int idir, unsigned uflag)
        {
          create_r2c_desc(dims,howmany,true);
        }

        /** plan for outplace, complex-to-real */
        void create_plan(int dims, int howmany, complex_type* in, real_type* out, int idir, unsigned uflag)
        {
          create_r2c_desc(dims,howmany,true);
        }

        inline void execute_fft(complex_type* inout)
        {
          DftiComputeForward(my_handle,mkl_mangle(inout));
        }

        inline void execute_ifft(complex_type* inout)
        {
          DftiComputeBackward(my_handle,mkl_mangle(inout));
        }

        inline void execute_fft(complex_type* in, complex_type* out)
        {
          DftiComputeForward(my_handle,mkl_mangle(in),mkl_mangle(out));
        }
        inline void execute_ifft(complex_type* in, complex_type* out)
        {
          DftiComputeBackward(my_handle,mkl_mangle(in),mkl_mangle(out));
        }

        inline void execute_fft(real_type* in, complex_type* out)
        {
          DftiComputeForward(my_handle,in,mkl_mangle(out));
        }
        inline void execute_ifft(complex_type* in, real_type* out)
        {
          DftiComputeBackward(my_handle,mkl_mangle(in),out);
        }

      private:
        /** create DFFI_DESCRIPTOR for complex-to-complex transformations */
        void create_c2c_desc(int dims, int howmany,bool outplace)
        {
          if(my_handle) return;
          //if(my_handle) DftiFreeDescriptor(&my_handle); 
          DFTI_CONFIG_VALUE my_precision=dfti_get_precision(real_type());
          DftiCreateDescriptor(&my_handle,my_precision,DFTI_COMPLEX,1,dims);
          if(outplace) DftiSetValue(my_handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
          DftiSetValue(my_handle,DFTI_NUMBER_OF_TRANSFORMS, howmany);
          DftiSetValue(my_handle,DFTI_OUTPUT_DISTANCE, dims);
          DftiCommitDescriptor(my_handle); 
        }

        /** create DFFI_DESCRIPTOR for real-to-complex/complex-to-real transformations */
        void create_r2c_desc(int dims, int howmany, bool outplace)
        {
          if(my_handle) return;
          //if(my_handle) DftiFreeDescriptor(&my_handle); 
          DFTI_CONFIG_VALUE my_precision=dfti_get_precision(real_type());
          DftiCreateDescriptor(&my_handle,my_precision,DFTI_REAL,1,dims);
          if(outplace) DftiSetValue(my_handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
          if(howmany>1) DftiSetValue(my_handle,DFTI_NUMBER_OF_TRANSFORMS, howmany);
          DftiSetValue(my_handle,DFTI_OUTPUT_DISTANCE, dims+2);
          DftiCommitDescriptor(my_handle); 
        }
    };

}
#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 3738 $   $Date: 2009-04-07 02:08:20 -0500 (Tue, 07 Apr 2009) $
 * $Id: fftw_engine.h 3738 2009-04-07 07:08:20Z jnkim $ 
 ***************************************************************************/

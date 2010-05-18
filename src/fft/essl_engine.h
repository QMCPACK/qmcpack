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
#ifndef QMCPLUSPLUS_ESSL_FFT_ENGINE_H
#define QMCPLUSPLUS_ESSL_FFT_ENGINE_H
#include <vector>
extern "C"
{
#include <essl.h>
}

namespace qmcplusplus
{
  const int ESSL_FFT_FORWARD=1;
  const int ESSL_FFT_BACKWARD=-1;

  /** convenient functions based on the function overloads */
  namespace essl
  {
    /** class to define a fft_plan similar to mkl descriptor
     * \tparam T real datatype
     *
     * essl uses double/float for the aux storage
     */
    template<typename T>
      struct fft_plan
      {
        int fft_size;
        int num_ffts;
        int in_displ;
        int out_displ;
        int in_stride;
        int out_stride;

        std::vector<T> vaux1;
        std::vector<T> vaux2;
        fft_plan():fft_size(0),num_ffts(0),in_displ(0),out_displ(0){}
        ~fft_plan() { }
        void resize(int* desc, int eflag)
        {
          fft_size=desc[FFT_LENGTH];
          num_ffts=desc[FFT_NUMBER_OF_TRANSFORMS];
          int naux1=0,naux2=0;
          if(eflag==ESSL_FFT_FORWARD)
          {
            in_displ=desc[FFT_IN_DISTANCE];
            out_displ=desc[FFT_OUT_DISTANCE];
            in_stride=desc[FFT_IN_STRIDE];
            out_stride=desc[FFT_OUT_STRIDE];
          }
          else
          {
            in_displ=desc[FFT_OUT_DISTANCE];
            in_stride=desc[FFT_OUT_STRIDE];
            out_displ=desc[FFT_IN_DISTANCE];
            out_stride=desc[FFT_IN_STRIDE];
          }
          if(desc[FFT_COMPLEX])
          {
            //in_displ=out_displ=fft_size;
            naux1=(fft_size>4096)?static_cast<int>(20000.+1.64*fft_size)+1:25000;
            naux2=naux1+(2*fft_size+256)*std::min(64,num_ffts);
          }
          else
          {
            naux1=(fft_size>2048)?static_cast<int>(20000.+2.28*fft_size)+1:20000;
            naux2=std::max(naux1,(2*fft_size+256)*std::min(64,num_ffts));
          }
          vaux1.resize(naux1);
          vaux2.resize(naux2);
        }

        inline int size1() const { return vaux1.size();}
        inline T* data1() { return &vaux1[0];}
        inline int size2() const{ return vaux2.size();}
        inline T* data2() { return &vaux2[0];}
      };

    /** initialize the aux storage for complex-to-complex 1d fft
     * @param aplan fft_plan
     * @param in address of the input data
     * @param out starting address of the output data
     * @param idir ESSL_FFT_FORWARD or ESSL_FFT_BACKWARD
     */
    inline void create_plan(fft_plan<double>& aplan, std::complex<double>* in, std::complex<double>* out, int eflag)
    {
      dcft(1 
          , in, aplan.in_stride, aplan.in_displ 
          , out, aplan.out_stride, aplan.out_displ
          , aplan.fft_size, aplan.num_ffts
          , eflag, 1.0
          , aplan.data1(), aplan.size1(), aplan.data2(), aplan.size2());
    }

    inline void fft(fft_plan<double>& aplan, std::complex<double>* in, std::complex<double>* out)
    {
      dcft(0 
          , in, aplan.in_stride, aplan.in_displ 
          , out, aplan.out_stride, aplan.out_displ
          , aplan.fft_size, aplan.num_ffts
          , ESSL_FFT_FORWARD, 1.0
          , aplan.data1(), aplan.size1(), aplan.data2(), aplan.size2());
    }

    inline void ifft(fft_plan<double>& aplan, std::complex<double>* in, std::complex<double>* out)
    {
      dcft(0 
          , in, aplan.in_stride, aplan.in_displ 
          , out, aplan.out_stride, aplan.out_displ
          , aplan.fft_size, aplan.num_ffts
          , ESSL_FFT_BACKWARD, 1.0
          , aplan.data1(), aplan.size1(), aplan.data2(), aplan.size2());
    }

    /** initialize the aux storage for real-to-complex 1d fft
     * @param aplan fft_plan
     * @param in address of the input data
     * @param out starting address of the output data
     */
    inline void create_plan(fft_plan<double>& aplan, double* in, std::complex<double>* out)
    {
      drcft(1 , in, aplan.in_displ , out, aplan.out_displ
          , aplan.fft_size, aplan.num_ffts
          , ESSL_FFT_FORWARD, 1.0
          , aplan.data1(), aplan.size1(), aplan.data2(), aplan.size2());
    }

    /** initialize the aux storage for omplex-to-real 1d fft
     * @param aplan fft_plan
     * @param in address of the input data
     * @param out starting address of the output data
     */
    inline void create_plan(fft_plan<double>& aplan, std::complex<double>* in, double* out)
    {
      dcrft(1 , in, aplan.in_displ , out, aplan.out_displ
          , aplan.fft_size, aplan.num_ffts
          , ESSL_FFT_BACKWARD, 1.0
          , aplan.data1(), aplan.size1(), aplan.data2(), aplan.size2());
    }

    inline void fft(fft_plan<double>& aplan, double* in, std::complex<double>* out)
    {
      drcft(0 , in, aplan.in_displ , out, aplan.out_displ
          , aplan.fft_size, aplan.num_ffts, ESSL_FFT_FORWARD, 1.0
          , aplan.data1(), aplan.size1(), aplan.data2(), aplan.size2());
    }
    inline void ifft(fft_plan<double>& aplan, std::complex<double>* in, double* out)
    {
      dcrft(0, in, aplan.in_displ, out, aplan.out_displ
          , aplan.fft_size, aplan.num_ffts, ESSL_FFT_BACKWARD, 1.0
          , aplan.data1(), aplan.size1(), aplan.data2(), aplan.size2());
    }
  }//end-of-essl namespace

  /** fft_engine_base<T,FFTESSL_ENG> to be specialized for ESSL
   *
   * DFT_DESCRIPTOR can be created only once.
   */
  template<typename T>
    class fft_engine_base<T,FFTESSL_ENG>
    {
      public:
        typedef T real_type;
        typedef std::complex<T> complex_type;
        typedef essl::fft_plan<T> fft_plan_type;
        fft_plan_type forward_plan;
        fft_plan_type backward_plan;

        ///default constructor
        fft_engine_base() { } 

        ///virtual destructor to clean up the plans
        virtual ~fft_engine_base() { }

        /** plan for inplace/outplace, complex-to-complex  transform
         */
        void create_plan(int* desc, real_type* in, real_type* out 
            , int dir, unsigned uflag)
        {
          if(in==out)
            std::cerr << "Creat plans for r2c inplace. Incorrect!!! "<< endl;
          else
            std::cerr << "Creat plans for r2c out place. Incorrect!!! "<< endl;
        }

        void create_plan(int* desc, complex_type* in, complex_type* out , int dir, unsigned uflag)
        {
          if(dir<0)
          {
            forward_plan.resize(desc,ESSL_FFT_FORWARD);
            backward_plan.resize(desc,ESSL_FFT_BACKWARD);
            essl::create_plan(forward_plan,in,out,ESSL_FFT_FORWARD);
            essl::create_plan(backward_plan,out,in,ESSL_FFT_BACKWARD);
          }
          else
          {
            if(dir==FFTW_FORWARD)
            {
              forward_plan.resize(desc,ESSL_FFT_FORWARD);
              essl::create_plan(forward_plan,in,out,ESSL_FFT_FORWARD);
            }
            else
            {
              backward_plan.resize(desc,ESSL_FFT_BACKWARD);
              essl::create_plan(backward_plan,out,in,ESSL_FFT_BACKWARD);
            }
          }
        }

        /** plan for outplace, real-to-complex */
        void create_plan(int* desc, real_type* in, complex_type* out, int idir, unsigned uflag)
        {
          forward_plan.resize(desc,ESSL_FFT_FORWARD);
          essl::create_plan(forward_plan,in,out);
        }

        /** plan for outplace, complex-to-real */
        void create_plan(int* desc, complex_type* in, real_type* out, int idir, unsigned uflag)
        {
          backward_plan.resize(desc,ESSL_FFT_BACKWARD);
          essl::create_plan(backward_plan,in,out);
        }

        inline void execute_fft(real_type* inout)
        {
          std::cerr << "fft r2c inplace is not working!!! "<< endl;
        }

        inline void execute_ifft(real_type* inout)
        {
          std::cerr << "ifft r2c inplace is not working!!! "<< endl;
        }

        inline void execute_fft(complex_type* inout)
        {
          essl::fft(forward_plan,inout,inout);
        }

        inline void execute_ifft(complex_type* inout)
        {
          essl::ifft(backward_plan,inout,inout);
        }

        inline void execute_fft(complex_type* in, complex_type* out)
        {
          essl::fft(forward_plan,in,out);
        }

        inline void execute_ifft(complex_type* in, complex_type* out)
        {
          essl::ifft(backward_plan,in,out);
        }

        inline void execute_fft(real_type* in, complex_type* out)
        {
          essl::fft(forward_plan,in,out);
        }

        inline void execute_ifft(complex_type* in, real_type* out)
        {
          essl::ifft(backward_plan,in,out);
        }
    };

}
#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 3738 $   $Date: 2009-04-07 02:08:20 -0500 (Tue, 07 Apr 2009) $
 * $Id: fftw_engine.h 3738 2009-04-07 07:08:20Z jnkim $ 
 ***************************************************************************/

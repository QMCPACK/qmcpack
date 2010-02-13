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
/** @file fftw_engine.h
 * @brief definition of fft_engine_base<fftw_real,FFTW_ENG> 
 */
#ifndef QMCPLUSPLUS_FFTW_ENGINE_BASE_H
#define QMCPLUSPLUS_FFTW_ENGINE_BASE_H

#include <fftw3.h>
typedef double fftw_real;

namespace qmcplusplus
{
  /** fftw_base to be specialized with the precision
   *
   * fftw_base class defines the interfaces to fftw libraries with matching argument types
   */
  template<>
    class fft_engine_base<fftw_real,FFTW_ENG>
    {
      template<typename T1>
        inline fftw_complex* mangle(T1* in)
        {
          return reinterpret_cast<fftw_complex*>(in);
        }

      public:
        typedef std::complex<fftw_real> complex_type;
        typedef fftw_plan fft_plan_type;
        fft_plan_type forward_plan;
        fft_plan_type backward_plan;

        ///default constructor
        fft_engine_base() :forward_plan(0), backward_plan(0) { }

        ///virtual destructor to clean up the plans
        virtual ~fft_engine_base() { clear(); }

        ///clean up the allocated data
        void clear()
        {
          if(forward_plan) {fftw_destroy_plan(forward_plan); forward_plan=0;}
          if(backward_plan) {fftw_destroy_plan(backward_plan); backward_plan=0;}
        }

        /** plan for outplace, complex-to-complex  transform
         */
        void create_plan(int dims, int howmany, complex_type* in, complex_type* out , int dir, unsigned uflag)
        {
          if(dir==FFTW_FORWARD)
          { 
            if(forward_plan) fftw_destroy_plan(forward_plan);
            forward_plan=create_plan_c2c(dims,howmany,mangle(in),mangle(out),dir,uflag);
          }
          else
          {
            if(backward_plan) fftw_destroy_plan(backward_plan);
            backward_plan=create_plan_c2c(dims,howmany,mangle(in),mangle(out),dir,uflag);
          }
        }

        /** plan for inplace or outplace, real-to-complex  transform
         */
        void create_plan(int dims, int howmany, fftw_real* in, fftw_real* out , int dir, unsigned uflag)
        {
          if(dir==FFTW_FORWARD)
          { 
            if(forward_plan) fftw_destroy_plan(forward_plan);
            forward_plan=create_plan_r2c(dims,howmany,in,mangle(out),uflag);
          }
          else
          {
            if(backward_plan) fftw_destroy_plan(backward_plan);
            backward_plan=create_plan_c2r(dims,howmany,mangle(in),out,uflag);
          }
        }

        /** plan for outplace, real-to-complex */
        void create_plan(int dims, int howmany, fftw_real* in, complex_type* out, int idir, unsigned uflag)
        {
          if(forward_plan) fftw_destroy_plan(forward_plan);
          forward_plan= create_plan_r2c(dims,howmany,in,mangle(out),uflag);
        }

        /** plan for outplace, complex-to-real */
        void create_plan(int dims, int howmany, complex_type* in, fftw_real* out, int idir, unsigned uflag)
        {
          if(backward_plan) fftw_destroy_plan(backward_plan);
          backward_plan=create_plan_c2r(dims,howmany,mangle(in),out,uflag);
        }

        inline void execute_fft(complex_type* inout)
        {
          fftw_execute_dft(forward_plan,mangle(inout),mangle(inout));
        }

        inline void execute_fft(complex_type* in, complex_type* out)
        {
          fftw_execute_dft(forward_plan,mangle(in),mangle(out));
        }

        inline void execute_ifft(complex_type* inout)
        {
          fftw_execute_dft(backward_plan,mangle(inout),mangle(inout));
        }

        inline void execute_ifft(complex_type* in, complex_type* out)
        {
          fftw_execute_dft(backward_plan,reinterpret_cast<fftw_complex*>(in),reinterpret_cast<fftw_complex*>(out));
        }

        inline void execute_fft(fftw_real* in, complex_type* out)
        {
          fftw_execute_dft_r2c(forward_plan,in,reinterpret_cast<fftw_complex*>(out));
        }

        inline void execute_ifft(complex_type* in, fftw_real* out)
        {
          fftw_execute_dft_c2r(backward_plan,reinterpret_cast<fftw_complex*>(in),out);
        }

      private:
        fft_plan_type create_plan_c2c(int dims, int howmany, fftw_complex* in, fftw_complex* out, int idir, unsigned uflag)
        {
          if(howmany>1)
          {
            fftw_iodim data_dims;
            fftw_iodim howmany_dims;
            data_dims.n=dims; data_dims.is=1; data_dims.os=1;
            howmany_dims.n=howmany; howmany_dims.is=dims; howmany_dims.os=dims;
            return fftw_plan_guru_dft(1,&data_dims, 1,&howmany_dims,in,out,idir,uflag);
          }
          else
          {
            return fftw_plan_dft_1d(dims,in,out,idir,uflag);
          }
        }
        fft_plan_type create_plan_r2c(int dims, int howmany, fftw_real* in, fftw_complex* out, unsigned uflag)
        {
          if(howmany>1)
          {
            fftw_iodim data_dims;
            fftw_iodim howmany_dims;
            data_dims.n=dims; data_dims.is=1; data_dims.os=1;
            howmany_dims.n=howmany; howmany_dims.is=dims+2; howmany_dims.os=dims/2+1;
            return fftw_plan_guru_dft_r2c(1,&data_dims, 1,&howmany_dims, in, out,uflag);
          }
          else
            return fftw_plan_dft_r2c_1d(dims,in,out,uflag);
        }

        /** create complex-to-real backward plan 
         * @param dims fft size
         * @param howmany number of concurrent ffts
         * @param in complex input data
         * @param out real output data
         * @param uflag fftw plan (measure ...)
         */
        fft_plan_type create_plan_c2r(int dims, int howmany, fftw_complex* in, fftw_real* out, unsigned uflag)
        {
          if(howmany>1)
          {
            fftw_iodim data_dims;
            fftw_iodim howmany_dims;
            data_dims.n=dims; data_dims.is=1; data_dims.os=1;
            howmany_dims.n=howmany; howmany_dims.is=dims/2+1; howmany_dims.os=dims+2;
            return fftw_plan_guru_dft_c2r(1,&data_dims, 1,&howmany_dims,in,out,uflag);
          }
          else
            return fftw_plan_dft_c2r_1d(dims,in,out,uflag);
        }

    };

}
#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 3738 $   $Date: 2009-04-07 02:08:20 -0500 (Tue, 07 Apr 2009) $
 * $Id: fftw_engine.h 3738 2009-04-07 07:08:20Z jnkim $ 
 ***************************************************************************/

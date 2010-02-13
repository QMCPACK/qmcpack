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
/** @file fft1d.h
 * @brief A master header file to define 1d fft interfaces
 */
#ifndef QMCPLUSPLUS_FFT_1D_H
#define QMCPLUSPLUS_FFT_1D_H

namespace qmcplusplus
{
  /** generic traits for fft1d engine property used by fft1d_engine
   */
  template<typename T1, typename T2> struct fft1d_property{};

  /** specialization of fft1d_property for complex-to-complex type
   */
  template<typename T>
    struct fft1d_property<T,T>
  {
    int num_ffts;
    int fft_size;
    int f_offset;
    int b_offset;
    inline fft1d_property(int dims=0):num_ffts(0),fft_size(dims),f_offset(dims),b_offset(dims){}
    inline void resize(int dims, int howmany=1)
    {
      fft_size=dims;
      b_offset=f_offset=dims;
      num_ffts=howmany;
    }
  };

  /** specialization of fft1d_property for real-to-complex type
   */
  template<typename T>
    struct fft1d_property<T,std::complex<T> >
  {
    int num_ffts;
    int fft_size;
    int f_offset;
    int b_offset;
    inline fft1d_property(int dims=0):num_ffts(0),fft_size(dims),f_offset(dims+2),b_offset(dims/2+1){}
    inline void resize(int dims, int howmany=1)
    {
      num_ffts=howmany;
      fft_size=dims;
      f_offset=dims+2;
      b_offset=dims/2+1;
    }
  };

  /** declaration of FFT_DFT 
   *
   * @tparam PREC precision of the data
   * @tparam T1  datatype for the real space
   * @tparam T2  datatype for the spectral space
   * @tparam ENG enumration of the FFT library that is employed
   */
  template<typename T1, typename T2, unsigned ENG>
    struct fft1d_engine
  {
    typedef T1 space_type;
    typedef T2 spectral_type;
    fft1d_property<T1,T2> my_property;
    fft_engine_base<typename scalar_traits<T1>::real_type,ENG> my_engine;

    /** constructor with a fft dimension
     * @param dims size of 1dfft
     * @param m number of ffts
     */
    explicit fft1d_engine(int dims=0):my_property(dims){ }

    /** return the size of 1d fft */
    inline int size() const { return my_property.fft_size;}
    inline int howmany() const { return my_property.num_ffts;}
    inline int offset(int idir)
    {
      return (idir==FFTW_FORWARD)?my_property.f_offset:my_property.b_offset;
    }

    /** create InPlace plan 
     * @param dims fft dimension
     * @param m number of ffts
     * @param in starting address
     * @param idir if idir<0, create both the forward and backward plans
     * @param uflag fftw-specific plan
     */
    void create(int dims, int m, T2* in, int idir=-1,unsigned uflag=FFTW_MEASURE)
    {
      my_property.resize(dims,m);
      if(idir<0)
      {
        my_engine.create_plan(dims,m,in,in,FFTW_FORWARD,uflag);
        my_engine.create_plan(dims,m,in,in,FFTW_BACKWARD,uflag);
      }
      else
        my_engine.create_plan(dims,m,in,in,idir,uflag);
    }

    /** create OutPlace plan 
     */
    void create(int dims, int m, T1* in, T2* out, int idir=-1, unsigned uflag=FFTW_MEASURE)
    {
      my_property.resize(dims,m);
#pragma omp critical
      {
        if(idir<0)
        {
          my_engine.create_plan(dims,m,in,out,FFTW_FORWARD,uflag);
          my_engine.create_plan(dims,m,out,in,FFTW_BACKWARD,uflag);
        }
        else
        {
          if(idir==FFTW_FORWARD) 
            my_engine.create_plan(dims,m,in,out,FFTW_FORWARD,uflag);
          else
            my_engine.create_plan(dims,m,out,in,FFTW_BACKWARD,uflag);
        }
      }
    }

    inline void fft_forward(T1* in) { my_engine.execute_fft(in); }
    inline void fft_backward(T2* in) { my_engine.execute_ifft(in); }
    inline void fft_forward(T1* in, T2* out) { my_engine.execute_fft(in,out); }
    inline void fft_backward(T2* in, T1* out) { my_engine.execute_ifft(in,out); }
  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 3738 $   $Date: 2009-04-07 02:08:20 -0500 (Tue, 07 Apr 2009) $
 * $Id: fftw_engine.h 3738 2009-04-07 07:08:20Z jnkim $ 
 ***************************************************************************/

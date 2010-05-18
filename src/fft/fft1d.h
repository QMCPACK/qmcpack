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
    typedef typename scalar_traits<T1>::real_type real_type;
    int mydesc[FFT_MAX];
    fft_engine_base<real_type,ENG> my_engine;

    /** default constructor
     */
    inline fft1d_engine()
    { 
      mydesc[FFT_COMPLEX]=is_complex2complex<T1,T2>::value;
    }

    inline double gflops() const
    {
      double nops=(std::log(static_cast<double>(mydesc[FFT_LENGTH]))/std::log(2.))*5.
        *mydesc[FFT_LENGTH]*mydesc[FFT_NUMBER_OF_TRANSFORMS]*1e-9;
      return mydesc[FFT_COMPLEX]? nops:nops*0.5;
    }

    /** operator to set the FFT properties */
    int& operator()(int i) { return mydesc[i];}
    /** operator to get the FFT properties */
    int operator()(int i) const { return mydesc[i];}

    /** return the size of 1d fft */
    inline int size() const { return mydesc[FFT_LENGTH];}
    inline int howmany() const { return mydesc[FFT_NUMBER_OF_TRANSFORMS];}
    inline int offset(int idir)
    {
      return (idir==FFTW_FORWARD)?mydesc[FFT_IN_DISTANCE]:mydesc[FFT_OUT_DISTANCE];
    }

    inline void transpose(int idir)
    {
      if(idir == FFTW_FORWARD)
      {
        mydesc[FFT_IN_DISTANCE]=1;
        mydesc[FFT_IN_STRIDE]=mydesc[FFT_NUMBER_OF_TRANSFORMS];
      }
      else
      {
        mydesc[FFT_OUT_DISTANCE]=1;
        mydesc[FFT_OUT_STRIDE]=mydesc[FFT_NUMBER_OF_TRANSFORMS];
      }
    }

    /** create InPlace plan 
     * @param dims fft dimension
     * @param m number of ffts
     * @param in starting address
     * @param idir if idir<0, create both the forward and backward plans
     * @param uflag fftw-specific plan
     */
    template<typename T>
    void create(int dims, int m, T* in, int idir=-1,unsigned uflag=FFTW_MEASURE)
    {
      set_defaults(dims,m);
      create(in,idir,uflag);
    }

    /** create OutPlace plan 
     */
    void create(int dims, int m, T1* in, T2* out, int idir=-1, unsigned uflag=FFTW_MEASURE)
    {
      set_defaults(dims,m);
      create(in,out,idir,uflag);
    }

    template<typename T>
      void create(T* in, int idir=-1,unsigned uflag=FFTW_MEASURE)
      {
        mydesc[FFT_INPLACE]=1;
        if(idir<0)
        {
          my_engine.create_plan(mydesc,in,in,FFTW_FORWARD,uflag);
          my_engine.create_plan(mydesc,in,in,FFTW_BACKWARD,uflag);
        }
        else
          my_engine.create_plan(mydesc,in,in,idir,uflag);
      }

    /** create OutPlace plan 
     */
    void create(T1* in, T2* out, int idir=-1, unsigned uflag=FFTW_MEASURE)
    {
      mydesc[FFT_INPLACE]=0;
      if(idir<0)
      {
        my_engine.create_plan(mydesc,in,out,FFTW_FORWARD,uflag);
        my_engine.create_plan(mydesc,out,in,FFTW_BACKWARD,uflag);
      }
      else
      {
        if(idir==FFTW_FORWARD) 
          my_engine.create_plan(mydesc,in,out,FFTW_FORWARD,uflag);
        else
          my_engine.create_plan(mydesc,out,in,FFTW_BACKWARD,uflag);
      }
    }

    template<typename T>
    inline void fft_forward(T* in) { my_engine.execute_fft(in); }
    template<typename T>
    inline void fft_backward(T* in) { my_engine.execute_ifft(in); }

    inline void fft_forward(T1* in, T2* out) { my_engine.execute_fft(in,out); }
    inline void fft_backward(T2* in, T1* out) { my_engine.execute_ifft(in,out); }

    inline void set_defaults(int dims, int m)
    {
      mydesc[FFT_COMPLEX]=is_complex2complex<T1,T2>::value;
      mydesc[FFT_LENGTH]=dims;
      mydesc[FFT_NUMBER_OF_TRANSFORMS]=m;
      if(mydesc[FFT_COMPLEX])
        mydesc[FFT_IN_DISTANCE]=mydesc[FFT_OUT_DISTANCE]=dims;
      else
      {
        mydesc[FFT_IN_DISTANCE]=dims+2;
        mydesc[FFT_OUT_DISTANCE]=dims/2+1;
      }
      mydesc[FFT_IN_STRIDE]=mydesc[FFT_OUT_STRIDE]=1;
    }
  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 3738 $   $Date: 2009-04-07 02:08:20 -0500 (Tue, 07 Apr 2009) $
 * $Id: fftw_engine.h 3738 2009-04-07 07:08:20Z jnkim $ 
 ***************************************************************************/

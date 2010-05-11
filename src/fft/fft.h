/** @file fft.h
 * @brief A master header file to define fft interfaces
 */
#ifndef QMCPLUSPLUS_FFT_MASTER_H
#define QMCPLUSPLUS_FFT_MASTER_H
#include <config.h>
#include <type_traits/scalar_traits.h>

namespace qmcplusplus
{
  ///enumeration to determine what FFT  engine is used
  enum {FFTW_ENG=0  /*< FFTW  */
    , FFTESSL_ENG=1 /*< FFT ESSL */
    , FFTMKL_ENG=2  /* FFT MKL */
    , MAX_FFT_ENG /* size of this enum */
  };

#if !defined(HAVE_LIBFFTW)
  const int FFTW_ESTIMATE=0;
  const int FFTW_MEASURE=0;
  const int FFTW_PATIENT=0;
#endif

  enum { FFT_INPLACE=0
    , FFT_COMPLEX=1
    , FFT_NUMBER_OF_TRANSFORMS=2
    , FFT_LENGTH=3
    , FFT_IN_DISTANCE=4
    , FFT_OUT_DISTANCE=5
    , FFT_IN_STRIDE=6
    , FFT_OUT_STRIDE=7
    , FFT_MAX
  };

  template<typename T1, typename T2>
    struct is_complex2complex
    {
      static const int value=1;
    };

  template<typename T>
    struct is_complex2complex<T,std::complex<T> >
    {
      static const int value=0;
    };


  /** dummy declaration of fft_engine_base */
  template<typename T, unsigned ENG> struct fft_engine_base{ };
}

#if defined(HAVE_LIBFFTW)
#include <fft/fftw_engine.h>
#endif
#if defined(HAVE_ESSL)
#include <fft/essl_engine.h>
#endif
#if defined(HAVE_MKL)
#include <fft/mkl_engine.h>
#endif
#include <fft/fft1d.h>
#endif

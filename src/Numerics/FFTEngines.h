//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



#ifndef FFTENGINES_H
#define FFTENGINES_H
#include <complex>
#include <fftw3.h>

namespace qmcplusplus
{

// dummy class to be specialized
template<unsigned dimensions, typename precision> class FFTWEngine { };

// partial specialization for std::complex<double>
template<unsigned dimensions>
class FFTWEngine<dimensions, std::complex<double> >
{
private:
  fftw_plan forwardPlan;
  fftw_plan backwardPlan;
  inline fftw_complex* mangle(std::complex<double>* p)
  {
    return (fftw_complex*)p;
  }
  FFTWEngine(const FFTWEngine&);
public:
  FFTWEngine() { }
  ~FFTWEngine()
  {
    if (forwardPlan)
    {
      fftw_destroy_plan(forwardPlan);
      forwardPlan = 0;
    }
    if (backwardPlan)
    {
      fftw_destroy_plan(backwardPlan);
      backwardPlan = 0;
    }
  }
  void initialize(const int* DimSizes, std::complex<double>* DataArray)
  {
    forwardPlan = fftw_plan_dft(dimensions, DimSizes, mangle(DataArray),
                                mangle(DataArray), FFTW_FORWARD, FFTW_PATIENT);
    backwardPlan = fftw_plan_dft(dimensions, DimSizes, mangle(DataArray),
                                 mangle(DataArray), FFTW_BACKWARD, FFTW_PATIENT);
  }
  inline void transformForward(std::complex<double>* DataArray)
  {
    fftw_execute_dft(forwardPlan, mangle(DataArray), mangle(DataArray));
  }
  inline void transformBackward(std::complex<double>* DataArray)
  {
    fftw_execute_dft(backwardPlan, mangle(DataArray), mangle(DataArray));
  }
};

// partial specialization for std::complex<long double>
template<unsigned dimensions>
class FFTWEngine<dimensions, std::complex<long double> >
{
private:
  fftwl_plan forwardPlan;
  fftwl_plan backwardPlan;
  inline fftwl_complex* mangle(std::complex<long double>* p)
  {
    return (fftwl_complex*)p;
  }
  FFTWEngine(const FFTWEngine&);
public:
  FFTWEngine() { }
  ~FFTWEngine()
  {
    if (forwardPlan)
    {
      fftwl_destroy_plan(forwardPlan);
      forwardPlan = 0;
    }
    if (backwardPlan)
    {
      fftwl_destroy_plan(backwardPlan);
      backwardPlan = 0;
    }
  }
  void initialize(const int* DimSizes, std::complex<long double>* DataArray)
  {
    forwardPlan = fftwl_plan_dft(dimensions, DimSizes, mangle(DataArray),
                                 mangle(DataArray), FFTW_FORWARD, FFTW_PATIENT);
    backwardPlan = fftwl_plan_dft(dimensions, DimSizes, mangle(DataArray),
                                  mangle(DataArray), FFTW_BACKWARD, FFTW_PATIENT);
  }
  inline void transformForward(std::complex<long double>* DataArray)
  {
    fftwl_execute_dft(forwardPlan, mangle(DataArray), mangle(DataArray));
  }
  inline void transformBackward(std::complex<long double>* DataArray)
  {
    fftwl_execute_dft(backwardPlan, mangle(DataArray), mangle(DataArray));
  }
};

// partial specialization for std::complex<float>
template<unsigned dimensions>
class FFTWEngine<dimensions, std::complex<float> >
{
private:
  fftwf_plan forwardPlan;
  fftwf_plan backwardPlan;
  inline fftwf_complex* mangle(std::complex<float>* p)
  {
    return (fftwf_complex*)p;
  }
  FFTWEngine(const FFTWEngine&);
public:
  FFTWEngine() { }
  ~FFTWEngine()
  {
    if (forwardPlan)
    {
      fftwf_destroy_plan(forwardPlan);
      forwardPlan = 0;
    }
    if (backwardPlan)
    {
      fftwf_destroy_plan(backwardPlan);
      backwardPlan = 0;
    }
  }
  void initialize(const int* DimSizes, std::complex<float>* DataArray)
  {
    forwardPlan = fftwf_plan_dft(dimensions, DimSizes, mangle(DataArray),
                                 mangle(DataArray), FFTW_FORWARD, FFTW_PATIENT);
    backwardPlan = fftwf_plan_dft(dimensions, DimSizes, mangle(DataArray),
                                  mangle(DataArray), FFTW_BACKWARD, FFTW_PATIENT);
  }
  inline void transformForward(std::complex<float>* DataArray)
  {
    fftwf_execute_dft(forwardPlan, mangle(DataArray), mangle(DataArray));
  }
  inline void transformBackward(std::complex<float>* DataArray)
  {
    fftwf_execute_dft(backwardPlan, mangle(DataArray), mangle(DataArray));
  }
};

}
#endif

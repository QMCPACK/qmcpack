//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



/**@file multidet.cpp
 * @brief Test codes for multidets
 */
#include <Configuration.h>
#include "Utilities/OhmmsInfo.h"
#include "Utilities/RandomGenerator.h"
#include "Utilities/Timer.h"
#include "Message/Communicate.h"
#include <fft/fft.h>

template<typename T>
inline T abs_diff(T a, T b)
{
  return std::abs(a-b);
}

template<typename T>
inline T abs_diff(std::complex<T>& a, std::complex<T>& b)
{
  return std::abs(a.real()-b.real())+std::abs(a.imag()-b.imag());
}

template<typename T>
inline T abs_diff(std::complex<T>& a, T b)
{
  return std::abs(a.real()-b);
}

int main(int argc, char** argv)
{
  using namespace qmcplusplus;
  OHMMS::Controller->initialize(argc,argv);
  Communicate* mycomm=OHMMS::Controller;
  OhmmsInfo Welcome(argc,argv,mycomm->rank());
  qmcplusplus::Random.init(0,1,11);
  int howmany=2;
  int niters=10;
  int fft_dim=16;
  typedef double real_type;
  typedef std::complex<double> complex_type;
  //const unsigned int fft_id=FFTW_ENG;
  const unsigned int fft_id=FFTMKL_ENG;
  typedef fft1d_engine<real_type,complex_type,fft_id> fft1d_engine_t;
  //typedef fft1d_engine<complex_type,complex_type,fft_id> fft1d_engine_t;
  typedef fft1d_engine_t::space_type space_type;
  typedef fft1d_engine_t::spectral_type spectral_type;
  fft1d_engine_t myfft;
  myfft.set_defaults(fft_dim,howmany);
  Matrix<space_type> in(myfft(FFT_NUMBER_OF_TRANSFORMS),myfft(FFT_IN_DISTANCE));
  Matrix<spectral_type> out(myfft(FFT_NUMBER_OF_TRANSFORMS),myfft(FFT_OUT_DISTANCE));
  std::vector<space_type> sample(fft_dim);
  real_type phase=2*M_PI/static_cast<real_type>(fft_dim);
  real_type norm=1.0/static_cast<real_type>(fft_dim);
  for(int i=0; i<sample.size(); ++i)
    sample[i]=0.5*sin(phase*i)+0.1*sin(phase*2*i)+0.3*sin(phase*3*i)+0.5*sin(phase*4.1*i);;
  //myfft.create(in.data(),out.data());
  myfft.create(in.data());
  for(int i=0; i<howmany; ++i)
    copy(sample.begin(),sample.end(),in[i]);
  myfft.fft_forward(in.data());
  myfft.fft_backward(in.data());
  //myfft.fft_forward(in.data(),out.data());
  //myfft.fft_backward(out.data(),in.data());
  real_type eps=10*numeric_limits<real_type>::epsilon();
  in *= norm;
  for(int i=0; i<howmany; ++i)
  {
    std::cout << "State index = " << i << std::endl;
    for(int k=0; k<fft_dim; ++k)
    {
      //cout << target(i,k) << " " << sample[k] << std::endl;
      if(abs_diff(in(i,k),sample[k])>eps)
        std::cerr << "WRONG " << in(i,k) << sample[k] <<  " " << std::endl;
    }
  }
  OHMMS::Controller->finalize();
  return 0;
}


//////////////////////////////////////////////////////////////////
// (c) Copyright 2008-  by Jeongnim Kim
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
/**@file multidet.cpp
 * @brief Test codes for multidets
 */
#include <Configuration.h>
#include "Utilities/OhmmsInfo.h"
#include "Utilities/RandomGenerator.h"
#include "Utilities/Timer.h"
#include "Message/Communicate.h"
//#define SINGLE_PREC
#include <fft/fft.h>

namespace qmcplusplus
{
template<typename T1, typename T2, unsigned fft_eng>
struct fft1d_test
{
  typedef fft1d_engine<T1,T2,fft_eng> fft1d_engine_t;
  typedef typename fft1d_engine_t::space_type space_type;
  typedef typename fft1d_engine_t::spectral_type spectral_type;
  typedef typename fft1d_engine_t::real_type real_type;

  fft1d_engine_t myfft;
  Matrix<space_type> in;
  Matrix<spectral_type> out;

  fft1d_test(int fft_dim, int howmany)
  {
    reset(fft_dim,howmany);
  }

  void reset(int fft_dim, int howmany)
  {
    myfft.set_defaults(fft_dim,howmany);
    in.resize(myfft(FFT_NUMBER_OF_TRANSFORMS),myfft(FFT_IN_DISTANCE));
    out.resize(myfft(FFT_NUMBER_OF_TRANSFORMS),myfft(FFT_OUT_DISTANCE));
    std::vector<space_type> sample(fft_dim);
    real_type phase=2*M_PI/static_cast<real_type>(fft_dim);
    for(int i=0; i<sample.size(); ++i)
      sample[i]=0.5*std::sin(phase*i)+0.1*std::sin(phase*2*i)+0.3*std::sin(phase*3*i)+0.5*std::sin(phase*4.1*i);;
    for(int i=0; i<howmany; ++i)
      std::copy(sample.begin(),sample.end(),in[i]);
    myfft.create(in.data(),out.data());
  }

  double gflops() const
  {
    return myfft.gflops();
  }

  inline void fft()
  {
    myfft.fft_forward(in.data(),out.data());
  }

  inline void ifft()
  {
    myfft.fft_backward(out.data(),in.data());
  }
};
}


inline void print_help(const string& msg)
{
  printf("%s -d fft_dim -m numer_of_fft -i iterations -p [d|s] -t [r2c|c2c] -e [fftw|mkl|essl] \n",msg.c_str());
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
  //accepted: r2c  or c2c
  string fft_type("r2c");
  string fft_eng("fftw");
  bool inplace=true;
  bool single_precision=false;
  int opt;
  while((opt = getopt(argc, argv, "hp:d:m:i:e:t:")) != -1)
  {
    switch(opt)
    {
    case 'h':
      print_help("Help Message");
      return 1;
    case 'd':
      fft_dim=atoi(optarg);
      break;
    case 'm':
      howmany=atoi(optarg);
      break;
    case 'i':
      niters=atoi(optarg);
      break;
    case 'p':
      if(optarg[0]=='s')
        single_precision=true;
      break;
    case 'e':
      fft_eng=optarg;
      break;
    case 't':
      fft_type=optarg;
      break;
    default:
      print_help("Unknown options");
      return 1;
    }
  }
  double fft_time=0.0;
  double ifft_time=0.0;
  double gflops=0.0;
  #pragma omp parallel  reduction(+:fft_time, ifft_time, gflops)
  {
    //fft1d_test<double,complex<double>,FFTW_ENG> fft_tester(fft_dim,howmany);
    //fft1d_test<float,complex<float>,FFTESSL_ENG> fft_tester(fft_dim,howmany);
    fft1d_test<complex<float>,complex<float>,FFTESSL_ENG> fft_tester(fft_dim,howmany);
    double t_f=0.0;
    double t_b=0.0;
    Timer clock;
    for(int i=0; i<niters; ++i)
    {
      clock.restart();
      fft_tester.fft();
      t_f += clock.elapsed();
      clock.restart();
      fft_tester.ifft();
      t_b += clock.elapsed();
    }
    fft_time+= t_f;
    ifft_time+= t_b;
    gflops += fft_tester.gflops();
  }
  gflops /= omp_get_max_threads();
  fft_time /= omp_get_max_threads()*niters;
  ifft_time /= omp_get_max_threads()*niters;
  cout << "gflops " << gflops << endl;
  cout << fft_time <<" " << ifft_time << " " << gflops/fft_time << " " << gflops/ifft_time << endl;
  OHMMS::Controller->finalize();
  return 0;
}

/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1770 $   $Date: 2007-02-17 17:45:38 -0600 (Sat, 17 Feb 2007) $
 * $Id: OrbitalBase.h 1770 2007-02-17 23:45:38Z jnkim $
 ***************************************************************************/

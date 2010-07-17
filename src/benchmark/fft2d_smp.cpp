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
/**@file fft2d.cpp
 * @brief Test codes for 2D FFT with OpenMP
 */
#include <Configuration.h>
#include "Utilities/OhmmsInfo.h"
#include "Utilities/RandomGenerator.h"
#include "Utilities/Timer.h"
#include "OhmmsPETE/OhmmsArray.h"
#include "Message/Communicate.h"
#include <benchmark/fft_help.h>
#include <benchmark/transpose.h>
#include <Numerics/OhmmsBlas.h>

#if defined(HAVE_ESSL)
#define TEST_FFT_ENG FFTESSL_ENG
#define TEST_TRANSPOSER ESSL_TRANSPOSER
#elif defined(HAVE_MKL)
#define TEST_FFT_ENG FFTMKL_ENG
#define TEST_TRANSPOSER MKL_TRANSPOSER
#else
#error "Only tested with ESSL and MKL library "
#endif

inline void print_help(const string& msg)
{
  printf("%s -d fft_dim -m numer_of_fft -i iterations -p [d|s] -t [r2c|c2c] -e [fftw|mkl|essl] \n",msg.c_str());
}

template<typename MT>
inline void transpose(const vector<MT*>& in, MT& out, int ip)
{
  const int np=in.size();
  for(int i=0, ii=ip*out.rows(); i<out.rows(); ++i,++ii)
  {
    for(int jp=0; jp<np; ++jp)
      BLAS::copy(in[jp]->rows() 
          ,in[jp]->data()+ii,in[jp]->cols() 
          ,out[i]+jp*in[jp]->rows(),1
          ); 
  }
}

int main(int argc, char** argv)
{

  using namespace qmcplusplus;
  OHMMS::Controller->initialize(argc,argv);
  Communicate* mycomm=OHMMS::Controller;
  OhmmsInfo Welcome(argc,argv,mycomm->rank());
  qmcplusplus::Random.init(0,1,11);

  int howmany=4;
  int niters=10;
  int nx=6;
  int ny=6;

  //accepted: r2c  or c2c
  string fft_type("r2c");
  string fft_eng("fftw");
  bool inplace=true;
  bool single_precision=false;
  bool debug=false;
  int opt;
  while((opt = getopt(argc, argv, "hdp:x:y:m:i:e:t:s:")) != -1) {
    switch(opt) {
      case 'h':
        print_help("Help Message");
        return 1;
      case 'd':
        debug=true;
        break;
      case 'x':
        nx=atoi(optarg);
        break;
      case 'y':
        ny=atoi(optarg);
        break;
      case 'm':
        howmany=atoi(optarg);
        break;
      case 'i':
        niters=atoi(optarg);
        break;
      case 'p':
        if(optarg[0]=='s') single_precision=true;
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


  typedef double real_type;
  typedef complex<real_type> complex_type;
  typedef fft1d_engine<complex_type, complex_type, TEST_FFT_ENG> fft1d_engine_t;
  typedef Matrix<complex_type> matrix_type;

  int np=omp_get_max_threads();
  vector<matrix_type*> in(omp_get_max_threads())
    ,in_t(omp_get_max_threads()) ,in_copy(omp_get_max_threads());
  vector<fft1d_engine_t*> fft_xy(np);
  vector<fft1d_engine_t*> fft_yx(np);

  int nx_thread=nx/np;
  int ny_thread=ny/np;

  Timer myclock;
#pragma omp parallel
  {
    int ip=omp_get_thread_num();
    in[ip]=new matrix_type(nx_thread,ny);
    in_copy[ip]=new matrix_type(nx_thread,ny);
    in_t[ip]=new matrix_type(ny_thread,nx);

    //for(int i=0,ii=ip*nx_thread; i<nx_thread; ++i,++ii)
    //  for(int jp=0; jp<np; ++jp)
    //    for(int jj=jp*ny_thread; jj<(jp+1)*ny_thread; ++jj)
    //    (*in_copy[ip])(i,jj)=complex_type(ii,jj);

    init_array(*in_copy[ip]);
    *in[ip]=*in_copy[ip];

    fft_xy[ip]=new fft1d_engine_t;
    fft_xy[ip]->set_defaults(ny,nx_thread);
    fft_xy[ip]->create(in[ip]->data());

    fft_yx[ip]=new fft1d_engine_t;
    fft_yx[ip]->set_defaults(nx,ny_thread);
    fft_yx[ip]->create(in_t[ip]->data());
  }

//DEBUG transpose
//  for(int ip=0; ip<np; ++ip)
//  {
//    cout << *in_copy[ip];
//  }
//
//#pragma omp parallel 
//  {
//    int ip=omp_get_thread_num();
//    transpose(in,*in_t[ip],ip);
//  }
//
//  cout << endl;
//  for(int ip=0; ip<np; ++ip)
//  {
//    cout << *in_t[ip];
//  }
//
//#pragma omp parallel 
//  {
//    int ip=omp_get_thread_num();
//    transpose(in_t,*in[ip],ip);
//  }
//
//  cout << endl;
//  for(int ip=0; ip<np; ++ip)
//  {
//    cout << *in[ip];
//  }
//
  double dt_f=0.0, dt_b=0.0;
#pragma omp parallel  reduction(+:dt_f,dt_b)
  {
    int ip=omp_get_thread_num();
    //init_array(*in_copy[ip]);
    //*in[ip]=*in_copy[ip];

    fft1d_engine_t& myfft_xy(*fft_xy[ip]);
    fft1d_engine_t& myfft_yx(*fft_yx[ip]);

    Timer clock;
    double dt_f_th=0.0, dt_b_th=0.0;
    for(int iter=0; iter<niters; ++iter)
    {
      clock.restart();
      myfft_xy.fft_forward(in[ip]->data());
#pragma omp barrier
      transpose(in,*in_t[ip],ip);
      myfft_yx.fft_forward(in_t[ip]->data());
      dt_f_th+=clock.elapsed();

      clock.restart();
      myfft_yx.fft_backward(in_t[ip]->data());
#pragma omp barrier
      transpose(in_t,*in[ip],ip);
      myfft_xy.fft_backward(in[ip]->data());
      dt_b_th+=clock.elapsed();

      dt_f += dt_f_th;
      dt_b += dt_b_th;
      if(debug && !iter)
        if(check_array(in[ip]->data(),in_copy[ip]->data(),in_copy[ip]->size(),1.0/static_cast<double>(nx*ny))) 
          cout << "We are good with 1D FFT+t(1D FFT)" << endl;
    }
  }

  double factor=1.0/static_cast<double>(niters*np);
  cout << dt_f*factor << " " << dt_b*factor << endl;
  OHMMS::Controller->finalize();
  return 0;
}

/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1770 $   $Date: 2007-02-17 17:45:38 -0600 (Sat, 17 Feb 2007) $
 * $Id: fft2d.cpp 1770 2007-02-17 23:45:38Z jnkim $ 
 ***************************************************************************/

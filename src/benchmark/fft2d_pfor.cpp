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
//#include <Numerics/OhmmsBlas.h>

#if defined(HAVE_ESSL)
#define TEST_FFT_ENG FFTESSL_ENG
#define TEST_TRANSPOSER ESSL_TRANSPOSER
#elif defined(HAVE_MKL)
#define TEST_FFT_ENG FFTMKL_ENG
#define TEST_TRANSPOSER MKL_TRANSPOSER
#else
#define TEST_FFT_ENG FFTW_ENG
#define TEST_TRANSPOSER DUMMY_TRANSPOSER
#endif

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
  int howmany=3;
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
  while((opt = getopt(argc, argv, "hdp:x:y:m:i:e:t:s:")) != -1)
  {
    switch(opt)
    {
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
  typedef double real_type;
  typedef complex<real_type> complex_type;
  typedef fft1d_engine<complex_type, complex_type, TEST_FFT_ENG> fft1d_engine_t;
  typedef Matrix<complex_type> matrix_type;
  int np=omp_get_max_threads();
  vector<matrix_type*> in(np*howmany),in_t(np*howmany) ,in_copy(np*howmany);
  vector<fft1d_engine_t*> fft_xy(np);
  vector<fft1d_engine_t*> fft_yx(np);
  int nx_thread=nx/np;
  int ny_thread=ny/np;
  #pragma omp parallel
  {
    int ip=omp_get_thread_num();
    for(int k=0,kk=ip*howmany; k<howmany; ++k,++kk)
    {
      in[kk]=new matrix_type(nx_thread,ny);
      in_copy[kk]=new matrix_type(nx_thread,ny);
      in_t[kk]=new matrix_type(ny_thread,nx);
    }
    //for(int i=0,ii=ip*nx_thread; i<nx_thread; ++i,++ii)
    //  for(int jp=0; jp<np; ++jp)
    //    for(int jj=jp*ny_thread; jj<(jp+1)*ny_thread; ++jj)
    //    (*in_copy[ip])(i,jj)=complex_type(ii,jj);
    int first=ip*howmany;
    fft_xy[ip]=new fft1d_engine_t;
    fft_xy[ip]->set_defaults(ny,nx_thread);
    fft_yx[ip]=new fft1d_engine_t;
    fft_yx[ip]->set_defaults(nx,ny_thread);
    #pragma omp critical
    {
      fft_xy[ip]->create(in[first]->data());
      fft_yx[ip]->create(in_t[first]->data());
    }
    for(int k=0,kk=ip*howmany; k<howmany; ++k,++kk)
    {
      init_array(*in_copy[kk]);
      *in[kk]=*in_copy[kk];
    }
  }
  double dt_f=0.0, dt_b=0.0, dt_trans=0.0;
  Timer clock_big, clock;
  clock_big.restart();
  for(int iter=0; iter<niters; ++iter)
  {
    clock.restart();
    #pragma omp parallel for
    for(int ip=0; ip<np; ++ip)
    {
      for(int k=0,kk=ip*howmany; k<howmany; ++k,++kk)
        fft_xy[ip]->fft_forward(in[kk]->data());
    }
    #pragma omp parallel for
    for(int ip=0; ip<np; ++ip)
    {
      for(int k=0,kk=ip*howmany; k<howmany; ++k,++kk)
      {
        transpose_block(in,*in_t[kk],ip,np,k,howmany);
        fft_yx[ip]->fft_forward(in_t[kk]->data());
      }
    }
    dt_f+=clock.elapsed();
    clock.restart();
    #pragma omp parallel for
    for(int ip=0; ip<np; ++ip)
    {
      for(int k=0,kk=ip*howmany; k<howmany; ++k,++kk)
        fft_yx[ip]->fft_backward(in_t[kk]->data());
    }
    #pragma omp parallel for
    for(int ip=0; ip<np; ++ip)
    {
      for(int k=0,kk=ip*howmany; k<howmany; ++k,++kk)
      {
        transpose_block(in_t,*in[kk],ip,np,k,howmany);
        fft_xy[ip]->fft_backward(in[kk]->data());
      }
    }
    dt_b+=clock.elapsed();
    if(debug && !iter)
      for(int ip=0; ip<np*howmany; ++ip)
        if(check_array(in[ip]->data(),in_copy[ip]->data(),in_copy[ip]->size(),1.0/static_cast<double>(nx*ny)))
          cout << "We are good with 1D FFT+t(1D FFT) "<< ip << endl;
  }
  double dt_t=clock_big.elapsed();
  double factor=1.0/static_cast<double>(niters);
  printf("tag nX  nY  M  OMP for back tot  \n");
  printf("fft2d %d %d %d %d %12.4e %12.4e %12.4e\n"
         , nx, ny, howmany, omp_get_max_threads()
         , dt_f*factor, dt_b*factor,dt_t*factor);
  OHMMS::Controller->finalize();
  return 0;
}

/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1770 $   $Date: 2007-02-17 17:45:38 -0600 (Sat, 17 Feb 2007) $
 * $Id: fft2d.cpp 1770 2007-02-17 23:45:38Z jnkim $
 ***************************************************************************/

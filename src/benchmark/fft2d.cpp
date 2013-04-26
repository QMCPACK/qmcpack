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
  while((opt = getopt(argc, argv, "hdp:x:y:m:i:e:t:")) != -1)
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
  typedef Matrix<complex_type> matrix_type;
  vector<matrix_type*> in,in_x,in_t;
  matrix_type incopy(nx,ny);
  for(int i=0; i<howmany; ++i)
  {
    in.push_back(new matrix_type(nx,ny));
    in_t.push_back(new matrix_type(ny,nx));
  }
  init_array(incopy);
  for(int i=0; i<howmany; ++i)
    *(in[i])=incopy;
  if(debug)
    niters=1; //reset
  if(debug && nx*ny<144)
  {
    cout << "Before transformation " << endl;
    print_array(*in[0]);
  }
  /// @typedef 1D FFT engine to be use
  typedef fft1d_engine<complex_type, complex_type, TEST_FFT_ENG> fft1d_engine_t;
  double dt_essl_f=0.0, dt_essl_b=0.0, dt_essl_t_f=0.0, dt_essl_t_b=0.0;
  //1D+zgetmo+1D
  {
    complex_type dummy;
    fft1d_engine_t myfft_xy;
    myfft_xy.set_defaults(ny,nx);
    myfft_xy.create(in[0]->data());
    fft1d_engine_t myfft_yx;
    myfft_yx.set_defaults(nx,ny);
    myfft_yx.create(in[0]->data());
    size_t nxy=nx*ny;
    Timer myclock;
    double my_dt_f=0.0, my_dt_b=0.0;
    for(int iter=0; iter<niters; ++iter)
    {
      for(int i=0; i<howmany; ++i)
      {
        complex_type* in_ptr=in[i]->data();
        complex_type* in_t_ptr=in_t[i]->data();
        myclock.restart();
        myfft_xy.fft_forward(in_ptr);
        Transpose2D<complex_type,TEST_TRANSPOSER>::apply(*in[i],*in_t[i]);
        myfft_yx.fft_forward(in_t_ptr);
        my_dt_f += myclock.elapsed();
        myclock.restart();
        myfft_yx.fft_backward(in_t_ptr);
        Transpose2D<complex_type,TEST_TRANSPOSER>::apply(*in_t[i],*in[i]);
        myfft_xy.fft_backward(in_ptr);
        my_dt_b += myclock.elapsed();
      }
    }
    if(debug)
      if(check_array(in[0]->data(),incopy.data(),incopy.size(),1.0/static_cast<double>(nx*ny)))
        cout << "We are good with 1D FFT+getmo+1D FFT" << endl;
    dt_essl_t_f=my_dt_f/static_cast<double>(niters);
    dt_essl_t_b=my_dt_b/static_cast<double>(niters);
  }
  for(int i=0; i<howmany; ++i)
    *(in[i])=incopy;
  //1D + transposed(1D), where the yx transpose is done on Tr(in) -> out
  {
    fft1d_engine_t myfft_xy, myfft_yx;
    myfft_xy.set_defaults(ny,nx);
    myfft_xy.create(in[0]->data());
    myfft_yx.set_defaults(nx,ny);
    myfft_yx.transpose(FFTW_FORWARD);
    myfft_yx.create(in[0]->data(),in_t[0]->data());
    size_t nxy=nx*ny;
    Timer myclock;
    double my_dt_f=0.0, my_dt_b=0.0;
    for(int iter=0; iter<niters; ++iter)
    {
      for(int i=0; i<howmany; ++i)
      {
        complex_type* in_ptr=in[i]->data();
        complex_type* in_t_ptr=in_t[i]->data();
        myclock.restart();
        myfft_xy.fft_forward(in_ptr);
        myfft_yx.fft_forward(in_ptr,in_t_ptr);
        my_dt_f += myclock.elapsed();
        myclock.restart();
        myfft_yx.fft_backward(in_t_ptr,in_ptr);
        myfft_xy.fft_backward(in_ptr);
        my_dt_b += myclock.elapsed();
        in_ptr += nxy;
        in_t_ptr += nxy;
      }
    }
    if(debug)
      if(check_array(in[0]->data(),incopy.data(),incopy.size(),1.0/static_cast<double>(nx*ny)))
        cout << "We are good with 1D FFT+t(1D FFT)" << endl;
    dt_essl_f=my_dt_f/static_cast<double>(niters);
    dt_essl_b=my_dt_b/static_cast<double>(niters);
  }
  for(int i=0; i<howmany; ++i)
    *(in[i])=incopy;
  //1D + transpose + 1D using transpose function
  double dt_f=0.0, dt_b=0.0;
  #pragma omp parallel
  {
    int np=omp_get_num_threads();
    int ip=omp_get_thread_num();
    complex_type dummy;
    fft1d_engine_t myfft_xy;
    myfft_xy.set_defaults(ny,nx/np);
    myfft_xy.create(&dummy);
    fft1d_engine_t myfft_yx;
    myfft_yx.set_defaults(nx,ny/np);
    myfft_yx.create(&dummy);
    int offset_x=ip*(nx/np)*ny;
    int offset_y=ip*(ny/np)*nx;
    int ny_b=ny/np;
    int nx_b=nx/np;
    int first_x=nx_b*ip;
    int last_x=nx_b*(ip+1);
    int first_y=ny_b*ip;
    int last_y=ny_b*(ip+1);
    size_t nxy=nx*ny;
    Timer myclock;
    double my_dt_f=0.0, my_dt_b=0.0;
    for(int iter=0; iter<niters; ++iter)
    {
      for(int i=0; i<howmany; ++i)
      {
        myfft_xy.fft_forward(in[i]->data()+offset_x);
        transpose_1(nx,ny,first_x,last_x,in[i]->data(), in_t[i]->data());
        //for(int k=0; k<ny; ++k)
        //  for(int j=first_x; j<last_x; ++j)
        //    in_x(i,k,j)=in(i,j,k);
      }
      #pragma omp barrier
      for(int i=0; i<howmany; ++i)
        myfft_yx.fft_forward(in_t[i]->data()+offset_y);
      my_dt_f+=myclock.elapsed();
      //////////////////////////////////////
      myclock.restart();
      for(int i=0; i<howmany; ++i)
      {
        myfft_yx.fft_backward(in_t[i]->data()+offset_y);
        transpose_1(ny,nx,first_y,last_y,in_t[i]->data(), in[i]->data());
        //for(int k=0; k<nx; ++k)
        //  for(int j=first_y; j<last_y; ++j)
        //    in(i,k,j)=in_x(i,j,k);
      }
      #pragma omp barrier
      for(int i=0; i<howmany; ++i)
        myfft_xy.fft_backward(in[i]->data()+offset_x);
      my_dt_b+=myclock.elapsed();
    }
    #pragma omp critical
    {
      dt_f += my_dt_f;
      dt_b += my_dt_b;
    }
  }
  dt_f /= static_cast<double>(niters*omp_get_max_threads());
  dt_b /= static_cast<double>(niters*omp_get_max_threads());
  if(debug)
    if(check_array(in[0]->data(),incopy.data(),incopy.size(),1.0/static_cast<double>(nx*ny)))
      cout << "We are good with 1D FFT+tfunc+ 1DFFT" << endl;
  printf("tag nX  nY  M  OMP 1D+t(1D) 1D+tfunc+1D 1D+getmo+1D  \n");
  printf("fft2d %d %d %d %d %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e\n"
         , nx, ny, howmany, omp_get_max_threads()
         , dt_essl_f, dt_essl_b, dt_f, dt_b, dt_essl_t_f, dt_essl_t_b);
  for(int i=0; i<howmany; ++i)
  {
    delete in[i];
    delete in_t[i];
  }
  OHMMS::Controller->finalize();
  return 0;
}

/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1770 $   $Date: 2007-02-17 17:45:38 -0600 (Sat, 17 Feb 2007) $
 * $Id: fft2d.cpp 1770 2007-02-17 23:45:38Z jnkim $
 ***************************************************************************/

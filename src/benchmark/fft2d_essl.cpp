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
#include <benchmark/fft_help.h>

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
  int howmany=1;
  int niters=10;
  int nx=6;
  int ny=6;
  //accepted: r2c  or c2c
  string fft_type("r2c");
  string fft_eng("fftw");
  bool inplace=true;
  bool single_precision=false;
  int opt;
  while((opt = getopt(argc, argv, "hp:x:y:m:i:e:t:")) != -1)
  {
    switch(opt)
    {
    case 'h':
      print_help("Help Message");
      return 1;
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
  int nx_b=nx, ny_b=ny;
  int nxy=nx*ny;
  int inc1x=1;
  int inc2x=ny;
  int inc1y=1;
  int inc2y=ny;
  //int inc1y=nx;
  //int inc2y=1;
  typedef complex<double> complex_type;
  Array<complex_type,3> in(howmany,nx,ny), out, incopy(howmany,nx,ny);
  if(inc1y==1)
    out.resize(howmany,nx,ny);
  else
    out.resize(howmany,ny,nx);
  init_array(in);
  cout << "Before transformation " << endl;
  if(nx*ny<64)
    print_array(in);
  //randomize(in.data(),in.size(),Random);
  incopy=in;
  int nbig=std::max(nx,ny);
  int naux1= (nbig<= 2048)?40000:static_cast<int>(40000+2.28*(nx+ny));
  int naux2= (nbig< 252)?20000:static_cast<int>(20000+(2*nbig+256)*(std::min(nx,ny)+2.28));
  Vector<double> aux1_f(naux1),aux1_b(naux1);
  Vector<double> aux2_f(naux2),aux2_b(naux2);
  double scale=1.0;
  int f_dir=ESSL_FFT_FORWARD;
  int b_dir=ESSL_FFT_BACKWARD;
  double iscale=1.0/static_cast<double>(nx*ny);
  dcft2 (1, in.data(), inc1x, inc2x, out.data(), inc1y, inc2y, ny, nx
         , f_dir, scale, aux1_f.data(), naux1, aux2_f.data(), naux2);
  dcft2 (1, out.data(), inc1y, inc2y, in.data(), inc1x, inc2x, ny, nx
         , b_dir, iscale, aux1_b.data(), naux1, aux2_b.data(), naux2);
  for(int i=0; i<howmany; ++i)
  {
    dcft2 (0, in.data()+i*nxy, inc1x, inc2x, out.data()+i*nxy, inc1y, inc2y, ny,nx
           , f_dir, scale, aux1_f.data(), naux1, aux2_f.data(), naux2);
  }
  cout << "After transformation " << endl;
  if(nx*ny<64)
    print_array(out);
  for(int i=0; i<howmany; ++i)
  {
    dcft2 (0, out.data()+i*nxy, inc1y, inc2y, in.data()+i*nxy, inc1x, inc2x, ny,nx
           , b_dir, iscale, aux1_b.data(), naux1, aux2_b.data(), naux2);
  }
  if(check_array(in.data(),incopy.data(),nx*ny,1.0))
    cout << "All is good " << endl;
  else
    return 0;
  Timer clock, clock_big;
  double dt_f=0.0, dt_b=0.0;
  for(int k=0; k<niters; ++k)
  {
    clock.restart();
    for(int i=0; i<howmany; ++i)
      dcft2 (0, in.data()+i*nxy, inc1x, inc2x, out.data()+i*nxy, inc1y, inc2y, ny,nx
             , f_dir, scale, aux1_f.data(), naux1, aux2_f.data(), naux2);
    dt_f += clock.elapsed();
    clock.restart();
    for(int i=0; i<howmany; ++i)
      dcft2 (0, out.data()+i*nxy, inc1y, inc2y, in.data()+i*nxy, inc1x, inc2x, ny,nx
             , b_dir, iscale, aux1_b.data(), naux1, aux2_b.data(), naux2);
    dt_b += clock.elapsed();
  }
  double t_norm=1.0/static_cast<double>(niters);
  cout << "Timer said = " << clock_big.elapsed()*t_norm << endl;
  printf("tag nX  nY  M  OMP 1D+t(1D) 1D+tfunc+1D 1D+getmo+1D  \n");
  printf("fft2d %d %d %d %d %12.4e %12.4e \n"
         , nx, ny, howmany, omp_get_max_threads()
         , dt_f*t_norm, dt_b*t_norm);
  OHMMS::Controller->finalize();
  return 0;
}


/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1770 $   $Date: 2007-02-17 17:45:38 -0600 (Sat, 17 Feb 2007) $
 * $Id: OrbitalBase.h 1770 2007-02-17 23:45:38Z jnkim $
 ***************************************************************************/

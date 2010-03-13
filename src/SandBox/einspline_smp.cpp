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
/**@file einspline_smp.cpp
 * @brief Benchmark einspline. Shared engine among threads.
 */
#include <SandBox/einspline_benchmark.h>
#include <Utilities/OhmmsInfo.h>
#include <Message/Communicate.h>
#include <getopt.h>
using namespace qmcplusplus;

int main(int argc, char** argv)
{

  OHMMS::Controller->initialize(argc,argv);
  Communicate* mycomm=OHMMS::Controller;
  OhmmsInfo Welcome(argc,argv,mycomm->rank());
  qmcplusplus::Random.init(0,1,11);

  int nx=48,ny=48,nz=48;
  int num_splines=128;
  int nsamples=512;
  int niters=10;

  int opt;
  while((opt = getopt(argc, argv, "hg:x:y:z:i:s:p:")) != -1) {
    switch(opt) {
      case 'h':
        printf("[-g grid| -x grid_x -y grid_y -z grid_z] -s states -p particles -i iterations\n");
        return 1;
      case 'g':
        nx=ny=nz=atoi(optarg);
        break;
      case 'x':
        nx=atoi(optarg);
        break;
      case 'y':
        ny=atoi(optarg);
        break;
      case 'z':
        nz=atoi(optarg);
        break;
      case 's':
        num_splines=atoi(optarg);
        break;
      case 'p':
        nsamples=atoi(optarg);
        break;
      case 'i':
        niters=atoi(optarg);
        break;
    }
  }

  if(mycomm->rank()==0)
  {
    cout << "#einspline benchmark grid = " << nx << " " << ny << " " << nz
      << " num_splines = " << num_splines << " num_samples = " << nsamples 
      << " iterations = " << niters<< endl;
    cout << "#MPI = " << mycomm->size() << "  OMP_NUM_THREADS = " << omp_get_max_threads() << endl;
    cout << "#datatype                       value             vgl              vgh " << endl;
  }

  typedef TinyVector<double,3> timer_type;
  timer_type d_timer_t,s_timer_t,z_timer_t,c_timer_t;

  einspline3d_benchmark<multi_UBspline_3d_z> z_bench;
  z_bench.set(nx,ny,nz,num_splines);

  einspline3d_benchmark<multi_UBspline_3d_d> d_bench;
  d_bench.set(nx,ny,nz,num_splines);

  einspline3d_benchmark<multi_UBspline_3d_s> s_bench;
  s_bench.set(nx,ny,nz,num_splines);

  einspline3d_benchmark<multi_UBspline_3d_c> c_bench;
  c_bench.set(nx,ny,nz,num_splines);

#pragma omp parallel 
  {
    random_position_generator<double> d_pos(nsamples,omp_get_thread_num());
    random_position_generator<float> s_pos(nsamples,omp_get_thread_num());
    timer_type d_timer,s_timer,z_timer,c_timer;
    for(int i=0; i<niters; ++i)
    {
      d_pos.randomize();
      d_timer+=d_bench.test_all(d_pos.Vpos, d_pos.VGLpos, d_pos.VGHpos);
      s_pos.randomize();
      s_timer+=s_bench.test_all(s_pos.Vpos, s_pos.VGLpos, s_pos.VGHpos);
      d_pos.randomize();
      z_timer+=z_bench.test_all(d_pos.Vpos, d_pos.VGLpos, d_pos.VGHpos);
      s_pos.randomize();
      c_timer+=c_bench.test_all(s_pos.Vpos, s_pos.VGLpos, s_pos.VGHpos);
    }
#pragma omp critical
    {
      d_timer_t += d_timer;
      s_timer_t += s_timer;
      z_timer_t += z_timer;
      c_timer_t += c_timer;
    }
  }

  double nops=num_splines*nsamples*niters*omp_get_max_threads();
  cout.setf(std::ios::scientific, std::ios::floatfield);
  cout.precision(6);

  cout << "evals "<< nops << setw(4) << omp_get_max_threads() <<  " double    "<< setw(15) << nops/d_timer_t[0] << setw(15) << nops/d_timer_t[1]<< setw(15) << nops/d_timer_t[2] << endl;
  cout << "evals "<< nops << setw(4) << omp_get_max_threads() <<  " single    "<< setw(15) << nops/s_timer_t[0] << setw(15) << nops/s_timer_t[1]<< setw(15) << nops/s_timer_t[2] << endl;
  cout << "evals "<< nops << setw(4) << omp_get_max_threads() <<  " d-complex "<< setw(15) << nops/z_timer_t[0] << setw(15) << nops/z_timer_t[1]<< setw(15) << nops/z_timer_t[2] << endl;
  cout << "evals "<< nops << setw(4) << omp_get_max_threads() <<  " s-complex "<< setw(15) << nops/c_timer_t[0] << setw(15) << nops/c_timer_t[1]<< setw(15) << nops/c_timer_t[2] << endl;

  double tfac=1.0/static_cast<double>(omp_get_max_threads());
  d_timer_t*=tfac;
  s_timer_t*=tfac;
  z_timer_t*=tfac;
  c_timer_t*=tfac;
  cout << "wtime " << nops << setw(4) << omp_get_max_threads() << " double    " << setw(15) << d_timer_t[0]  << setw(15) << d_timer_t[1] << setw(15) << d_timer_t[2] << endl;
  cout << "wtime " << nops << setw(4) << omp_get_max_threads() << " single    " << setw(15) << s_timer_t[0]  << setw(15) << s_timer_t[1] << setw(15) << s_timer_t[2] << endl;
  cout << "wtime " << nops << setw(4) << omp_get_max_threads() << " d-complex " << setw(15) << z_timer_t[0]  << setw(15) << z_timer_t[1] << setw(15) << z_timer_t[2] << endl;
  cout << "wtime " << nops << setw(4) << omp_get_max_threads() << " s-complex " << setw(15) << c_timer_t[0]  << setw(15) << c_timer_t[1] << setw(15) << c_timer_t[2] << endl;

  OHMMS::Controller->finalize();
  return 0;
}

/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1770 $   $Date: 2007-02-17 17:45:38 -0600 (Sat, 17 Feb 2007) $
 * $Id: OrbitalBase.h 1770 2007-02-17 23:45:38Z jnkim $ 
 ***************************************************************************/

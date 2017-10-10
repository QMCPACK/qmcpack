//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



/**@file multidet.cpp
 * @brief Test codes for multidets
 */
#include <Configuration.h>
#include <SandBox/einspline_benchmark.h>
#include <Message/Communicate.h>
#include <mpi/collectives.h>
#include <getopt.h>
using namespace qmcplusplus;

int main(int argc, char** argv)
{
  OHMMS::Controller->initialize(argc,argv);
  Communicate* mycomm=OHMMS::Controller;
  OhmmsInfo Welcome("einspline",mycomm->rank());
  Random.init(0,1,11);
  int nx=48,ny=48,nz=48;
  int num_splines=128;
  int nsamples=512;
  int niters=10;
  int opt;
  while((opt = getopt(argc, argv, "hg:x:y:z:i:s:p:")) != -1)
  {
    switch(opt)
    {
    case 'h':
      printf("[-g grid| -x grid_x -y grid_y -z grid_z] -s states -p particles -i iterations -t [d|s|z|c] \n");
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
  app_log() << "#einspline benchmark grid = " << nx << " " << ny << " " << nz
            << " num_splines = " << num_splines << " num_samples = " << nsamples
            << " iterations = " << niters
            << " number of operations in millions " << std::endl;
  app_log() << "#MPI = " << mycomm->size() << "  OMP_NUM_THREADS = " << omp_get_max_threads() << std::endl;
  app_log() << "#   mpi   openmp    datatype     "
            << "  value_op         vgl_op              vgh_op         value_time       vgl_time         vgh_time" << std::endl;
  app_log().flush();
  typedef TinyVector<double,3> timer_type;
  timer_type d_timer_t(0.0),s_timer_t(0.0),z_timer_t(0.0),c_timer_t(0.0);
  #pragma omp parallel
  {
    einspline3d_benchmark<multi_UBspline_3d_d> d_bench;
    d_bench.set(nx,ny,nz,num_splines);
    random_position_generator<double> d_pos(nsamples,omp_get_thread_num());
    timer_type d_timer(0.0);
    for(int i=0; i<niters; ++i)
    {
      d_pos.randomize();
      d_timer+=d_bench.test_all(d_pos.Vpos, d_pos.VGLpos, d_pos.VGHpos);
    }
    #pragma omp critical
    {
      d_timer_t += d_timer;
    }
  }
//#pragma omp parallel
//  {
//    einspline3d_benchmark<multi_UBspline_3d_z> d_bench;
//    d_bench.set(nx,ny,nz,num_splines);
//    random_position_generator<double> d_pos(nsamples,omp_get_thread_num());
//    timer_type d_timer(0.0);
//    for(int i=0; i<niters; ++i)
//    {
//      d_pos.randomize();
//      d_timer+=d_bench.test_all(d_pos.Vpos, d_pos.VGLpos, d_pos.VGHpos);
//    }
//#pragma omp critical
//    {
//      z_timer_t += d_timer;
//    }
//  }
  #pragma omp parallel
  {
    einspline3d_benchmark<multi_UBspline_3d_s> d_bench;
    d_bench.set(nx,ny,nz,num_splines);
    random_position_generator<float> d_pos(nsamples,omp_get_thread_num());
    timer_type d_timer(0.0);
    for(int i=0; i<niters; ++i)
    {
      d_pos.randomize();
      d_timer+=d_bench.test_all(d_pos.Vpos, d_pos.VGLpos, d_pos.VGHpos);
    }
    #pragma omp critical
    {
      s_timer_t += d_timer;
    }
  }
//#pragma omp parallel
//  {
//    einspline3d_benchmark<multi_UBspline_3d_c> d_bench;
//    d_bench.set(nx,ny,nz,num_splines);
//    random_position_generator<float> d_pos(nsamples,omp_get_thread_num());
//    timer_type d_timer(0.0);
//    for(int i=0; i<niters; ++i)
//    {
//      d_pos.randomize();
//      d_timer+=d_bench.test_all(d_pos.Vpos, d_pos.VGLpos, d_pos.VGHpos);
//    }
//#pragma omp critical
//    {
//      c_timer_t += d_timer;
//    }
//  }
  mpi::reduce(*mycomm,d_timer_t);
  mpi::reduce(*mycomm,z_timer_t);
  mpi::reduce(*mycomm,s_timer_t);
  mpi::reduce(*mycomm,c_timer_t);
  double nops=num_splines*nsamples*1.e-6;
  double tfac=1.0/static_cast<double>(mycomm->size()*omp_get_max_threads()*niters);
  d_timer_t*=tfac;
  s_timer_t*=tfac;
  z_timer_t*=tfac;
  c_timer_t*=tfac;
  app_log().setf(std::ios::scientific, std::ios::floatfield);
  app_log().precision(6);
  app_log() << "einspline "<< std::setw(4) << mycomm->size()<< std::setw(4) << omp_get_max_threads() <<  " double    "<< nops/d_timer_t << d_timer_t << std::endl;
  app_log() << "einspline "<< std::setw(4) << mycomm->size()<< std::setw(4) << omp_get_max_threads() <<  " d-complex "<< nops/z_timer_t << z_timer_t << std::endl;
  app_log() << "einspline "<< std::setw(4) << mycomm->size()<< std::setw(4) << omp_get_max_threads() <<  " single    "<< nops/s_timer_t << s_timer_t << std::endl;
  app_log() << "einspline "<< std::setw(4) << mycomm->size()<< std::setw(4) << omp_get_max_threads() <<  " s-complex "<< nops/c_timer_t << c_timer_t << std::endl;
  OHMMS::Controller->finalize();
  return 0;
}


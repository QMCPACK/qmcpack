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
    
    



#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include "Numerics/MatrixOperators.h"
#include "Utilities/RandomGenerator.h"
#include "Utilities/Timer.h"
using namespace qmcplusplus;

int main(int argc, char** argv)
{
  Random.init(0,1,-1);
  typedef double scalar_t;
  int N=100;
  int M=2000;
  int L=5;
  Matrix<scalar_t> A(N,M), B(M,L), C(N,L);
  Vector<scalar_t> y(M),z(N);
  double gemm_t=0.0;
  double gemv_t=0.0;
  Timer myTimer;
  int niter=100;
  for(int iter=0; iter<niter; iter++)
  {
    for(int i=0; i<N; i++)
      for(int j=0; j<M; j++)
        A(i,j)=Random();
    //A(i,j)=scalar_t(Random(),Random());
    for(int i=0; i<M; i++)
      for(int j=0; j<L; j++)
        B(i,j)=Random();
    //B(i,j)=scalar_t(Random(),Random());
    for(int i=0; i<M; i++)
      y(i)=Random();
    //y(i)=scalar_t(Random(),Random());
    myTimer.restart();
    MatrixOperators::product(A,B,C);
    gemm_t+=myTimer.elapsed();
    myTimer.restart();
    MatrixOperators::product(A,y,z.data());
    gemv_t+=myTimer.elapsed();
    // std::cout << "<<<<<< GEMM TEST " << std::endl;
    // for(int i=0; i<N; i++) {
    //   for(int j=0; j<L; j++) {
    //     scalar_t v=0.0;
    //     for(int k=0; k<M; k++) v+=A(i,k)*B(k,j);
    //     std::cout << i << "," << j << " " << v-C(i,j) << " " << v << std::endl;
    //   }
    // }
    // std::cout << "<<<<<< GEMV TEST " << std::endl;
    // for(int i=0; i<N; i++) {
    //   scalar_t v=0.0;
    //   for(int j=0; j<M; j++) v+=A(i,j)*y(j);
    //   std::cout << i << " " << v-z[i] << " " << v << std::endl;
    // }
  }
  std::cout << "|"<<N << "|"<<M << "|" << L << "|" << gemm_t/niter << "|" << gemv_t/niter<< "|" << std::endl;
}

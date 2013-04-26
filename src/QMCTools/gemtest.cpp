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
    // cout << "<<<<<< GEMM TEST " << endl;
    // for(int i=0; i<N; i++) {
    //   for(int j=0; j<L; j++) {
    //     scalar_t v=0.0;
    //     for(int k=0; k<M; k++) v+=A(i,k)*B(k,j);
    //     cout << i << "," << j << " " << v-C(i,j) << " " << v << endl;
    //   }
    // }
    // cout << "<<<<<< GEMV TEST " << endl;
    // for(int i=0; i<N; i++) {
    //   scalar_t v=0.0;
    //   for(int j=0; j<M; j++) v+=A(i,j)*y(j);
    //   cout << i << " " << v-z[i] << " " << v << endl;
    // }
  }
  cout << "|"<<N << "|"<<M << "|" << L << "|" << gemm_t/niter << "|" << gemv_t/niter<< "|" << endl;
}

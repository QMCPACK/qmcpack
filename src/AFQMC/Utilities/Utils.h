#ifndef QMCPLUSPLUS_AFQMC_UTILS_H
#define QMCPLUSPLUS_AFQMC_UTILS_H

#include<iostream>
#include<cstdlib>
#include<vector>
#include<string>
#include<algorithm>

#include<AFQMC/config.h>

namespace qmcplusplus { 

namespace afqmc {


// careful 
// FIX FIX FIX
// this routine returns interchanged (i,j)/(k,l), so it is wrong due to a std::complex conjugation for std::complex matrix elements 
int cntExcitations(int NAEA, int NAEB, std::vector<IndexType>& DL, std::vector<IndexType>& DR, IndexType& n0, IndexType& n1, IndexType& n2, IndexType& n3, std::vector<IndexType>& occ, RealType& sg);

template<class Vec,
        class RandomNumberGenerator_,
        typename = typename std::enable_if_t<std::decay<Vec>::type::dimensionality==1>
        > 
void sampleGaussianFields(Vec&& V, RandomNumberGenerator_& rng)
{
  size_t n = V.size();
  for (int i=0; i+1<n; i+=2)
  {
    RealType temp1=1-0.9999999999*rng(), temp2=rng();
    RealType mag = std::sqrt(-2.0*std::log(temp1));
    V[i]  =mag*std::cos(6.283185306*temp2);
    V[i+1]=mag*std::sin(6.283185306*temp2);
  }
  if (n%2==1)
  {
    RealType temp1=1-0.9999999999*rng(), temp2=rng();
    V[n-1]=std::sqrt(-2.0*std::log(temp1))*std::cos(6.283185306*temp2);
  }
}

template<class Mat,
        class RandomNumberGenerator_,
        typename = typename std::enable_if_t<(std::decay<Mat>::type::dimensionality>1)>,
        typename = void
        >
void sampleGaussianFields(Mat&& M, RandomNumberGenerator_& rng)
{
  for(int i=0, iend=M.shape()[0]; i<iend; ++i)
    sampleGaussianFields(M[i],rng);
}

template<class T,
        class RandomNumberGenerator_>
void sampleGaussianFields_n(T* V, int n, RandomNumberGenerator_& rng)
{
  for (int i=0; i+1<n; i+=2)
  {
    RealType temp1=1-0.9999999999*rng(), temp2=rng();
    RealType mag = std::sqrt(-2.0*std::log(temp1));
    V[i]  =mag*std::cos(6.283185306*temp2);
    V[i+1]=mag*std::sin(6.283185306*temp2);
  }
  if (n%2==1)
  {
    RealType temp1=1-0.9999999999*rng(), temp2=rng();
    V[n-1]=std::sqrt(-2.0*std::log(temp1))*std::cos(6.283185306*temp2);
  }
}



}

}

#endif

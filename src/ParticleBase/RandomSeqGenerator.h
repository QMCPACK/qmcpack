//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#ifndef QMCPLUSPLUS_RANDOMSEQUENCEGENERATOR_H
#define QMCPLUSPLUS_RANDOMSEQUENCEGENERATOR_H
#include <algorithm>
#if defined(HAVE_LIBBLITZ)
#include <blitz/array.h>
#endif
#include "OhmmsPETE/OhmmsMatrix.h"
#include "ParticleBase/ParticleAttrib.h"
#include "Utilities/RandomGenerator.h"

/*!\fn template<class T> void assignGaussRand(T* restrict a, unsigned n)
  *\param a the starting pointer
  *\param n the number of type T to be assigned
  *\brief Assign Gaussian distributed random numbers using Box-Mueller algorithm. Called by overloaded funtions makeGaussRandom
  */
//template<class T>
//inline void assignGaussRand(T* restrict a, unsigned n) {
//  for (int i=0; i+1<n; i+=2) {
//    T temp1=1-0.9999999999*Random(), temp2=Random();
//    a[i]  =sqrt(-2.0*log(temp1))*cos(6.283185306*temp2);
//    a[i+1]=sqrt(-2.0*log(temp1))*sin(6.283185306*temp2);
//  }
//  if (n%2==1) {
//    T temp1=1-0.9999999999*Random(), temp2=Random();
//    a[n-1]=sqrt(-2.0*log(temp1))*cos(6.283185306*temp2);
//  }
//}

namespace qmcplusplus
{

template<class T, class RG>
inline void assignGaussRand(T* restrict a, unsigned n, RG& rng)
{
  OHMMS_PRECISION_FULL slightly_less_than_one = 1.0 -
    std::numeric_limits<OHMMS_PRECISION_FULL>::epsilon();
  int nm1=n-1;
  OHMMS_PRECISION_FULL temp1, temp2;
  for (int i=0; i<nm1; i+=2)
  {
    temp1=std::sqrt(-2.0*std::log(1.0-slightly_less_than_one*rng()));
    temp2=2.0*M_PI*rng();
    a[i]  =temp1*std::cos(temp2);
    a[i+1]=temp1*std::sin(temp2);
  }
  if (n%2==1)
  {
    temp1=std::sqrt(-2.0*std::log(1.0-slightly_less_than_one*rng()));
    temp2=2.0*M_PI*rng();
    a[nm1]=temp1*std::cos(temp2);
  }
}

/*!\fn template<class T> void assignUniformRand(T* restrict a, unsigned n)
  *\param a the starting pointer
  *\param n the number of type T to be assigned
  *\brief Assign unifor distributed random numbers [0,1)
  */
//template<class T>
//inline void assignUniformRand(T* restrict a, unsigned n) {
//  for(int i=0; i<n;i++) a[i] = Random();
//  //generate(a,a+n,Random);
//}
template<class T, class RG>
inline void assignUniformRand(T* restrict a, unsigned n, RG& rng)
{
  for(int i=0; i<n; i++)
    a[i] = rng();
}

#if defined(HAVE_LIBBLITZ)
///specialized functions: stick to overloading
template<typename T, unsigned D>
inline void makeGaussRandom(blitz::Array<TinyVector<T,D>, 2>& a)
{
  assignGaussRand(&(a(0,0)[0]), a.size()*D, Random);
}
#endif

template<typename T, unsigned D>
inline void makeGaussRandom(std::vector<TinyVector<T,D> >& a)
{
  assignGaussRand(&(a[0][0]), a.size()*D, Random);
}

///specialized functions: stick to overloading
template<typename T, unsigned D>
inline void makeGaussRandom(Matrix<TinyVector<T,D> >& a)
{
  assignGaussRand(&(a(0,0)[0]), a.size()*D, Random);
}

template<typename T, unsigned D>
inline void makeGaussRandom(ParticleAttrib<TinyVector<T,D> >& a)
{
  assignGaussRand(&(a[0][0]), a.size()*D, Random);
}

template<typename T>
inline void makeGaussRandom(ParticleAttrib<T>& a)
{
  assignGaussRand(&(a[0]), a.size(), Random);
}

template<typename T, unsigned D>
inline void makeUniformRandom(ParticleAttrib<TinyVector<T,D> >& a)
{
  assignUniformRand(&(a[0][0]), a.size()*D, Random);
}

template<typename T>
inline void makeUniformRandom(ParticleAttrib<T>& a)
{
  assignUniformRand(&(a[0]), a.size(), Random);
}

template<typename T>
inline void makeSphereRandom(ParticleAttrib<TinyVector<T,3> >& a)
{
  for(int i=0; i<a.size(); i++)
  {
    bool failed=true;
    while(failed)
    {
      T x=1.0-2.0*Random();
      T y=1.0-2.0*Random();
      T z=1.0-2.0*Random();
      T sep=std::sqrt(x*x+y*y+z*z);
      if(sep<1)
      {
        T rinv=1.0/sep;
        a[i][0]=x*rinv;
        a[i][1]=y*rinv;
        a[i][2]=z*rinv;
        failed=false;
      }
    }
  }
}

template<typename T>
inline void makeSphereRandom(ParticleAttrib<TinyVector<T,2> >& a)
{
  for(int i=0; i<a.size(); i++)
  {
    bool failed=true;
    while(failed)
    {
      T x=1.0-2.0*Random();
      T y=1.0-2.0*Random();
      T sep=std::sqrt(x*x+y*y);
      if(sep<1)
      {
        T rinv=1.0/sep;
        a[i][0]=x*rinv;
        a[i][1]=y*rinv;
        failed=false;
      }
    }
  }
}

template<typename T, unsigned D, class RG>
inline void makeGaussRandomWithEngine(ParticleAttrib<TinyVector<T,D> >& a, RG& rng)
{
  assignGaussRand(&(a[0][0]), a.size()*D, rng);
}

template<typename T, unsigned D, class RG>
inline void makeGaussRandomWithEngine(std::vector<TinyVector<T,D> >& a, RG& rng)
{
  assignGaussRand(&(a[0][0]), a.size()*D, rng);
}


}

#endif

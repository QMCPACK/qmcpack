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
template<class T> 
inline void assignGaussRand(T* restrict a, unsigned n) {
  for (int i=0; i+1<n; i+=2) {
    T temp1=1-0.9999999999*Random(), temp2=Random();
    a[i]  =sqrt(-2.0*log(temp1))*cos(6.283185306*temp2);
    a[i+1]=sqrt(-2.0*log(temp1))*sin(6.283185306*temp2);
  }
  if (n%2==1) {
    T temp1=1-0.9999999999*Random(), temp2=Random();
    a[n-1]=sqrt(-2.0*log(temp1))*cos(6.283185306*temp2);
  }
}

/*!\fn template<class T> void assignUniformRand(T* restrict a, unsigned n)
  *\param a the starting pointer
  *\param n the number of type T to be assigned
  *\brief Assign unifor distributed random numbers [0,1)
  */
template<class T>
inline void assignUniformRand(T* restrict a, unsigned n) {
  for(int i=0; i<n;i++) a[i] = Random();
  //generate(a,a+n,Random);
}

#if defined(HAVE_LIBBLITZ)
///specialized functions: stick to overloading
inline void makeGaussRandom(blitz::Array<TinyVector<double,3>, 2>& a) {
  assignGaussRand(&(a(0,0)[0]), a.size()*3);
}
#endif

///specialized functions: stick to overloading
inline void makeGaussRandom(Matrix<TinyVector<double,3> >& a) {
  assignGaussRand(&(a(0,0)[0]), a.size()*3);
}

inline void makeGaussRandom(ParticleAttrib<TinyVector<double,3> >& a) {
  assignGaussRand(&(a[0][0]), a.size()*3);
}

inline void makeGaussRandom(ParticleAttrib<double>& a) {
  assignGaussRand(&(a[0]), a.size());
}

inline void makeUniformRandom(ParticleAttrib<TinyVector<double,3> >& a) {
  assignUniformRand(&(a[0][0]), a.size()*3);
}

inline void makeUniformRandom(ParticleAttrib<double>& a) {
  assignUniformRand(&(a[0]), a.size());
}
#endif

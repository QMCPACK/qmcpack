#ifndef QMCPLUSPLUS_COULOMB_TYPE_HPP
#define QMCPLUSPLUS_COULOMB_TYPE_HPP
#include <config.h>
//model to implement the mixed precision
#if (__cplusplus >= 201103L)

#define DECLARE_COULOMB_TYPES \
  using pRealType=OHMMS_PRECISION;              \
  using mRealType=double;                       \
  using pComplexType=std::complex<pRealType>;   \
  using mComplexType=std::complex<double>;      \
  using PosType=TinyVector<pRealType,OHMMS_DIM>; 

#else

#define DECLARE_COULOMB_TYPES \
  typedef OHMMS_PRECISION pRealType;               \
  typedef double mRealType;                        \
  typedef std::complex<pRealType> pComplexType;    \
  typedef std::complex<double> mComplexType;       \
  typedef TinyVector<pRealType,OHMMS_DIM> PosType;

#endif

#endif

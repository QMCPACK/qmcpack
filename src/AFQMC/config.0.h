#ifndef AFQMC_CONFIG_0_H 
#define AFQMC_CONFIG_0_H 

#define ADD_TESTS_TIMERS

#define AFQMC_DEBUG 3 
#define AFQMC_TIMER 

#define MAXIMUM_EMPLACE_BUFFER_SIZE 102400 

// maximum size in Bytes for a dataset with walker data on WalkerIO
#define WALKER_HDF_BLOCK_SIZE 100000000 

// maximum size in Bytes for a block of data in CSR matrix HDF IO 
#define CSR_HDF_BLOCK_SIZE 2000000 

// careful here that RealType is consistent with this!!!
#define MKL_INT         int
#define MKL_Complex8    std::complex<float> 
#define MKL_Complex16   std::complex<double> 

#define byRows   999
#define byCols   111

#define PsiT_IN_SHM

namespace qmcplusplus
{
namespace afqmc
{

  typedef OHMMS_INDEXTYPE                 IndexType;
  typedef OHMMS_INDEXTYPE                 OrbitalType;
  typedef OHMMS_PRECISION_FULL            RealType;
  typedef OHMMS_PRECISION                 SPRealType;

#if defined(QMC_COMPLEX)
  typedef std::complex<RealType>         ValueType;
  typedef std::complex<SPRealType>       SPValueType;
#else
  typedef RealType                       ValueType;
  typedef SPRealType                     SPValueType;
#endif
  typedef std::complex<RealType>         ComplexType;
  typedef std::complex<SPRealType>       SPComplexType;

}
}


#endif

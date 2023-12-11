#ifndef AFQMC_CONFIG_0_H
#define AFQMC_CONFIG_0_H

#if defined __INTEL_COMPILER
#pragma warning disable 2196
#endif

#define BOOST_NO_AUTO_PTR

#define ADD_TESTS_TIMERS

#define AFQMC_DEBUG 3
#define AFQMC_TIMER

#define MAXIMUM_EMPLACE_BUFFER_SIZE 102400

// maximum size in Bytes for a dataset with walker data on WalkerIO
#define WALKER_HDF_BLOCK_SIZE 100000000

// maximum size in Bytes for a block of data in CSR matrix HDF IO
#define CSR_HDF_BLOCK_SIZE 2000000

// careful here that RealType is consistent with this!!!
#define MKL_INT int
#define MKL_Complex8 std::complex<float>
#define MKL_Complex16 std::complex<double>

#define byRows 999
#define byCols 111

// guard with directive that checks if boost version is >=1.65
// uncomment to enable stacktrace
//#define ENABLE_STACKTRACE
#include <boost/version.hpp>
#if (BOOST_VERSION >= 106500) && defined(ENABLE_STACKTRACE)
#include <boost/stacktrace.hpp>
#define print_stacktrace std::cout << boost::stacktrace::stacktrace();
#else
#define print_stacktrace std::cout << "stacktrace not enabled.\n";
#endif

namespace qmcplusplus
{
namespace afqmc
{
using IndexType   = OHMMS_INDEXTYPE;
using OrbitalType = OHMMS_INDEXTYPE;
using RealType    = OHMMS_PRECISION_FULL;
using SPRealType  = OHMMS_PRECISION;

#if defined(QMC_COMPLEX)
using ValueType   = std::complex<RealType>;
using SPValueType = std::complex<SPRealType>;
#else
using ValueType   = RealType;
using SPValueType = SPRealType;
#endif
using ComplexType   = std::complex<RealType>;
using SPComplexType = std::complex<SPRealType>;

} // namespace afqmc
} // namespace qmcplusplus


#endif

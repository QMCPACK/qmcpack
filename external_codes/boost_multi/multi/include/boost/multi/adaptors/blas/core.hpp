// Copyright 2019-2024 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_MULTI_ADAPTORS_BLAS_CORE_HPP
#define BOOST_MULTI_ADAPTORS_BLAS_CORE_HPP
#pragma once

// https://software.intel.com/en-us/articles/intel-mkl-link-line-advisor

#include <array>
#include<cassert>
#include<complex>
#include<cstdint>      // int64_t
#include<cstring>      // std::memcpy
#include<iostream>     // for debug
#include<limits>       // numeric_limits
#include<type_traits>  // is_convertible

// #include "../../config/MARK.hpp"

#include <boost/multi/adaptors/blas/traits.hpp>  // IWYU pragma: export

#if ! defined(NDEBUG)
	#include<stdexcept>
	#include<string>
	#define BOOST_MULTI_ASSERT1(ExpR)              (void)((ExpR)?0:throw std::logic_error("\n" __FILE__ ":"+std::to_string(__LINE__)+"::\n"+std::string(__PRETTY_FUNCTION__)+"\nLogic assertion `" #ExpR "' failed.")) /*NOLINT(fuchsia-default-arguments-calls,cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay)*/
	#define BOOST_MULTI_ASSERT2(ExpR, DescriptioN) (void)((ExpR)?0:throw std::DescriptioN("\n" __FILE__ ":"+std::to_string(__LINE__)+"::\n"+std::string(__PRETTY_FUNCTION__)+"\nLogic assertion `" #ExpR "' failed."))
#else
	#define BOOST_MULTI_ASSERT1(ExpR)              assert(ExpR)
	#define BOOST_MULTI_ASSERT2(ExpR, DescriptioN) assert(EXpR)
#endif

#ifdef CBLAS_H
#define BLAS(NamE) cblas_##NamE
#else
#define BLAS(NamE) NamE##_

#ifndef MULTI_BLAS_INT
#if defined(__INTPTR_WIDTH__)
	#define MULTI_BLAS_INT __INTPTR_WIDTH__
#endif
#endif

#define s float
#define d double
#define c std::complex<s>
#define z std::complex<d>
#define v void

// cppcheck-suppress unusedStructMember
using Complex_float  = struct { float  real; float imag; };
// cppcheck-suppress unusedStructMember
using Complex_double = struct { double real; double imag; };

#define C Complex_float   // _Complex s
#define Z Complex_double  // _Complex d

#if defined(MULTI_BLAS_INT)
	#if   MULTI_BLAS_INT==32
		using INT = std::int32_t;  // #define INT int32_t
	#elif MULTI_BLAS_INT==64
		using INT = std::int64_t;  // #define INT int64_t
	#else
		using INT = std::int32_t;  // #define INT int32_t  // 32bit safe? pesimistic?
	#endif
#else
	using INT = std::int32_t;  // #define INT int32_t  // 32bit safe? pesimistic?
#endif

namespace core {
	using size_t = INT;
	using ssize_t = std::make_signed_t<size_t>;
}  // end namespace core

extern "C" {

#define INTEGER INT const&
#define N INTEGER n
#define INCX INTEGER incx
#define INCY INTEGER incy

static_assert(sizeof(INT)==32/8 || sizeof(INT)==64/8, "please set MULTI_BLAS_INT to int32_t or int64_t");

// indented declarations like in https://www.netlib.org/lapack/lug/node145.html

#define xROTG(T1, T2)     v       T1##rotg ##_ (                                                                T1 const*, T1 const*, T2*, T1*                          )  // NOLINT(bugprone-macro-parentheses,readability-identifier-length) macro arg expands to type
#define xROTMG(T)         v       T ##rotmg##_ (                                                        T*, T*, T*       , T  const&,                     T(&param)[5]  )  // NOLINT(bugprone-macro-parentheses,readability-identifier-length) macro arg expands to type
#define xROT(TT, T, S)    v       TT##rot  ##_ (    N,              T       *x, INCX, T       *y, INCY,                               S const&, S const&                )  // NOLINT(bugprone-macro-parentheses,readability-identifier-length) macro arg expands to type
#define xROTM(T)          v       T ##rotm ##_ (    N,              T       *x, INCX, T       *y, INCY,                                                   T const(&p)[5])  // NOLINT(bugprone-macro-parentheses,readability-identifier-length) macro arg expands to type
#define xSWAP(T)          v       T ##swap ##_ (    N,              T       *x, INCX, T       *y, INCY                                                                  )  // NOLINT(bugprone-macro-parentheses,readability-identifier-length) macro arg expands to type
#define xSCAL(TT, TA, TX) v       TT##scal ##_ (    N, TA const& a, TX      *x, INCX                                                                                    )  // NOLINT(bugprone-macro-parentheses,readability-identifier-length) macro arg expands to type
#define xCOPY(T)          v       T ##copy ##_ (    N,              T const *x, INCX, T       *y, INCY                                                                  )  // NOLINT(bugprone-macro-parentheses,readability-identifier-length) macro arg expands to type
#define xAXPY(T)          v       T ##axpy ##_ (    N,  T const* a, T const *x, INCX, T       *y, INCY                                                                  )  // NOLINT(bugprone-macro-parentheses,readability-identifier-length) macro arg expands to type

// in MKL and vanilla BLAS, OpenBLAS, real dot always return by stack
#define xDOT(R, TT, T)    auto    TT##dot  ##_ (    N,              T const *x, INCX, T const *y, INCY) -> R  // NOLINT(readability-identifier-length) conventional BLAS naming

// PGI/NVC++ compiler uses a blas version that needs -DRETURN_BY_STACK
#if defined(BLAS_DOT_RETURNS_VOID)
//#if defined(RETURN_BY_STACK) || (defined(FORTRAN_COMPLEX_FUNCTIONS_RETURN_VOID) && FORTRAN_COMPLEX_FUNCTIONS_RETURN_VOID)
//#define xDOT(R, TT, T)    v       TT##dot  ##_ (R*, N,              T const *x, INCX, T const *y, INCY)
#define xDOTU(R, T)       v       T ##dotu ##_ (R*, N,              T const * /*x*/, INCX, T const * /*y*/, INCY)  // NOLINT(bugprone-macro-parentheses) : macro arg expands to type
#define xDOTC(R, T)       v       T ##dotc ##_ (R*, N,              T const * /*x*/, INCX, T const * /*y*/, INCY)  // NOLINT(bugprone-macro-parentheses) : macro arg expands to type
#else
#define xDOTU(R, T)       auto    T ##dotu ##_ (    N,              T const *x, INCX, T const *y, INCY) -> R  // NOLINT(readability-identifier-length) conventional BLAS naming
#define xDOTC(R, T)       auto    T ##dotc ##_ (    N,              T const *x, INCX, T const *y, INCY) -> R  // NOLINT(readability-identifier-length) conventional BLAS naming
//#define xxDOT(TT, T)      auto    TT##dot  ##_ (    N,  T const& a, T const *x, INCX, T const *y, INCY) -> T
#endif

#define xNRM2(R, TT, T)   auto    TT##nrm2##_ (    N,               T const *x, INCX) -> R    // NOLINT(readability-identifier-length) conventional BLAS naming
#define xASUM(R, TT, T)   auto    TT##asum##_ (    N,               T const *x, INCX) -> R    // NOLINT(readability-identifier-length) conventional BLAS naming
#define IxAMAX(T)         auto i##T ##amax##_ (    N,               T const* x, INCX) -> INT  // NOLINT(readability-identifier-length) conventional BLAS naming

xROTG(s, s)   ; xROTG(d,d)    ;  // TODO(correaa) MKL extension for "(c, s)" and "(z, d)"?
xROTMG(s)     ; xROTMG(d)     ;
xROT(s, s, s) ; xROT(d, d, d) ;                 xROT(cs, c, s); xROT(zd, z, d);
xROTM(s)      ; xROTM(d)      ;
xSWAP(s)      ; xSWAP(d)      ; xSWAP(c)      ; xSWAP(z);
xSCAL(s, s, s); xSCAL(d, d, d); xSCAL(c, c, c); xSCAL(z, z, z); xSCAL(zd, d, z); xSCAL(cs, s, c);
xCOPY(s)      ; xCOPY(d)      ; xCOPY(c)      ; xCOPY(z)      ;
xAXPY(s)      ; xAXPY(d)      ; xAXPY(c)      ; xAXPY(z)      ;

xDOT(s, s, s) ; xDOT(d, d, d) ;                                   xDOT(d, ds, s);
xDOTU(C, c); xDOTU(Z, z);  // TODO(correaa) MKL extension for "(c, c)" and "(z, z)"?
xDOTC(C, c); xDOTC(Z, z);  // TODO(correaa) MKL extension for "(sds, s)"

xNRM2(s, s, s); xNRM2(d, d, d); xNRM2(s, sc, c); xNRM2(d, dz, z);
xASUM(s, s, s); xASUM(d, d, d); xASUM(s, sc, c); xASUM(d, dz, z);
IxAMAX(s); IxAMAX(d); IxAMAX(c); IxAMAX(z);

#define TRANS const char& trans
#define NR INTEGER nr
#define NC INTEGER nc
#define LDA INTEGER lda
#define UPLO const char& uplo
#define DIAG const char& diag

#define xGEMV(T) void  T## gemv ##_ (      TRANS,       NR, NC, T const& a, T const* A, LDA, T const* X, INCX, T const& beta, T*       Y, INCY           )  // NOLINT(bugprone-macro-parentheses,readability-identifier-length) macro arg expands to type
#define xGER( T)  void T## ger  ##_ (                   NR, NC, T const& a,                  T const* X, INCX,                T const* Y, INCY, T* A, LDA)  // NOLINT(bugprone-macro-parentheses,readability-identifier-length) macro arg expands to type
#define xGERU(T) void  T## geru ##_ (                   NR, NC, T const& a,                  T const* X, INCX,                T const* Y, INCY, T* A, LDA)  // NOLINT(bugprone-macro-parentheses,readability-identifier-length) macro arg expands to type
#define xGERC(T) void  T## gerc ##_ (                   NR, NC, T const& a,                  T const* X, INCX,                T const* Y, INCY, T* A, LDA)  // NOLINT(bugprone-macro-parentheses,readability-identifier-length) macro arg expands to type
#define xTRSV(T) void  T## trsv ##_ (UPLO, TRANS, DIAG, N,                  T const* A, LDA, T* X      , INCX                                            )  // NOLINT(bugprone-macro-parentheses,readability-identifier-length) macro arg expands to type

xGEMV(s); xGEMV(d); xGEMV(c); xGEMV(z);
xGER(s) ; xGER(d) ;
xGERU(c); xGERU(z);
xGERC(c); xGERC(z);
xTRSV(s); xTRSV(d); xTRSV(c); xTRSV(z);

#define TRANSA const char& transa
#define TRANSB const char& transb
#define NK  INTEGER nk
#define LDB INTEGER ldb
#define LDC INTEGER ldc

#define SIDE const char& side

#define xGEMM(T)     void T ##gemm ##_ (            TRANSA, TRANSB,       NR, NC, NK, T  const& a, T const* A, LDA, T const* B, LDB, T  const& b     , T const* CC, LDC)  // NOLINT(readability-identifier-length) conventional BLAS naming
#define xSYRK(T)     void T ##syrk ##_ (      UPLO, TRANSA,               NR, NK,     T  const& a, T const* A, LDA,                  T  const& b     , T*       CC, LDC)  // NOLINT(bugprone-macro-parentheses,readability-identifier-length) macro arg expands to type
#define xHERK(TT, T) void T ##herk ##_ (      UPLO, TRANSA,               NR, NK,     TT const& a, T const* A, LDA,                  TT const& b     , T*       CC, LDC)  // NOLINT(bugprone-macro-parentheses,readability-identifier-length) macro arg expands to type
#define xTRSM(T)     void T ##trsm ##_ (SIDE, UPLO, TRANSA,         DIAG, NR, NK,     T  const& a, T const* A, LDA,                  T  const* B, LDB                  )  // NOLINT(readability-identifier-length) conventional BLAS naming

xGEMM(s); xGEMM(d); xGEMM(c)   ; xGEMM(z)   ;
xSYRK(s); xSYRK(d); xSYRK(c)   ; xSYRK(z)   ;
                    xHERK(s, c); xHERK(d, z);
xTRSM(s); xTRSM(d); xTRSM(c)   ; xTRSM(z)   ;

#undef TRANS
#undef UPLO
#undef SIDE
#undef DIAG
#undef xROTG
#undef xROTMG
#undef xROT
#undef xROTM
#undef xSCAL
#undef xSWAP
#undef xCOPY
#undef xAXPY

#undef xDOT
#undef xDOTU
#undef xDOTC

#undef xNRM2
#undef xASUM
#undef IxAMAX
#undef xGEMV
#undef xGER
#undef xGERU
#undef xGERC
#undef xGEMM
#undef xHERK
#undef xTRSM

#undef s
#undef d
#undef c
#undef z
#undef C
#undef Z
#undef v
#undef INTEGER
#undef N
#undef INCX
#undef INCY
#undef TRANSA
#undef TRANSB
#undef LDA
#undef LDB
#undef LDC
}  // end of extern "C"
#endif

namespace boost::multi::blas {

// Boundary Checked value
#define BC(value) [](auto checked) {assert(checked >= std::numeric_limits<INT>::min() && checked < std::numeric_limits<INT>::max()); return checked;}(value)  /*NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay)*/

namespace core {

using s = float;
using d = double;
using c = std::complex<float>;
using z = std::complex<double>;
using v = void;

using std::enable_if_t;
using std::is_assignable;

using ::core::ssize_t;

// TODO(correaa) implement xrotg, xrotmg, xrot, xrotm

// NOLINTBEGIN(modernize-use-constraints) TODO(correaa) for C++20
template<class SX, class SY, enable_if_t<is_s<SX>{} && is_s<SY>{} && is_assignable<SY&, SX&>{},int> =0> void swap(ssize_t n, SX* x, ptrdiff_t incx, SY* y, ptrdiff_t incy) noexcept {BLAS(sswap)(n, reinterpret_cast<             float   *>(x), incx, reinterpret_cast<             float  *>(y), incy);}  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast,cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays,google-readability-casting,readability-identifier-length) // NOSONAR
template<class DX, class DY, enable_if_t<is_d<DX>{} && is_d<DY>{} && is_assignable<DY&, DX&>{},int> =0> void swap(ssize_t n, DX* x, ptrdiff_t incx, DY* y, ptrdiff_t incy) noexcept {BLAS(dswap)(n, reinterpret_cast<             double *>(x), incx, reinterpret_cast<             double *>(y), incy);}  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast,cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays,google-readability-casting,readability-identifier-length) // NOSONAR
template<class CX, class CY, enable_if_t<is_c<CX>{} && is_c<CY>{} && is_assignable<CY&, CX&>{},int> =0> void swap(ssize_t n, CX* x, ptrdiff_t incx, CY* y, ptrdiff_t incy) noexcept {BLAS(cswap)(n, reinterpret_cast<std::complex<float  >*>(x), incx, reinterpret_cast<std::complex<float >*>(y), incy);}  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast,cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays,google-readability-casting,readability-identifier-length) // NOSONAR
template<class ZX, class ZY, enable_if_t<is_z<ZX>{} && is_z<ZY>{} && is_assignable<ZY&, ZX&>{},int> =0> void swap(ssize_t n, ZX* x, ptrdiff_t incx, ZY* y, ptrdiff_t incy) noexcept {BLAS(zswap)(n, reinterpret_cast<std::complex<double>*>(x), incx, reinterpret_cast<std::complex<double>*>(y), incy);}  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast,cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays,google-readability-casting,readability-identifier-length) // NOSONAR

template<class SX, class SY, enable_if_t<is_s<SX>{} && is_s<SY>{} && is_assignable<SY&, SX&>{},int> =0> void copy(ssize_t n, SX* x, ptrdiff_t incx, SY* y, ptrdiff_t incy) {BLAS(scopy)(n, reinterpret_cast<             float    const*>(x), incx, reinterpret_cast<             float   *>(y), incy);}  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast,cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays,google-readability-casting,readability-identifier-length)
template<class DX, class DY, enable_if_t<is_d<DX>{} && is_d<DY>{} && is_assignable<DY&, DX&>{},int> =0> void copy(ssize_t n, DX* x, ptrdiff_t incx, DY* y, ptrdiff_t incy) {BLAS(dcopy)(n, reinterpret_cast<             double  const*>(x), incx, reinterpret_cast<             double *>(y), incy);}  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast,cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays,google-readability-casting,readability-identifier-length)
template<class CX, class CY, enable_if_t<is_c<CX>{} && is_c<CY>{} && is_assignable<CY&, CX&>{},int> =0> void copy(ssize_t n, CX* x, ptrdiff_t incx, CY* y, ptrdiff_t incy) {BLAS(ccopy)(n, reinterpret_cast<std::complex<float  > const*>(x), incx, reinterpret_cast<std::complex<float  >*>(y), incy);}  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast,cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays,google-readability-casting,readability-identifier-length)
template<class ZX, class ZY, enable_if_t<is_z<ZX>{} && is_z<ZY>{} && is_assignable<ZY&, ZX&>{},int> =0> void copy(ssize_t n, ZX* x, ptrdiff_t incx, ZY* y, ptrdiff_t incy) {BLAS(zcopy)(n, reinterpret_cast<std::complex<double> const*>(x), incx, reinterpret_cast<std::complex<double>*>(y), incy);}  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast,cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays,google-readability-casting,readability-identifier-length)
// NOLINTEND(modernize-use-constraints) TODO(correaa) for C++20

// TODO(correaa) : add mixed-type scal (zdscal, csscal)
// NOLINTBEGIN(modernize-use-constraints) TODO(correaa) for C++20
template<class ALPHAP, class SXP, class SX = typename std::pointer_traits<SXP>::element_type, class ALPHA = typename std::pointer_traits<ALPHAP>::element_type, enable_if_t<is_s<SX>{} && is_s<ALPHA>{} && is_assignable<SX&, decltype(*ALPHAP{}*SX{})>{}>* = nullptr> void scal(ssize_t n, ALPHAP a, SXP xp, ptrdiff_t incx) {BLAS(sscal)(n, *reinterpret_cast<             float    const*>(a), reinterpret_cast<             float   *>(xp), incx);}  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast,cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays,google-readability-casting,readability-identifier-length)
template<class ALPHAP, class DXP, class DX = typename std::pointer_traits<DXP>::element_type, class ALPHA = typename std::pointer_traits<ALPHAP>::element_type, enable_if_t<is_d<DX>{} && is_d<ALPHA>{} && is_assignable<DX&, decltype(*ALPHAP{}*DX{})>{}>* = nullptr> void scal(ssize_t n, ALPHAP a, DXP xp, ptrdiff_t incx) {BLAS(dscal)(n, *reinterpret_cast<             double  const*>(a), reinterpret_cast<             double *>(xp), incx);}  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast,cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays,google-readability-casting,readability-identifier-length)
template<class ALPHAP, class CXP, class CX = typename std::pointer_traits<CXP>::element_type, class ALPHA = typename std::pointer_traits<ALPHAP>::element_type, enable_if_t<is_c<CX>{} && is_c<ALPHA>{} && is_assignable<CX&, decltype(*ALPHAP{}*CX{})>{}>* = nullptr> void scal(ssize_t n, ALPHAP a, CXP xp, ptrdiff_t incx) {BLAS(cscal)(n, *reinterpret_cast<std::complex<float  > const*>(a), reinterpret_cast<std::complex<float  >*>(xp), incx);}  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast,cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays,google-readability-casting,readability-identifier-length)
template<class ALPHAP, class ZXP, class ZX = typename std::pointer_traits<ZXP>::element_type, class ALPHA = typename std::pointer_traits<ALPHAP>::element_type, enable_if_t<is_z<ZX>{} && is_z<ALPHA>{} && is_assignable<ZX&, decltype(*ALPHAP{}*ZX{})>{}>* = nullptr> void scal(ssize_t n, ALPHAP a, ZXP xp, ptrdiff_t incx) {BLAS(zscal)(n, *reinterpret_cast<std::complex<double> const*>(a), reinterpret_cast<std::complex<double>*>(xp), incx);}  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast,cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays,google-readability-casting,readability-identifier-length)
// NOLINTEND(modernize-use-constraints) TODO(correaa) for C++20

using std::pointer_traits;
using std::enable_if_t;
using std::is_convertible_v;

#define xaxpy(T) \
template<class ALPHA, class SXP, class SX = typename pointer_traits<SXP>::element_type, class SYP, class SY = typename pointer_traits<SYP>::element_type, enable_if_t<  /* NOLINT(modernize-use-constraints) */ \
	is_##T<ALPHA>{} && is_##T<SX>{} && is_##T<SY>{} && is_assignable<SY&, decltype(ALPHA{}*SX{})>{} \
	&& is_convertible_v<SXP, SX*> && is_convertible_v<SYP, SY*> \
, int> =0> \
void axpy(size_t n, ALPHA const* a, SXP x, size_t incx, SYP y, size_t incy) {BLAS(T##axpy)(n, reinterpret_cast<T const *>(a), reinterpret_cast<T const*>(static_cast<SX*>(x)), incx, reinterpret_cast<T*>(static_cast<SY*>(y)), incy);}  /*NOLINT(readability-identifier-length,cppcoreguidelines-pro-type-reinterpret-cast,bugprone-macro-parentheses) conventional BLAS name*/

xaxpy(s)       xaxpy(d)       xaxpy(c)       xaxpy(z)
#undef  xaxpy
}  // end namespace core

// #undef xrotg
// #undef xrot
// #undef xswap
// #undef xscal
// #undef xcopy
// #undef xaxpy
// #undef xdot

#ifndef CBLAS_H

namespace core {

using std::enable_if_t;
using std::is_assignable;

// NOLINTBEGIN(modernize-use-constraints) TODO(correaa) for C++20
template<class XP, class X = typename std::pointer_traits<XP*>::element_type, class YP, class Y = typename std::pointer_traits<YP*>::element_type, class RP, class R = typename std::pointer_traits<RP*>::element_type, enable_if_t<is_s<X>{} && is_s<Y>{} && is_assignable<R&, decltype(0.0F+X{}*Y{}+X{}*Y{})>{}, int> =0> void dot (ssize_t n, XP* xp, ptrdiff_t incx, YP* yp, ptrdiff_t incy, RP* rp) {
	// Apple Accelerate BLAS is known to have bugs in single precission function
	// `sdot` and `smrm2`, be careful:
	// https://stackoverflow.com/a/77017238/225186,
	// https://fortran-lang.discourse.group/t/how-many-blas-libraries-have-this-error/4454/23,
	// https://forums.developer.apple.com/forums/thread/717757
#ifdef MULTI_BLAS_USE_SDOT  // disable workararound for Apple Accelerate framework bug
	auto const rr = BLAS(sdot )(n, reinterpret_cast<s const*>(static_cast<X*>(xp)), incx, reinterpret_cast<s const*>(static_cast<Y*>(yp)), incy); std::memcpy(reinterpret_cast<float  *    >(static_cast<R*>(r)), &rr, sizeof(rr)); static_assert(sizeof(rr)==sizeof(*r));
#else
	BLAS(sgemv)('N', 1, n, 1.0F, reinterpret_cast<s const*>(static_cast<X*>(xp)), incx, reinterpret_cast<s const*>(static_cast<Y*>(yp)), incy, 0.0F, reinterpret_cast<s*>(static_cast<R*>(rp)), 1);  // NOLINT(readability-suspicious-call-argument,cppcoreguidelines-pro-type-reinterpret-cast) 
#endif
}
template<class XP, class X = typename std::pointer_traits<XP*>::element_type, class YP, class Y = typename std::pointer_traits<YP*>::element_type, class RP, class R = typename std::pointer_traits<RP*>::element_type, enable_if_t<is_d<X>{} && is_d<Y>{} && is_assignable<R&, decltype(0.0 +X{}*Y{}+X{}*Y{})>{}, int> =0> void dot (ssize_t n, XP* xp, ptrdiff_t incx, YP* yp, ptrdiff_t incy, RP* r) {auto const rr = BLAS(ddot )(n, reinterpret_cast<d const*>(static_cast<X*>(xp)), incx, reinterpret_cast<d const*>(static_cast<Y*>(yp)), incy); std::memcpy(reinterpret_cast<double *    >(static_cast<R*>(r)), &rr, sizeof(rr)); static_assert(sizeof(rr)==sizeof(*r));} // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast,cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays,google-readability-casting,readability-identifier-length)  // NOSONAR
// NOLINTEND(modernize-use-constraints) TODO(correaa) for C++20

// PGI/NVC++ compiler uses a blas version that needs -DRETURN_BY_STACK
// NOLINTBEGIN(modernize-use-constraints) TODO(correaa) for C++20
#if defined(BLAS_DOT_RETURNS_VOID)
template<class XP, class X = typename std::pointer_traits<XP>::element_type, class YP, class Y = typename std::pointer_traits<YP>::element_type, class RP, class R = typename std::pointer_traits<RP>::element_type, enable_if_t<is_c<X>{} && is_c<Y>{} && is_assignable<R&, decltype(0.0F+X{}*Y{}+X{}*Y{})>{}, int> =0> void dotu(ssize_t n, XP xp, ptrdiff_t incx, YP yp, ptrdiff_t incy, RP rp) {
	[[maybe_unused]] static bool const check = []{
		std::array<std::complex<float>, 3> const v1 = {std::complex<float>{1.0F, 2.0F}, std::complex<float>{3.0F,  4.0F}, std::complex<float>{ 5.0F,  6.0F}};
		std::array<std::complex<float>, 3> const v2 = {std::complex<float>{7.0F, 8.0F}, std::complex<float>{9.0F, 10.0F}, std::complex<float>{11.0F, 12.0F}};
		Complex_float rr{-1.0F, -2.0F};
		BLAS(cdotu)(&rr, 3, v1.data(), 1, v2.data(), 1);
		if( std::abs(rr.real - std::real(v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2])) > 1.0e-8 ) { throw std::logic_error("[real] cdotu should be configured as non-void returing"); }
		if( std::abs(rr.imag - std::imag(v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2])) > 1.0e-8 ) { throw std::logic_error("[imag] cdotu should be configured as non-void returing"); }
		return true;
	}();
	// BLAS(cdotu)(reinterpret_cast<Complex_float *>(rp), n, reinterpret_cast<c const*>(static_cast<X*>(xp)), incx, reinterpret_cast<c const*>(static_cast<Y*>(yp)), incy);  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast)
	BLAS(cgemv)('N', 1, n, std::complex<float>{1.0F, 0.0F}, reinterpret_cast<c const*>(static_cast<X*>(xp)), incx, reinterpret_cast<c const*>(static_cast<Y*>(yp)), incy, std::complex<float>{0.0F, 0.0F}, reinterpret_cast<c*>(static_cast<R*>(rp)), 1);  // NOLINT(readability-suspicious-call-argument,cppcoreguidelines-pro-type-reinterpret-cast)
}                                                                                                                          // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast,google-readability-casting) : adapt types
template<class XP, class X = typename std::pointer_traits<XP>::element_type, class YP, class Y = typename std::pointer_traits<YP>::element_type, class RP, class R = typename std::pointer_traits<RP>::element_type, enable_if_t<is_z<X>{} && is_z<Y>{} && is_assignable<R&, decltype(0.0 +X{}*Y{}+X{}*Y{})>{}, int> =0> void dotu(ssize_t n, XP xp, ptrdiff_t incx, YP yp, ptrdiff_t incy, RP rp) {                BLAS(zdotu)(reinterpret_cast<Complex_double*>(rp), n, reinterpret_cast<z const*>(static_cast<X*>(xp)), incx, reinterpret_cast<z const*>(static_cast<Y*>(yp)), incy);}                                                                                                                          // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast,google-readability-casting) : adapt types

template<class XP, class X = typename std::pointer_traits<XP>::element_type, class YP, class Y = typename std::pointer_traits<YP>::element_type, class RP, class R = typename std::pointer_traits<RP>::element_type, enable_if_t<is_c<X>{} && is_c<Y>{} && is_assignable<R&, decltype(0.0F+X{}*Y{}+X{}*Y{})>{}, int> =0> void dotc(ssize_t n, XP xp, ptrdiff_t incx, YP yp, ptrdiff_t incy, RP rp) {
	std::clog << "using cdotc void\n";
	BLAS(cdotc)(reinterpret_cast<Complex_float *>(rp), n, reinterpret_cast<c const*>(static_cast<X*>(xp)), incx, reinterpret_cast<c const*>(static_cast<Y*>(yp)), incy);  // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast)
	// BLAS(cgemv)('C', n, 1, std::complex<float>{1.0F, 0.0F}, reinterpret_cast<c const*>(static_cast<X*>(xp)), incx, reinterpret_cast<c const*>(static_cast<Y*>(yp)), incy, std::complex<float>{0.0F, 0.0F}, reinterpret_cast<c*>(static_cast<R*>(rp)), 1);  // NOLINT(readability-suspicious-call-argument)
}                                                                                                                          // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast,google-readability-casting) : adapt types
template<class XP, class X = typename std::pointer_traits<XP>::element_type, class YP, class Y = typename std::pointer_traits<YP>::element_type, class RP, class R = typename std::pointer_traits<RP>::element_type, enable_if_t<is_z<X>{} && is_z<Y>{} && is_assignable<R&, decltype(0.0 +X{}*Y{}+X{}*Y{})>{}, int> =0> void dotc(ssize_t n, XP xp, ptrdiff_t incx, YP yp, ptrdiff_t incy, RP rp) {                BLAS(zdotc)(reinterpret_cast<Complex_double*>(rp), n, reinterpret_cast<z const*>(static_cast<X*>(xp)), incx, reinterpret_cast<z const*>(static_cast<Y*>(yp)), incy);}                                                                                                                          // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast,google-readability-casting) : adapt types
#else
// NOLINTBEGIN(cppcoreguidelines-pro-type-reinterpret-cast,cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays,google-readability-casting,readability-identifier-length)
// TODO(correaa) implement workaround for bug in Apple Accelerate BLAS ? https://stackoverflow.com/a/77017238/225186
template<class XP, class X = typename std::pointer_traits<XP>::element_type, class YP, class Y = typename std::pointer_traits<YP>::element_type, class RP, class R = typename std::pointer_traits<RP>::element_type, enable_if_t<is_c<X>{} && is_c<Y>{} && is_assignable<R&, decltype(/*0.0F+*/ X{}*Y{}+X{}*Y{})>{}, int> =0> void dotu(ssize_t n, XP xp, ptrdiff_t incx, YP yp, ptrdiff_t incy, RP rp) {
	// [[maybe_unused]] static bool const use_cdotu = []{
	//  std::array<std::complex<float>, 3> const v1 = {std::complex<float>{1.0F, 2.0F}, std::complex<float>{3.0F,  4.0F}, std::complex<float>{ 5.0F,  6.0F}};
	//  std::array<std::complex<float>, 3> const v2 = {std::complex<float>{7.0F, 8.0F}, std::complex<float>{9.0F, 10.0F}, std::complex<float>{11.0F, 12.0F}};

	//  Complex_float rr{-1.0F, -2.0F};
	//  rr = BLAS(cdotu)(3, v1.data(), 1, v2.data(), 1);

	//  if( !(std::abs(rr.real - std::real(v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2])) < 1.0e-8)
	//    || !(std::abs(rr.imag - std::imag(v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2])) < 1.0e-8) ) {
	//    std::clog
	//      << "multi::blas setup warning: when using cdotu that returns non-void,\n"
	//      << "cdotu returned (" << rr.real << ", " << rr.imag << ", it should return " << v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2] << '\n'
	//      << "This problem appears with BLAS and OpenBLAS 32bit.\n"
	//      << "... falling back to cgemv\n";
	//    {
	//      std::complex<float> gemv_rr{-12.345F, -54.321F};
	//      BLAS(cgemv)('N', 1, v1.size(), std::complex<float>{1.0F, 0.0F}, v1.data(), 1, v2.data(), 1, std::complex<float>{0.0F, 0.0F}, &gemv_rr, 1);
	//      std::clog << "cgemv gives " << gemv_rr << ", it should give " << v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2] << '\n';
			
	//      if( !(std::abs(gemv_rr - (v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2])) < 1.0e-8) ) {
	//        std::clog << "gemv also failed" << '\n';
	//      }
	//      return false;  // dot not use cdotu
	//    }
	//  }
	//  return true;  // use cdotu
	// }();
	// if(use_cdotu) {
	//  Complex_float const rr = BLAS(cdotu)(                                      n, reinterpret_cast<c const*>(static_cast<X*>(x)), incx, reinterpret_cast<c const*>(static_cast<Y*>(y)), incy); std::memcpy(reinterpret_cast<std::array<float , 2>*>(static_cast<R*>(rp))->data(), &rr, sizeof(rr)); static_assert(sizeof(rr)==sizeof(*rp));
	// } else {
		BLAS(cgemv)('N', 1, n, std::complex<float>{1.0F, 0.0F}, reinterpret_cast<c const*>(static_cast<X*>(xp)), incx, reinterpret_cast<c const*>(static_cast<Y*>(yp)), incy, std::complex<float>{0.0F, 0.0F}, reinterpret_cast<c*>(static_cast<R*>(rp)), 1);  // NOLINT(readability-suspicious-call-argument)
	// }
}
template<class XP, class X = typename std::pointer_traits<XP>::element_type, class YP, class Y = typename std::pointer_traits<YP>::element_type, class RP, class R = typename std::pointer_traits<RP>::element_type, enable_if_t<is_z<X>{} && is_z<Y>{} && is_assignable<R&, decltype(/*0.0 +*/ X{}*Y{}+X{}*Y{})>{}, int> =0> void dotu(ssize_t n, XP xp, ptrdiff_t incx, YP yp, ptrdiff_t incy, RP rp) {
	// auto const rr = BLAS(zdotu)(                                      n, reinterpret_cast<z const*>(static_cast<X*>(xp)), incx, reinterpret_cast<z const*>(static_cast<Y*>(yp)), incy); std::memcpy(reinterpret_cast<std::array<double, 2>*>(static_cast<R*>(rp))->data(), &rr, sizeof(rr)); static_assert(sizeof(rr)==sizeof(*rp));
	BLAS(zgemv)('N', 1, n, std::complex<double>{1.0, 0.0}, reinterpret_cast<z const*>(static_cast<X*>(xp)), incx, reinterpret_cast<z const*>(static_cast<Y*>(yp)), incy, std::complex<double>{0.0, 0.0}, reinterpret_cast<z*>(static_cast<R*>(rp)), 1);  // NOLINT(readability-suspicious-call-argument)
}

template<class XP, class X = typename std::pointer_traits<XP>::element_type, class YP, class Y = typename std::pointer_traits<YP>::element_type, class RP, class R = typename std::pointer_traits<RP>::element_type, enable_if_t<is_c<X>{} && is_c<Y>{} && is_assignable<R&, decltype(/*0.0F+*/ X{}*Y{}+X{}*Y{})>{}, int> =0> void dotc(ssize_t n, XP xp, ptrdiff_t incx, YP yp, ptrdiff_t incy, RP rp) {
	// std::clog << "using cdotc non void\n";
	// c
	auto const rr = BLAS(cdotc)(                                      n, reinterpret_cast<c const*>(static_cast<X*>(xp)), incx, reinterpret_cast<c const*>(static_cast<Y*>(yp)), incy); std::memcpy(reinterpret_cast<std::array<float ,2>*>(static_cast<R*>(rp))->data(), &rr, sizeof(rr)); static_assert(sizeof(rr)==sizeof(*rp));
	// BLAS(cgemv)('N', 1, n, std::complex<float>{1.0F, 0.0F}, reinterpret_cast<c const*>(static_cast<X*>(xp)), incx, reinterpret_cast<c const*>(static_cast<Y*>(yp)), incy, std::complex<float>{0.0F, 0.0F}, reinterpret_cast<c*>(static_cast<R*>(rp)), 1);  // NOLINT(readability-suspicious-call-argument)
}

template<class XP, class X = typename std::pointer_traits<XP>::element_type, class YP, class Y = typename std::pointer_traits<YP>::element_type, class RP, class R = typename std::pointer_traits<RP>::element_type, enable_if_t<is_z<X>{} && is_z<Y>{} && is_assignable<R&, decltype(/*0.0 +*/ X{}*Y{}+X{}*Y{})>{}, int> =0> void dotc(ssize_t n, XP x, ptrdiff_t incx, YP y, ptrdiff_t incy, RP r) {auto const rr = BLAS(zdotc)(                                      n, reinterpret_cast<z const*>(static_cast<X*>(x)), incx, reinterpret_cast<z const*>(static_cast<Y*>(y)), incy); std::memcpy(reinterpret_cast<std::array<double,2>*>(static_cast<R*>(r))->data(), &rr, sizeof(rr)); static_assert(sizeof(rr)==sizeof(*r));}  // NOSONAR
// NOLINTEND(cppcoreguidelines-pro-type-reinterpret-cast,cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays,google-readability-casting,readability-identifier-length)
#endif
// NOLINTEND(modernize-use-constraints)

} // end namespace core
#else
// TODO(correaa) : make cblas version
#define xdotu(T) template<class S> v dotu(S n, add_const_ptr_t<T> x, S incx, add_const_ptr_t<T> y, S incy, add_ptr_t<T> r){BLAS(T##dotu_sub)(BC(n), x, BC(incx), y, BC(incy), r);}
#define xdotc(T) template<class S> v dotc(S n, add_const_ptr_t<T> x, S incx, add_const_ptr_t<T> y, S incy, add_ptr_t<T> r){BLAS(T##dotc_sub)(BC(n), x, BC(incx), y, BC(incy), r);}

namespace core {
	xdotu(c) xdotu(z)
	xdotc(c) xdotc(z)
}

#undef xdotu
#undef xdotc
#endif

namespace core {
	template<class S> auto dot(S n, s const& b, s const* x, S incx, s const* y, S incy) -> s {return BLAS(sdsdot)(BC(n), b, x, BC(incx), y, BC(incy));}  // NOLINT(readability-identifier-length) conventional BLAS name
} // end namespace core

#define xasum(T, TT)    template<class S> auto asum (S n, T const* x, S incx){return BLAS(TT##asum  )(BC(n), x, BC(incx))    ;}  // NOLINT(readability-identifier-length) conventional BLAS name
#define ixamax(T)       template<class S> auto iamax(S n, T const* x, S incx){return BLAS(i##T##amax)(BC(n), x, BC(incx)) - 1;}  // NOLINT(readability-identifier-length) conventional BLAS name

namespace core {

// NOLINTBEGIN(modernize-use-constraints) TODO(correaa) for C++20
template<class XP, class X = typename std::pointer_traits<XP>::element_type, class RP, class R = typename std::pointer_traits<RP>::element_type, enable_if_t<is_s<X>{} && is_s<R>{} && std::is_assignable<R&, decltype(X{})>{}           , int> =0> void asum(ssize_t n, XP x, ptrdiff_t incx, RP r) {auto rr = BLAS(sasum) (n, reinterpret_cast<s const*>(static_cast<X*>(x)), incx); std::memcpy(reinterpret_cast<s*>(static_cast<R*>(r)), &rr, sizeof(s));} // NOLINT(google-readability-casting,readability-identifier-length,cppcoreguidelines-pro-type-reinterpret-cast)
template<class XP, class X = typename std::pointer_traits<XP>::element_type, class RP, class R = typename std::pointer_traits<RP>::element_type, enable_if_t<is_d<X>{} && is_d<R>{} && std::is_assignable<R&, decltype(X{})>{}           , int> =0> void asum(ssize_t n, XP x, ptrdiff_t incx, RP r) {auto rr = BLAS(dasum) (n, reinterpret_cast<d const*>(static_cast<X*>(x)), incx); std::memcpy(reinterpret_cast<s*>(static_cast<R*>(r)), &rr, sizeof(d));} // NOLINT(google-readability-casting,readability-identifier-length,cppcoreguidelines-pro-type-reinterpret-cast)

template<class XP, class X = typename std::pointer_traits<XP>::element_type, class RP, class R = typename std::pointer_traits<RP>::element_type, enable_if_t<is_c<X>{} && is_s<R>{} && std::is_assignable<R&, decltype(norm(X{}))>{}, int> =0> void asum(ssize_t n, XP x, ptrdiff_t incx, RP r) {auto rr   = BLAS(scasum)(n, reinterpret_cast<c const*>(static_cast<X*>(x)), incx); std::memcpy(reinterpret_cast<s*>(static_cast<R*>(r)), &rr, sizeof(s));} // NOLINT(google-readability-casting,readability-identifier-length,cppcoreguidelines-pro-type-reinterpret-cast)
template<class XP, class X = typename std::pointer_traits<XP>::element_type, class RP, class R = typename std::pointer_traits<RP>::element_type, enable_if_t<is_z<X>{} && is_d<R>{} && std::is_assignable<R&, decltype(norm(X{}))>{}, int> =0> void asum(ssize_t n, XP x, ptrdiff_t incx, RP r) {auto rr = BLAS(dzasum)(n, reinterpret_cast<z const*>(static_cast<X*>(x)), incx); std::memcpy(reinterpret_cast<s*>(static_cast<R*>(r)), &rr, sizeof(d));} // NOLINT(google-readability-casting,readability-identifier-length,cppcoreguidelines-pro-type-reinterpret-cast)

// TODO(correaa) implement workaround for bug in Apple Accelerate BLAS ? https://stackoverflow.com/a/77017238/225186
template<class XP, class X = typename std::pointer_traits<XP>::element_type, class RP, class R = typename std::pointer_traits<RP>::element_type, enable_if_t<is_s<X>{} && is_s<R>{} && std::is_assignable<R&, decltype(X{})>{}           , int> =0> void nrm2(ssize_t n, XP x, ptrdiff_t incx, RP r) {auto rr = BLAS(snrm2) (n, reinterpret_cast<s const*>(static_cast<X*>(x)), incx); std::memcpy(reinterpret_cast<s*>(static_cast<R*>(r)), &rr, sizeof(s));} // NOLINT(google-readability-casting,readability-identifier-length,cppcoreguidelines-pro-type-reinterpret-cast)
template<class XP, class X = typename std::pointer_traits<XP>::element_type, class RP, class R = typename std::pointer_traits<RP>::element_type, enable_if_t<is_d<X>{} && is_d<R>{} && std::is_assignable<R&, decltype(X{})>{}           , int> =0> void nrm2(ssize_t n, XP x, ptrdiff_t incx, RP r) {auto rr = BLAS(dnrm2) (n, reinterpret_cast<d const*>(static_cast<X*>(x)), incx); std::memcpy(reinterpret_cast<s*>(static_cast<R*>(r)), &rr, sizeof(d));} // NOLINT(google-readability-casting,readability-identifier-length,cppcoreguidelines-pro-type-reinterpret-cast)

template<class XP, class X = typename std::pointer_traits<XP>::element_type, class RP, class R = typename std::pointer_traits<RP>::element_type, enable_if_t<is_c<X>{} && is_s<R>{} && std::is_assignable<R&, decltype(norm(X{}))>{}, int> =0> void nrm2(ssize_t n, XP x, ptrdiff_t incx, RP r) {auto rr = BLAS(scnrm2)(n, reinterpret_cast<c const*>(static_cast<X*>(x)), incx); std::memcpy(reinterpret_cast<s*>(static_cast<R*>(r)), &rr, sizeof(s));} // NOLINT(google-readability-casting,readability-identifier-length,cppcoreguidelines-pro-type-reinterpret-cast)
template<class XP, class X = typename std::pointer_traits<XP>::element_type, class RP, class R = typename std::pointer_traits<RP>::element_type, enable_if_t<is_z<X>{} && is_d<R>{} && std::is_assignable<R&, decltype(norm(X{}))>{}, int> =0> void nrm2(ssize_t n, XP x, ptrdiff_t incx, RP r) {auto rr = BLAS(dznrm2)(n, reinterpret_cast<z const*>(static_cast<X*>(x)), incx); std::memcpy(reinterpret_cast<s*>(static_cast<R*>(r)), &rr, sizeof(d));} // NOLINT(google-readability-casting,readability-identifier-length,cppcoreguidelines-pro-type-reinterpret-cast)
// NOLINTEND(modernize-use-constraints) TODO(correaa) for C++20

	ixamax(s)      ixamax(d)       ixamax(c)       ixamax(z)
} // end namespace core

#undef xasum
#undef ixamax

///////////////////////////////////////////////////////////////////////////////
// LEVEL2
//#define xgemv(T) template<class C, class S> v gemv(C trans, S m, S n, T const& a, T const* A, S lda, T const* X, S incx, T beta, T*       Y, S incy             ){BLAS(T##gemv)(trans, BC(m), BC(n), a, A, BC(lda), X, BC(incx), beta, Y, BC(incy)            );}
//#define xger(T)  template<         class S> v ger (         S m, S n, T const& a,                    T const* X, S incx,         T const* Y, S incy, T* A, S lda){BLAS(T##ger )(       BC(m), BC(n), a,             X, BC(incx),       Y, BC(incy), A, BC(lda));}
//                 template<         class S> v ger (         S m, S n, c const& a,                    c const* X, S incx,         c const* Y, S incy, c* A, S lda){BLAS(cgeru  )(       BC(m), BC(n), a,             X, BC(incx),       Y, BC(incy), A, BC(lda));}
//                template<         class S> v ger (         S m, S n, z const& a,                    z const* X, S incx,         z const* Y, S incy, z* A, S lda){BLAS(zgeru  )(       BC(m), BC(n), a,             X, BC(incx),       Y, BC(incy), A, BC(lda));}
//#define xgeru(T) template<         class S> v geru(         S m, S n, T const& a,                    T const* X, S incx,         T const* Y, S incy, T* A, S lda){BLAS(T##geru)(       BC(m), BC(n), a,             X, BC(incx),       Y, BC(incy), A, BC(lda));}
//#define xgerc(T) template<         class S> v gerc(         S m, S n, T const& a,                    T const* X, S incx,         T const* Y, S incy, T* A, S lda){BLAS(T##gerc)(       BC(m), BC(n), a,             X, BC(incx),       Y, BC(incy), A, BC(lda));}

namespace core {

using std::enable_if_t;
using std::is_assignable;

using size_t = INT;
using ssize_t = std::make_signed_t<size_t>;

// NOLINTBEGIN(modernize-use-constraints) TODO(correaa) for C++20
template<class A, class M, class X, class B, class Y, enable_if_t<is_s<M>{} && is_s<X>{} && is_s<Y>{} && is_assignable<Y&, decltype(A{}*M{}*X{}+B{}*Y{})>{}, int> =0> void gemv(char trans, size_t m, size_t n, A const* a, M* ma, size_t lda, X* x, size_t incx, B const* b, Y* y, size_t incy) {BLAS(sgemv)(trans, m, n, *a, reinterpret_cast<s const*>(ma), lda, reinterpret_cast<s const*>(x), incx, *b, reinterpret_cast<s*>(y), incy);}  // NOLINT(google-readability-casting,readability-identifier-length,cppcoreguidelines-pro-type-reinterpret-cast) // NOSONAR wrapped func has 11 params
template<class A, class M, class X, class B, class Y, enable_if_t<is_d<M>{} && is_d<X>{} && is_d<Y>{} && is_assignable<Y&, decltype(A{}*M{}*X{}+B{}*Y{})>{}, int> =0> void gemv(char trans, size_t m, size_t n, A const* a, M* ma, size_t lda, X* x, size_t incx, B const* b, Y* y, size_t incy) {BLAS(dgemv)(trans, m, n, *a, reinterpret_cast<d const*>(ma), lda, reinterpret_cast<d const*>(x), incx, *b, reinterpret_cast<d*>(y), incy);}  // NOLINT(google-readability-casting,readability-identifier-length,cppcoreguidelines-pro-type-reinterpret-cast) // NOSONAR wrapped func has 11 params
template<class A, class M, class X, class B, class Y, enable_if_t<is_c<M>{} && is_c<X>{} && is_c<Y>{} && is_assignable<Y&, decltype(A{}*M{}*X{}+B{}*Y{})>{}, int> =0> void gemv(char trans, size_t m, size_t n, A const* a, M* ma, size_t lda, X* x, size_t incx, B const* b, Y* y, size_t incy) {BLAS(cgemv)(trans, m, n, *a, reinterpret_cast<c const*>(ma), lda, reinterpret_cast<c const*>(x), incx, *b, reinterpret_cast<c*>(y), incy);}  // NOLINT(google-readability-casting,readability-identifier-length,cppcoreguidelines-pro-type-reinterpret-cast) // NOSONAR wrapped func has 11 params
template<class A, class M, class X, class B, class Y, enable_if_t<is_z<M>{} && is_z<X>{} && is_z<Y>{} && is_assignable<Y&, decltype(std::declval<A const&>()*std::declval<M const&>()*std::declval<X const&>()+std::declval<B const&>()*std::declval<Y const&>())>{}, int> =0> void gemv(char trans, size_t m, size_t n, A const* a, M* ma, size_t lda, X* x, size_t incx, B const* b, Y* y, size_t incy) {  // NOLINT(google-readability-casting,readability-identifier-length,cppcoreguidelines-pro-type-reinterpret-cast)  //NOSONAR wrapped func has 11 params
	BLAS(zgemv)(trans, m, n, *a, reinterpret_cast<z const*>(ma), lda, reinterpret_cast<z const*>(x), incx, *b, reinterpret_cast<z*>(y), incy);  // NOLINT(fuchsia-default-arguments-calls,google-readability-casting,readability-identifier-length,cppcoreguidelines-pro-type-reinterpret-cast)  // NOSONAR
}
// NOLINTEND(modernize-use-constraints) TODO(correaa) for C++20

// TODO(correaa) implement get, geru, gerc

using s = float;
using d = double;
using c = std::complex<float>;
using z = std::complex<double>;
using v = void;

template<class T>
struct blas2 {};

template<> struct blas2<s> {template<class... As> static v    trsv(As... args)                                   {BLAS(strsv)(args...);}};
template<> struct blas2<d> {template<class... As> static v    trsv(As... args)                                   {BLAS(dtrsv)(args...);}};
template<> struct blas2<c> {template<class... As> static v    trsv(As... args)                                   {BLAS(ctrsv)(args...);}};
template<> struct blas2<z> {template<class... As> static auto trsv(As... args) -> decltype(BLAS(ztrsv)(args...)) {BLAS(ztrsv)(args...);}};

}  // end namespace core

///////////////////////////////////////////////////////////////////////////////
// LEVEL 3

namespace core {

using std::is_convertible_v;
using std::pointer_traits;
using std::enable_if_t;
using std::max;

#define xsyrk(T) \
template<class UL, class C, class S, class ALPHA, class AAP, class AA = typename pointer_traits<AAP>::element_type, class BETA, class CCP, class CC = typename pointer_traits<CCP>::element_type, \
enable_if_t<  /* NOLINT(modernize-use-constraints) TODO(correaa) for C++20 */                                                                                                                                                                                    \
	is_##T<AA>{} && is_##T<CC>{} && is_assignable<CC&, decltype(ALPHA{}*AA{}*AA{})>{} &&                                                                                                       \
	is_convertible_v<AAP, AA*> && is_convertible_v<CCP, CC*>                                                                                                                                     \
, int> =0>                                                                                                                                                                                        \
v syrk(        UL uplo, C transA,             S n, S k, ALPHA const* alpha, AAP aa, S lda,             BETA const* beta, CCP cc, S ldc)  /*NOLINT(bugprone-easily-swappable-parameters,readability-identifier-length)*/      \
/*=delete;*/                                                                                                                                                                                      \
{                                                                                                                                                                                                 \
	if(transA == 'N' || transA == 'n') { BOOST_MULTI_ASSERT1( lda >= max(S{1}, n) ); }                                                                                                                     \
	if(transA != 'N' && transA != 'n') { BOOST_MULTI_ASSERT1( lda >= max(S{1}, k) ); }                                                                                                                     \
	BOOST_MULTI_ASSERT1( ldc >= max(S{1}, n) );                                                                                                                                                           \
	BLAS(T##syrk)(      uplo, transA,            BC(n), BC(k), *reinterpret_cast<T const*>(alpha), aa, BC(lda),        *reinterpret_cast<T const*>(beta), cc, BC(ldc));  /*NOLINT(cppcoreguidelines-pro-type-reinterpret-cast)*/                                                               \
}                                                                                                                                                                                                 \

#define xherk(T) \
template<class UL, class C, class S, class ALPHA, class AAP, class AA = typename pointer_traits<AAP>::element_type, class BETA, class CCP, class CC = typename pointer_traits<CCP>::element_type, class Real = typename T::value_type, \
enable_if_t<  /* NOLINT(modernize-use-constraints) TODO(correaa) for C++20 */                                                                                                                                                                                                                         \
	is_##T<AA>{} && is_##T<CC>{} && is_assignable<CC&, decltype(ALPHA{}*AA{}*AA{})>{} &&                                                                                                                                            \
	is_convertible_v<AAP, AA*> && is_convertible_v<CCP, CC*>                                                                                                                                                                          \
, int> =0>                                                                                                                                                                                                                             \
v herk(        UL uplo, C transA,             S n, S k, ALPHA const* alpha, AAP aa, S lda,             BETA const* beta, CCP cc, S ldc)  /*NOLINT(bugprone-easily-swappable-parameters,readability-identifier-length)*/                \
/*=delete;*/                                                                                                                                                                                                                           \
{                                                                                                                                                                                                                                      \
	if(transA == 'N' ||  transA == 'n') {BOOST_MULTI_ASSERT1( lda >= max(S{1}, n) );}  /* NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay)*/                                                                   \
	if(transA != 'N' && transA != 'n') {BOOST_MULTI_ASSERT1( lda >= max(S{1}, k) );}                                                                                                                                                          \
	BOOST_MULTI_ASSERT1( ldc >= max(S{1}, n) );  /* NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay)*/                                                                                                         \
	/*BOOST_MULTI_MARK_SCOPE("cpu_herk");*/                                                                                                                                                                                                      \
	BLAS(T##herk)(      uplo, transA,            BC(n), BC(k), *reinterpret_cast<Real const*>(alpha), aa, BC(lda),        *reinterpret_cast<Real const*>(beta), cc, BC(ldc));  /*NOLINT(cppcoreguidelines-pro-type-reinterpret-cast)*/                                                                                            \
}                                                                                                                                                                                                                                      \

#define xgemm(T) \
template<class ALPHA, class AAP, class AA = typename pointer_traits<AAP>::element_type, class BBP, class BB = typename pointer_traits<BBP>::element_type, class BETA, class CCP, class CC = typename pointer_traits<CCP>::element_type, \
enable_if_t<  /* NOLINT(modernize-use-constraints) TODO(correaa) for C++20 */                                                                                                                                                                                                                            \
	is_##T<AA>{} && is_##T<BB>{} && is_##T<CC>{} &&                                                                                                                           \
	is_convertible_v<AAP, AA*> && is_convertible_v<BBP, BB*> && is_convertible_v<CCP, CC*>                                                                                                                                            \
, int> =0>                                                                                                                                                                                                                              \
v gemm(char transA, char transB, ssize_t m, ssize_t n, ssize_t k, ALPHA const* alpha, AAP aa, ssize_t lda, BBP bb, ssize_t ldb, BETA const* beta, CCP cc, ssize_t ldc) {  /*NOLINT(bugprone-easily-swappable-parameters)*/              \
	using std::max;                                                                                                                                                                                                                     \
	BOOST_MULTI_ASSERT1((transA != 'N') || (lda >= max(ssize_t{1}, m)));                                                                                                                                                                               \
	BOOST_MULTI_ASSERT1((transA == 'N') || (lda >= max(ssize_t{1}, k)));                                                                                                                                                                               \
	BOOST_MULTI_ASSERT1((transB != 'N') || (ldb >= max(ssize_t{1}, k)));                                                                                                                                                                               \
	BOOST_MULTI_ASSERT1((transB == 'N') || (ldb >= max(ssize_t{1}, n)));                                                                                                                                                                               \
\
	BOOST_MULTI_ASSERT1( aa != cc );                                                                                                                                                                                                          \
	BOOST_MULTI_ASSERT1( bb != cc );                                                                                                                                                                                                          \
\
	if(!( ldc >= max(ssize_t{1}, m) )) {throw std::logic_error("failed 'ldc >= max(1, m)' with ldc = "+ std::to_string(ldc) +" and m = "+ std::to_string(m));}                                                                               \
	if(*beta != 0.0) {BOOST_MULTI_ASSERT1((is_assignable<CC&, decltype(std::declval<ALPHA>()*std::declval<AA>()*std::declval<BB>() + std::declval<BETA>()*std::declval<CC>())> {}));}                                                          \
	BLAS(T##gemm)(transA, transB, BC(m), BC(n), BC(k), *reinterpret_cast<T const*>(alpha), reinterpret_cast<T const*>(static_cast<AA*>(aa)), BC(lda), reinterpret_cast<T const*>(static_cast<BB*>(bb)), BC(ldb), *reinterpret_cast<T const*>(beta), reinterpret_cast<T*>(static_cast<CC*>(cc)) /*NOLINT(cppcoreguidelines-pro-type-reinterpret-cast,bugprone-macro-parentheses)*/ /*TODO(correaa) check constness*/, BC(ldc)); \
}                                                                                                                                                                                                                        \

// NOLINTNEXTLINE(readability-identifier-length) conventional BLAS name
xgemm(s) xgemm(d) xgemm(c) xgemm(z)  // NOLINT(modernize-use-constraints,readability-function-cognitive-complexity) : 36 of 25
#undef xgemm

#define xtrsm(T) \
template<class ALPHA, class AAP, class AA = typename pointer_traits<AAP>::element_type, class BBP, class BB = typename pointer_traits<BBP>::element_type,                                                         \
enable_if_t<  /* NOLINT(modernize-use-constraints) TODO(correaa) for C++20 */                                                                                                                                                                                                      \
	is_##T<AA>{} && is_##T<BB>{} && is_assignable<BB&, decltype(AA{}*BB{}/ALPHA{})>{} && is_assignable<BB&, decltype(ALPHA{}*BB{}/AA{})>{} &&                                                                 \
	is_convertible_v<AAP, AA*> && is_convertible_v<BBP, BB*>                                                                                                                                                     \
,int> =0>                                                                                                                                                                                                         \
v trsm(char side, char uplo, char transA, char diag, ssize_t m, ssize_t n, ALPHA alpha, AAP aa, ssize_t lda, BBP bb, ssize_t ldb) { /*NOLINT(bugprone-easily-swappable-parameters,readability-identifier-length)*/  \
	assert( side   == 'L' || side   == 'R' );                   /* NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay)*/                                                             \
	assert( uplo   == 'U' || uplo   == 'L' );                   /* NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay)*/                                                             \
	assert( transA == 'N' || transA == 'T' || transA == 'C' );  /* NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay)*/                                                             \
	assert( diag   == 'U' || diag   == 'N' );                   /* NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay)*/                                                             \
	BOOST_MULTI_ASSERT1( m >= 0 && n >= 0 );                                                                                                                                                                           \
	using std::max;                                                                                                                                                                                           \
	if(side == 'L') {BOOST_MULTI_ASSERT1( lda >= max(ssize_t{1}, m) );}   /* NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay)*/                                                         \
	if(side == 'R') {BOOST_MULTI_ASSERT1( lda >= max(ssize_t{1}, n) );}                                                                                                                                                 \
	BOOST_MULTI_ASSERT1( ldb >= max(ssize_t{1}, m) );                                                                                                                                                                   \
	BLAS(T##trsm)(side, uplo, transA, diag, BC(m), BC(n), alpha, reinterpret_cast<T const*>(static_cast<AA*>(aa)), BC(lda), reinterpret_cast<T*>(static_cast<BB*>(bb)), BC(ldb));   /*NOLINT(cppcoreguidelines-pro-type-reinterpret-cast,bugprone-macro-parentheses)*/                                                                  \
}                                                                                                                                                                                                                 \

xtrsm(s) xtrsm(d) xtrsm(c) xtrsm(z)  // NOLINT(modernize-use-constraints,readability-function-cognitive-complexity) : 29 of 25
#undef xtrsm

xsyrk(s) xsyrk(d) xsyrk(c) xsyrk(z)  // NOLINT(modernize-use-constraints) TODO(correaa) for C++20
#undef xsyrk
	              xherk(c) xherk(z)  // NOLINT(modernize-use-constraints) TODO(correaa) for C++20

} // end namespace core

#undef xherk

#undef BC

struct context { // stateless (and thread safe)
	template<class... As>
	static auto scal(As... args)
	->decltype(core::scal(args...)) {
		return core::scal(args...); }

	template<class... As>
	static auto copy(As... args) noexcept
	->decltype(core::copy(args...)) {
		return core::copy(args...); }

	template<class... As>
	static auto swap(As... args) noexcept
	->decltype(core::swap(args...)) {
		return core::swap(args...); }

	template<class... As>
	static auto axpy(As... args)
	->decltype(core::axpy(args...)) {
		return core::axpy(args...); }

	template<class... As>
	static auto dot(As... args)
	->decltype(core::dot(args...)) {
		return core::dot(args...); }

	template<class... As>
	static auto dotc(As... args)
	->decltype(core::dotc(args...)) {
		return core::dotc(args...); }

	template<class... As>
	static auto dotu(As... args)
	->decltype(core::dotu(args...)) {
		return core::dotu(args...); }

	template<class... As>
	static auto gemm(As&&... args)
	->decltype(core::gemm(std::forward<As>(args)...)) {
		return core::gemm(std::forward<As>(args)...); }

	template<class... As>
	static auto gemv(As&&... args)
	->decltype(core::gemv(std::forward<As>(args)...)) {
	    return core::gemv(std::forward<As>(args)...); }

	template<class... As>
	static auto asum(As... args)
	->decltype(core::asum(args...)) {
		return core::asum(args...); }

	template<class... As>
	static auto nrm2(As... args)
	-> decltype(auto)
	//->decltype(core::nrm2(args...)) {
	{   return core::nrm2(args...); }

	template<class... As>
	static auto trsm(As&&... args)  // TODO(correaa) remove &&
	->decltype(core::trsm(std::forward<As>(args)...)) {
		return core::trsm(std::forward<As>(args)...); }

	template<class... As>
	static auto herk(As&&... args)
	->decltype(core::herk(std::forward<As>(args)...)) {
		return core::herk(std::forward<As>(args)...); }
};

template<class Context> struct is_context    : std::false_type {};

template<> struct is_context<context>        : std::true_type  {};
template<> struct is_context<context&&>      : std::true_type  {};
template<> struct is_context<context&>       : std::true_type  {};
template<> struct is_context<context const&> : std::true_type  {};

template<> struct is_context<void*&> : std::true_type {};

template<class TPtr,
	std::enable_if_t<std::is_convertible_v<TPtr, typename std::pointer_traits<TPtr>::element_type*>, int> =0>  // NOLINT(modernize-use-constraints) TODO(correaa) for C++20
auto default_context_of(TPtr const& /*unused*/) -> blas::context* {
	static blas::context dc;
	return &dc;
}

} // end namespace boost::multi::blas

#endif

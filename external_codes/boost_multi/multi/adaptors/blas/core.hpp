#ifndef MULTI_ADAPTORS_BLAS_CORE_HPP // -*-indent-tabs-mode:t;c-basic-offset:4;tab-width:4;autowrap:nil;-*-
#define MULTI_ADAPTORS_BLAS_CORE_HPP
// © Alfredo A. Correa 2019-2021

// https://software.intel.com/en-us/articles/intel-mkl-link-line-advisor

//#include <cblas/cblas.h> // consider being replaceable by cblas.h

//#include<iostream> // debug
#include<cassert>
#include<complex>
#include<cstdint> // int64_t
#include<cstring> // std::memcpy
#include<limits> // numeric_limits
#include<type_traits> // is_convertible

#include "../../config/MARK.hpp"

#include "../blas/traits.hpp"

#if 0
	#define MULTI_ASSERT1(ExpR)              assert       (ExpR)
	#define MULTI_ASSERT2(ExpR, DescriptioN) MULTI_ASSERT1(ExpR && ##DescriptioN)
#else
	#if not defined(NDEBUG)
		#include<stdexcept>
		#include<string>
		#define MULTI_ASSERT1(ExpR)              (void)((ExpR)?0:throw std::logic_error("\n" __FILE__ ":"+std::to_string(__LINE__)+"::\n"+std::string(__PRETTY_FUNCTION__)+"\nLogic assertion `" #ExpR "' failed.")) /*NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay)*/
		#define MULTI_ASSERT2(ExpR, DescriptioN) (void)((ExpR)?0:throw std::DescriptioN("\n" __FILE__ ":"+std::to_string(__LINE__)+"::\n"+std::string(__PRETTY_FUNCTION__)+"\nLogic assertion `" #ExpR "' failed."))
	#else
		#define MULTI_ASSERT1(ExpR)              assert(ExpR)
		#define MULTI_ASSERT2(ExpR, DescriptioN) assert(EXpR)
	#endif
#endif

#ifdef CBLAS_H
#define BLAS(NamE) cblas_##NamE
#else
#define BLAS(NamE) NamE##_
extern "C"{

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


using Complex_float  = struct { float  real, imag; };
using Complex_double = struct { double real, imag; };
//typedef struct { float  real, imag; } Complex_float ;
//typedef struct { double real, imag; } Complex_double;

#define C Complex_float // _Complex s
#define Z Complex_double // _Complex d

#if defined(MULTI_BLAS_INT)
	#if   MULTI_BLAS_INT==32
		#define INT int32_t
	#elif MULTI_BLAS_INT==64
		#define INT int64_t
	#else
		#define INT int32_t // 32bit safe? pesimistic?
	#endif
#else
	#define INT int32_t // 32bit safe? pesimistic?
#endif

namespace core{
	using size_t = INT;
	using ssize_t = std::make_signed_t<size_t>;
} // end namespace core

#define INTEGER INT const&
#define N INTEGER n
#define INCX INTEGER incx
#define INCY INTEGER incy

static_assert(sizeof(INT)==32/8 or sizeof(INT)==64/8, "please set MULTI_BLAS_INT to int32_t or int64_t");

// TODO(correaa) indent declarations like here https://www.netlib.org/lapack/lug/node145.html

#define xROTG(T1, T2)     v BLAS(   T1##rotg)(                           T1 const*, T1 const*, T2*, T1*)              // NOLINT(bugprone-macro-parentheses) : macro arg expands to type
#define xROTMG(T)         v BLAS(   T##rotmg)(                           T*, T*, T*, T const&, T(&param)[5])          // NOLINT(bugprone-macro-parentheses) : macro arg expands to type
#define xROT(TT, T, S)    v BLAS(  TT##rot  )(N,              T       *x, INCX, T       *y, INCY, S const&, S const&) // NOLINT(bugprone-macro-parentheses) : macro arg expands to type
#define xROTM(T)          v BLAS(   T##rotm )(N, T* x, INCX, T* y, INCY, T const(&p)[5])                              // NOLINT(bugprone-macro-parentheses) : macro arg expands to type
#define xSWAP(T)          v T ##swap##_ (N,              T       *x, INCX, T       *y, INCY)                          // NOLINT(bugprone-macro-parentheses) : macro arg expands to type
#define xSCAL(TT, TA, TX) v TT##scal##_ (N, TA const& a, TX      *x, INCX                  )                          // NOLINT(bugprone-macro-parentheses) : macro arg expands to type
#define xCOPY(T)          v T ##copy##_ (N,              T const *x, INCX, T       *y, INCY)                          // NOLINT(bugprone-macro-parentheses) : macro arg expands to type
#define xAXPY(T)          v T ##axpy##_ (N,  T const* a, T const *x, INCX, T       *y, INCY)                          // NOLINT(bugprone-macro-parentheses) : macro arg expands to type
#define xDOT(R, TT, T)    auto BLAS(  TT##dot  )(N,              T const *x, INCX, T const *y, INCY) -> R
// PGI/NVC++ compiler uses a blas version that needs -DRETURN_BY_STACK
#if defined(RETURN_BY_STACK) || (defined(FORTRAN_COMPLEX_FUNCTIONS_RETURN_VOID) && FORTRAN_COMPLEX_FUNCTIONS_RETURN_VOID)
#define xDOTU(R, T)       v    T##dotu ##_ (R*, N,              T const *x, INCX, T const *y, INCY)
#define xDOTC(R, T)       v    T##dotc ##_ (R*, N,              T const *x, INCX, T const *y, INCY)
#else
#define xDOTU(R, T)       auto    T ##dotu##_ (    N,              T const *x, INCX, T const *y, INCY) -> R
#define xDOTC(R, T)       auto    T ##dotc##_ (    N,              T const *x, INCX, T const *y, INCY) -> R
#endif
#define xxDOT(TT, T)    auto    TT##dot ##_ (    N,  T const& a, T const *x, INCX, T const *y, INCY) -> T
#define xNRM2(R, TT, T) auto    TT##nrm2##_ (    N,              T const *x, INCX                  ) -> R
#define xASUM(R, TT, T) auto    TT##asum##_ (    N,              T const *x, INCX                  ) -> R
#define IxAMAX(T)       auto i##T ##amax##_ (    N,              T const* x, INCX                  ) -> INT

xROTG(s, s)   ; xROTG(d,d)    ;// MKL extension xROTG(c, s); xROTG(z, d);
xROTMG(s)     ; xROTMG(d)     ;
xROT(s, s, s) ; xROT(d, d, d) ;                 xROT(cs, c, s); xROT(zd, z, d);
xROTM(s)      ; xROTM(d)      ;
xSWAP(s)      ; xSWAP(d)      ; xSWAP(c)      ; xSWAP(z);
xSCAL(s, s, s); xSCAL(d, d, d); xSCAL(c, c, c); xSCAL(z, z, z); xSCAL(zd, d, z); xSCAL(cs, s, c);
xCOPY(s)      ; xCOPY(d)      ; xCOPY(c)      ; xCOPY(z)      ;
xAXPY(s)      ; xAXPY(d)      ; xAXPY(c)      ; xAXPY(z)      ;
xDOT(s, s, s); xDOT(d, d, d);                                   xDOT(d, ds, s);

xDOTU(C, c); xDOTU(Z, z);
//xDOTU(c, c); xDOTU(z, z);

xDOTC(C, c); xDOTC(Z, z);
xxDOT(sds, s);
xNRM2(s, s, s); xNRM2(d, d, d); xNRM2(s, sc, c); xNRM2(d, dz, z);
xASUM(s, s, s); xASUM(d, d, d); xASUM(s, sc, c); xASUM(d, dz, z);
IxAMAX(s); IxAMAX(d); IxAMAX(c); IxAMAX(z);

#define TRANS const char& trans
#define NR INTEGER nr
#define NC INTEGER nc
#define LDA INTEGER lda
#define UPLO const char& uplo
#define DIAG const char& diag

#define xGEMV(T) void  T## gemv ##_ (      TRANS,       NR, NC, T const& a, T const* A, LDA, T const* X, INCX, T const& beta, T*       Y, INCY           ) // NOLINT(bugprone-macro-parentheses) : macro arg expands to type
#define xGER(T)  void  T## ger  ##_ (                   NR, NC, T const& a,                  T const* X, INCX,                T const* Y, INCY, T* A, LDA) // NOLINT(bugprone-macro-parentheses) : macro arg expands to type
#define xGERU(T) void  T## geru ##_ (                   NR, NC, T const& a,                  T const* X, INCX,                T const* Y, INCY, T* A, LDA) // NOLINT(bugprone-macro-parentheses) : macro arg expands to type
#define xGERC(T) void  T## gerc ##_ (                   NR, NC, T const& a,                  T const* X, INCX,                T const* Y, INCY, T* A, LDA) // NOLINT(bugprone-macro-parentheses) : macro arg expands to type
#define xTRSV(T) void  T## trsv ##_ (UPLO, TRANS, DIAG, N,                  T const* A, LDA, T* X      , INCX                                            ) // NOLINT(bugprone-macro-parentheses) : macro arg expands to type

xGEMV(s); xGEMV(d); xGEMV(c); xGEMV(z);
xGER(s); xGER(d);
xGERU(c); xGERU(z);
xGERC(c); xGERC(z);
xTRSV(s); xTRSV(d); xTRSV(c); xTRSV(z);

#define TRANSA const char& transa
#define TRANSB const char& transb
#define NK  INTEGER nk
#define LDB INTEGER ldb
#define LDC INTEGER ldc

#define SIDE const char& side

#define xGEMM(T)     void T ##gemm ##_ (            TRANSA, TRANSB,       NR, NC, NK, T  const& a, T const* A, LDA, T const* B, LDB, T  const& b     , T const* CC, LDC)
#define xSYRK(T)     void T ##syrk ##_ (      UPLO, TRANSA,               NR, NK,     T  const& a, T const* A, LDA,                  T  const& b     , T*       CC, LDC) // NOLINT(bugprone-macro-parentheses) : macro arg expands to type
#define xHERK(TT, T) void T ##herk ##_ (      UPLO, TRANSA,               NR, NK,     TT const& a, T const* A, LDA,                  TT const& b     , T*       CC, LDC) // NOLINT(bugprone-macro-parentheses) : macro arg expands to type
#define xTRSM(T)     void T ##trsm ##_ (SIDE, UPLO, TRANSA,         DIAG, NR, NK,     T  const& a, T const* A, LDA,                  T  const* B, LDB                  )

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
#undef xxDOT
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
}
#endif

namespace boost{
namespace multi{
namespace blas{

//template<class T> struct complex_ptr{
//	std::complex<T>* impl_;
//	template<class TT, class=std::enable_if_t<sizeof(*TT{})==sizeof(std::complex<T>) and sizeof(*TT{})==sizeof(TT{}->real())+sizeof(TT{}->imag())>>
//	explicit complex_ptr(TT tt) : impl_{reinterpret_cast<std::complex<T>*>(tt)}{}
//	complex_ptr(complex_ptr const&) = delete;
//	operator std::complex<T>*() const{return   impl_;}
//	std::complex<T>& operator*() const{return *impl_;}
//};

//template<class T> struct complex_const_ptr{
//	std::complex<T> const* impl_;
//	template<class TT, class=std::enable_if_t<sizeof(*TT{})==sizeof(std::complex<T>) and sizeof(*TT{})==sizeof(TT{}->real())+sizeof(TT{}->imag())>>
//	explicit complex_const_ptr(TT tt) : impl_{reinterpret_cast<std::complex<T> const*>(tt)}{} // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast)
//	complex_const_ptr(complex_const_ptr const&) = delete;
//	explicit operator std::complex<T> const*() const{return impl_;}
//	auto operator*() const -> std::complex<T> const&{return *impl_;}
//};

//template<class T> struct add_ptr{using type = T*;};
//template<class T> struct add_const_ptr{using type = T const*;};

//template<class T> struct add_ptr<std::complex<T>>{using type = complex_ptr<T>;};
//template<class T> struct add_const_ptr<std::complex<T>>{using type = complex_const_ptr<T>;};

//template<class T> using add_ptr_t = typename add_ptr<T>::type;
//template<class T> using add_const_ptr_t = typename add_const_ptr<T>::type;

//namespace t{
	using s = float;
	using d = double;
	using c = std::complex<s>; //using C = Complex_float ;
	using z = std::complex<d>; //using Z = Complex_double;
	using v = void;
//} // end namespace types

//using namespace types;

#define BC(x) [](auto xx){assert(xx>=std::numeric_limits<INT>::min() and xx<std::numeric_limits<INT>::max()); return xx;}(x) /*NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay)*/

//#define xrotg(T1, T2)                       v   rotg (T1 const& a, T1 const& b, T2& cc, T1& ss                                       ){     BLAS(T1##rotg )(const_cast<T1*>(&a), const_cast<T1*>(&b), &cc, &ss);  }
//#define xrotmg(T)                           v   rotmg(T& d1, T& d2, T& A, T const& B, T(&p)[5]                                       ){     BLAS( T##rotmg)(&d1, &d2, &A, B, p);                                  }
//#define xrot(T, TT, CS)   template<class S> v   rot  (S n,       T       *x, S incx, T       *y, S incy, CS const& cos, CS const& sin){     BLAS(TT##rot  )(BC(n),    x, BC(incx), y, BC(incy), cos, sin);    }
//#define xrotm(T)          template<class S> v   rotm (S n,       T       *x, S incx, T       *y, S incy, T const(&p)[5]              ){     BLAS( T##rotm )(BC(n),    x, BC(incx), y, BC(incy), p);               }
//#define xswap(T)          template<class S> v   swap (S n,       T       *x, S incx, T       *y, S incy                              ){     BLAS( T##swap )(BC(n),    x, BC(incx), y, BC(incy));                  }
//#define xscal(XX, TA, TX) template<class S> TX* scal (INT n, TA const* a, TX       *x, S incx                                        ){     BLAS(XX##scal )(BC(n), *a, x, BC(incx)             ); return x+n*incx;}
//#define xcopy(T)          v   copy (INT n,              T  const *x, INT incx, T       *y, INT incy                              ){     BLAS( T##copy )(BC(n),    x, BC(incx), y, BC(incy));                  }
//#define xaxpy(T)          template<class S> T*  axpy (S n, T  a, T const *x, S incx, T       *y, S incy                          ){     BLAS( T##axpy )(BC(n), a, x, BC(incx), y, BC(incy)); return y+n*incy; }
#if 0
#define xdot(R, TT, T)    template<class S> v   dot  (S n,       T const* x, S incx, T const* y, S incy, R* r                    ){\
		MULTI_MARK_SCOPE("cpu_dot"); *r = BLAS(TT##dot  )(BC(n),    x, BC(incx), y, BC(incy));                  }
#endif

//xrotg(s, s)    xrotg(d, d) //MKL extension xrotg(c, s); xrotg(z, d);
//xrotmg(s)      xrotmg(d)
//xrot(s, s, s)  xrot(d, d, d)  xrot(c, cs, s) xrot(z, zd, d)
//xrotm(s)       xrotm(d)
//xswap(s)       xswap(d)       xswap(c)       xswap(z)

namespace core{

//xscal(s, s, s) xscal(d, d, d) xscal(c, c, c) xscal(z, z, z) xscal(zd, d, z) xscal(cs, s, c)

using std::enable_if_t;
using std::is_assignable;

template<class SXP, class SYP, class SX = typename std::pointer_traits<SXP>::element_type, class SY = typename std::pointer_traits<SYP>::element_type, enable_if_t<is_s<SX>{} and is_s<SY>{} and is_assignable<SY&, SX&>{},int> =0> void swap(ssize_t n, SX* x, ptrdiff_t incx, SY* y, ptrdiff_t incy){BLAS(sswap)(n, (             float  *)(x), incx, (             float  *)(y), incy);} // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast,cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays,google-readability-casting)
template<class DXP, class DYP, class DX = typename std::pointer_traits<DXP>::element_type, class DY = typename std::pointer_traits<DYP>::element_type, enable_if_t<is_d<DX>{} and is_d<DY>{} and is_assignable<DY&, DX&>{},int> =0> void swap(ssize_t n, DX* x, ptrdiff_t incx, DY* y, ptrdiff_t incy){BLAS(dswap)(n, (             double *)(x), incx, (             double *)(y), incy);} // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast,cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays,google-readability-casting)
template<class CXP, class CYP, class CX = typename std::pointer_traits<CXP>::element_type, class CY = typename std::pointer_traits<CYP>::element_type, enable_if_t<is_c<CX>{} and is_c<CY>{} and is_assignable<CY&, CX&>{},int> =0> void swap(ssize_t n, CX* x, ptrdiff_t incx, CY* y, ptrdiff_t incy){BLAS(cswap)(n, (std::complex<float >*)(x), incx, (std::complex<float >*)(y), incy);} // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast,cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays,google-readability-casting)
template<class ZXP, class ZYP, class ZX = typename std::pointer_traits<ZXP>::element_type, class ZY = typename std::pointer_traits<ZYP>::element_type, enable_if_t<is_z<ZX>{} and is_z<ZY>{} and is_assignable<ZY&, ZX&>{},int> =0> void swap(ssize_t n, ZX* x, ptrdiff_t incx, ZY* y, ptrdiff_t incy){BLAS(zswap)(n, (std::complex<double>*)(x), incx, (std::complex<double>*)(y), incy);} // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast,cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays,google-readability-casting)

template<class SX, class SY, enable_if_t<is_s<SX>{} and is_s<SY>{} and is_assignable<SY&, SX&>{},int> =0> void copy(ssize_t n, SX* x, ptrdiff_t incx, SY* y, ptrdiff_t incy){BLAS(scopy)(n, (             float   const*)(x), incx, (             float  *)(y), incy);} // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast,cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays,google-readability-casting)
template<class DX, class DY, enable_if_t<is_d<DX>{} and is_d<DY>{} and is_assignable<DY&, DX&>{},int> =0> void copy(ssize_t n, DX* x, ptrdiff_t incx, DY* y, ptrdiff_t incy){BLAS(dcopy)(n, (             double  const*)(x), incx, (             double *)(y), incy);} // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast,cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays,google-readability-casting)
template<class CX, class CY, enable_if_t<is_c<CX>{} and is_c<CY>{} and is_assignable<CY&, CX&>{},int> =0> void copy(ssize_t n, CX* x, ptrdiff_t incx, CY* y, ptrdiff_t incy){BLAS(ccopy)(n, (std::complex<float > const*)(x), incx, (std::complex<float >*)(y), incy);} // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast,cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays,google-readability-casting)
template<class ZX, class ZY, enable_if_t<is_z<ZX>{} and is_z<ZY>{} and is_assignable<ZY&, ZX&>{},int> =0> void copy(ssize_t n, ZX* x, ptrdiff_t incx, ZY* y, ptrdiff_t incy){BLAS(zcopy)(n, (std::complex<double> const*)(x), incx, (std::complex<double>*)(y), incy);} // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast,cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays,google-readability-casting)

// TODO(correaa) : add mixed-type scal (zdscal, csscal)
template<class ALPHAP, class SXP, class SX = typename std::pointer_traits<SXP>::element_type, class ALPHA = typename std::pointer_traits<ALPHAP>::element_type, enable_if_t<is_s<SX>{} and is_s<ALPHA>{} and is_assignable<SX&, decltype(*ALPHAP{}*SX{})>{}>* = nullptr> void scal(ssize_t n, ALPHAP a, SXP xp, ptrdiff_t incx){BLAS(sscal)(n, *(             float   const*)a, (             float  *)xp, incx);} // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast,cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays,google-readability-casting)
template<class ALPHAP, class DXP, class DX = typename std::pointer_traits<DXP>::element_type, class ALPHA = typename std::pointer_traits<ALPHAP>::element_type, enable_if_t<is_d<DX>{} and is_d<ALPHA>{} and is_assignable<DX&, decltype(*ALPHAP{}*DX{})>{}>* = nullptr> void scal(ssize_t n, ALPHAP a, DXP xp, ptrdiff_t incx){BLAS(dscal)(n, *(             double  const*)a, (             double *)xp, incx);} // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast,cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays,google-readability-casting)
template<class ALPHAP, class CXP, class CX = typename std::pointer_traits<CXP>::element_type, class ALPHA = typename std::pointer_traits<ALPHAP>::element_type, enable_if_t<is_c<CX>{} and is_c<ALPHA>{} and is_assignable<CX&, decltype(*ALPHAP{}*CX{})>{}>* = nullptr> void scal(ssize_t n, ALPHAP a, CXP xp, ptrdiff_t incx){BLAS(cscal)(n, *(std::complex<float > const*)a, (std::complex<float >*)xp, incx);} // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast,cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays,google-readability-casting)
template<class ALPHAP, class ZXP, class ZX = typename std::pointer_traits<ZXP>::element_type, class ALPHA = typename std::pointer_traits<ALPHAP>::element_type, enable_if_t<is_z<ZX>{} and is_z<ALPHA>{} and is_assignable<ZX&, decltype(*ALPHAP{}*ZX{})>{}>* = nullptr> void scal(ssize_t n, ALPHAP a, ZXP xp, ptrdiff_t incx){BLAS(zscal)(n, *(std::complex<double> const*)a, (std::complex<double>*)xp, incx);} // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast,cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays,google-readability-casting)

//xdot(s, s, s)  xdot(d, d, d)                                xdot(d, ds, s)

//template<class SXP, class SYP, class SRP, enable_if_t<is_s<SX>{} and is_s<SY>{} and is_assignable<SY&, SX&>{}>* = nullptr> void dot(size_t n, SXP xp, ptrdiff_t incx, SYP yp, ptrdiff_t incy, RP rp

using std::pointer_traits;
using std::enable_if_t;
using std::is_convertible_v;

#define xaxpy(T) \
template<class ALPHA, class SXP, class SX = typename pointer_traits<SXP>::element_type, class SYP, class SY = typename pointer_traits<SYP>::element_type, enable_if_t< \
	is_##T<ALPHA>{} and is_##T<SX>{} and is_##T<SY>{} and is_assignable<SY&, decltype(ALPHA{}*SX{})>{} \
	and is_convertible_v<SXP, SX*> and is_convertible_v<SYP, SY*> \
, int> =0> \
void axpy(size_t n, ALPHA const* a, SXP x, size_t incx, SYP y, size_t incy){BLAS(T##axpy)(n, (T const *)a, (T const*)static_cast<SX*>(x), incx, (T*)static_cast<SY*>(y), incy);}

xaxpy(s)       xaxpy(d)       xaxpy(c)       xaxpy(z)
#undef  xaxpy
//template<class A, class SX, class SY, enable_if_t<is_s<SX>{} and is_s<SY>{} and is_assignable<SY&, decltype(A{}*SX{})>{}, int> =0> void axpy(size_t n, A a, SX* x, size_t incx, SY* y, size_t incy){BLAS(saxpy)(n, a, (s const*)(x), incx, (s*)(y), incy);}
//template<class A, class DX, class DY, enable_if_t<is_d<DX>{} and is_d<DY>{} and is_assignable<DY&, decltype(A{}*DX{})>{}, int> =0> void axpy(size_t n, A a, DX* x, size_t incx, DY* y, size_t incy){BLAS(daxpy)(n, a, (d const*)(x), incx, (d*)(y), incy);}
//template<class A, class CX, class CY, enable_if_t<is_c<CX>{} and is_c<CY>{} and is_assignable<CY&, decltype(A{}*CX{})>{}, int> =0> void axpy(size_t n, A a, CX* x, size_t incx, CY* y, size_t incy){BLAS(caxpy)(n, a, (c const*)(x), incx, (c*)(y), incy);}
//template<class A, class ZX, class ZY, enable_if_t<is_z<ZX>{} and is_z<ZY>{} and is_assignable<ZY&, decltype(A{}*ZX{})>{}, int> =0> void axpy(size_t n, A a, ZX* x, size_t incx, ZY* y, size_t incy){BLAS(zaxpy)(n, a, (z const*)(x), incx, (z*)(y), incy);}

} // end namespace core

//template<class R, class S, class T> R dot(S n, T const* x, S incx, T const* y, S incy){
//	R ret;
//	dot(n, x, incx, y, incy, &ret);
//	return ret;
//}

//template<class S, class T> T dot(S n, T const* x, S incx, T const* y, S incy){
//	return dot<T, S, T>(n, x, incx, y, incy);
//}

#undef xrotg
#undef xrot
#undef xswap
#undef xscal
#undef xcopy
#undef xaxpy
#undef xdot

#ifndef CBLAS_H

namespace core{

using std::enable_if_t;
using std::is_assignable;

// PGI/NVC++ compiler uses a blas version that needs -DRETURN_BY_STACK
#if defined(RETURN_BY_STACK) || (defined(FORTRAN_COMPLEX_FUNCTIONS_RETURN_VOID) && FORTRAN_COMPLEX_FUNCTIONS_RETURN_VOID)
template<class X, class Y, class R, enable_if_t<is_s<X>{} and is_s<Y>{} and is_assignable<R&, decltype(0.+X{}*Y{}+X{}*Y{})>{}, int> =0> void dotu(size_t n, X* x, size_t incx, Y* y, size_t incy, R* r){BLAS(sdot )((float *)r, n, (s const*)x, incx, (s const*)y, incy);}
template<class X, class Y, class R, enable_if_t<is_d<X>{} and is_d<Y>{} and is_assignable<R&, decltype(0.+X{}*Y{}+X{}*Y{})>{}, int> =0> void dotu(size_t n, X* x, size_t incx, Y* y, size_t incy, R* r){BLAS(ddot )((double*)r, n, (d const*)x, incx, (d const*)y, incy);}

template<class X, class Y, class R, enable_if_t<is_c<X>{} and is_c<Y>{} and is_assignable<R&, decltype(0.+X{}*Y{}+X{}*Y{})>{}, int> =0> void dotu(size_t n, X* x, size_t incx, Y* y, size_t incy, R* r){BLAS(cdotu)((Complex_float *)r, n, (c const*)x, incx, (c const*)y, incy);}
template<class X, class Y, class R, enable_if_t<is_z<X>{} and is_z<Y>{} and is_assignable<R&, decltype(0.+X{}*Y{}+X{}*Y{})>{}, int> =0> void dotu(size_t n, X* x, size_t incx, Y* y, size_t incy, R* r){BLAS(zdotu)((Complex_double*)r, n, (z const*)x, incx, (z const*)y, incy);}

template<class X, class Y, class R, enable_if_t<is_c<X>{} and is_c<Y>{} and is_assignable<R&, decltype(0.+X{}*Y{}+X{}*Y{})>{}, int> =0> void dotc(size_t n, X* x, size_t incx, Y* y, size_t incy, R* r){BLAS(cdotc)(reinterpret_cast<Complex_float *>(r), n, (c const*)x, incx, (c const*)y, incy);}
template<class X, class Y, class R, enable_if_t<is_z<X>{} and is_z<Y>{} and is_assignable<R&, decltype(0.+X{}*Y{}+X{}*Y{})>{}, int> =0> void dotc(size_t n, X* x, size_t incx, Y* y, size_t incy, R* r){BLAS(zdotc)(reinterpret_cast<Complex_double*>(r), n, (z const*)x, incx, (z const*)y, incy);} // NOLINT(ppcoreguidelines-pro-type-reinterpret-cast) : adapt types
#else
template<class XP, class X = typename std::pointer_traits<XP>::element_type, class YP, class Y = typename std::pointer_traits<YP>::element_type, class RP, class R = typename std::pointer_traits<RP>::element_type, enable_if_t<is_s<X>{} and is_s<Y>{} and is_assignable<R&, decltype(0.+X{}*Y{}+X{}*Y{})>{}, int> =0> void dot (ssize_t n, XP x, ptrdiff_t incx, YP y, ptrdiff_t incy, RP r){auto rr = BLAS(sdot )(n, (s const*)static_cast<X*>(x), incx, (s const*)static_cast<Y*>(y), incy); std::memcpy(reinterpret_cast<float *>(static_cast<R*>(r)), &rr, sizeof(rr)); static_assert(sizeof(rr)==sizeof(*r));} // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast,cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays,google-readability-casting)
template<class XP, class X = typename std::pointer_traits<XP>::element_type, class YP, class Y = typename std::pointer_traits<YP>::element_type, class RP, class R = typename std::pointer_traits<RP>::element_type, enable_if_t<is_d<X>{} and is_d<Y>{} and is_assignable<R&, decltype(0.+X{}*Y{}+X{}*Y{})>{}, int> =0> void dot (ssize_t n, XP x, ptrdiff_t incx, YP y, ptrdiff_t incy, RP r){auto rr = BLAS(ddot )(n, (d const*)static_cast<X*>(x), incx, (d const*)static_cast<Y*>(y), incy); std::memcpy(reinterpret_cast<double*>(static_cast<R*>(r)), &rr, sizeof(rr)); static_assert(sizeof(rr)==sizeof(*r));} // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast,cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays,google-readability-casting)

template<class XP, class X = typename std::pointer_traits<XP>::element_type, class YP, class Y = typename std::pointer_traits<YP>::element_type, class RP, class R = typename std::pointer_traits<RP>::element_type, enable_if_t<is_c<X>{} and is_c<Y>{} and is_assignable<R&, decltype(0.+X{}*Y{}+X{}*Y{})>{}, int> =0> void dotu(ssize_t n, XP x, ptrdiff_t incx, YP y, ptrdiff_t incy, RP r){auto rr = BLAS(cdotu)(n, (c const*)static_cast<X*>(x), incx, (c const*)static_cast<Y*>(y), incy); std::memcpy(reinterpret_cast<float (*)[2]>(static_cast<R*>(r)), &rr, sizeof(rr)); static_assert(sizeof(rr)==sizeof(*r));} // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast,cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays,google-readability-casting)
template<class XP, class X = typename std::pointer_traits<XP>::element_type, class YP, class Y = typename std::pointer_traits<YP>::element_type, class RP, class R = typename std::pointer_traits<RP>::element_type, enable_if_t<is_z<X>{} and is_z<Y>{} and is_assignable<R&, decltype(0.+X{}*Y{}+X{}*Y{})>{}, int> =0> void dotu(ssize_t n, XP x, ptrdiff_t incx, YP y, ptrdiff_t incy, RP r){auto rr = BLAS(zdotu)(n, (z const*)static_cast<X*>(x), incx, (z const*)static_cast<Y*>(y), incy); std::memcpy(reinterpret_cast<double(*)[2]>(static_cast<R*>(r)), &rr, sizeof(rr)); static_assert(sizeof(rr)==sizeof(*r));} // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast,cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays,google-readability-casting)

template<class XP, class X = typename std::pointer_traits<XP>::element_type, class YP, class Y = typename std::pointer_traits<YP>::element_type, class RP, class R = typename std::pointer_traits<RP>::element_type, enable_if_t<is_c<X>{} and is_c<Y>{} and is_assignable<R&, decltype(0.+X{}*Y{}+X{}*Y{})>{}, int> =0> void dotc(ssize_t n, XP x, ptrdiff_t incx, YP y, ptrdiff_t incy, RP r){auto rr = BLAS(cdotc)(n, (c const*)static_cast<X*>(x), incx, (c const*)static_cast<Y*>(y), incy); std::memcpy(reinterpret_cast<float (*)[2]>(static_cast<R*>(r)), &rr, sizeof(rr)); static_assert(sizeof(rr)==sizeof(*r));} // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast,cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays,google-readability-casting)
template<class XP, class X = typename std::pointer_traits<XP>::element_type, class YP, class Y = typename std::pointer_traits<YP>::element_type, class RP, class R = typename std::pointer_traits<RP>::element_type, enable_if_t<is_z<X>{} and is_z<Y>{} and is_assignable<R&, decltype(0.+X{}*Y{}+X{}*Y{})>{}, int> =0> void dotc(ssize_t n, XP x, ptrdiff_t incx, YP y, ptrdiff_t incy, RP r){auto rr = BLAS(zdotc)(n, (z const*)static_cast<X*>(x), incx, (z const*)static_cast<Y*>(y), incy); std::memcpy(reinterpret_cast<double(*)[2]>(static_cast<R*>(r)), &rr, sizeof(rr)); static_assert(sizeof(rr)==sizeof(*r));} // NOLINT(cppcoreguidelines-pro-type-reinterpret-cast,cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays,google-readability-casting)
#endif

} // end namespace core
#else
// TODO: make cblas version
#define xdotu(T) template<class S> v dotu(S n, add_const_ptr_t<T> x, S incx, add_const_ptr_t<T> y, S incy, add_ptr_t<T> r){BLAS(T##dotu_sub)(BC(n), x, BC(incx), y, BC(incy), r);}
#define xdotc(T) template<class S> v dotc(S n, add_const_ptr_t<T> x, S incx, add_const_ptr_t<T> y, S incy, add_ptr_t<T> r){BLAS(T##dotc_sub)(BC(n), x, BC(incx), y, BC(incy), r);}

namespace core{
	xdotu(c) xdotu(z)
	xdotc(c) xdotc(z)
}

#undef xdotu
#undef xdotc
#endif

namespace core{
	template<class S> auto dot(S n, s const& b, s const* x, S incx, s const* y, S incy) -> s{return BLAS(sdsdot)(BC(n), b, x, BC(incx), y, BC(incy));}

//template<class S> void dot(S n, s const& b, s const* x, S incx, s const* y, S incy, s* result){*result = BLAS(sdsdot)(BC(n), b, x, BC(incx), y, BC(incy));}
} // end namespace core

//#define xnrm2(R, T, TT) template<class S>    v nrm2 (S n, add_const_ptr_t<T> x, S incx, R* r){*r = BLAS(TT##nrm2  )(BC(n), x, BC(incx));}

#define xasum(T, TT)    template<class S> auto asum (S n, T const* x, S incx){return BLAS(TT##asum  )(BC(n), x, BC(incx));}
#define ixamax(T)       template<class S> auto iamax(S n, T const* x, S incx){return BLAS(i##T##amax)(BC(n), x, BC(incx)) - 1;}
xasum(s, s)    xasum(d, d)                        xasum (c, sc)                  xasum(z, dz)

namespace core{
//	xnrm2(s, s, s) xnrm2(d, d, d)  xnrm2(s, c, sc) xnrm2(d, z, dz)

template<class XP, class X = typename std::pointer_traits<XP>::element_type, class RP, class R = typename std::pointer_traits<RP>::element_type, enable_if_t<is_s<X>{} and is_s<R>{} and std::is_assignable<R&, decltype(X{})>{}           , int> =0> void nrm2(ssize_t n, XP x, ptrdiff_t incx, RP r){auto rr = BLAS(snrm2) (n, (s const*)static_cast<X*>(x), incx); std::memcpy((s*)static_cast<R*>(r), &rr, sizeof(s));} // NOLINT(google-readability-casting)
template<class XP, class X = typename std::pointer_traits<XP>::element_type, class RP, class R = typename std::pointer_traits<RP>::element_type, enable_if_t<is_d<X>{} and is_d<R>{} and std::is_assignable<R&, decltype(X{})>{}           , int> =0> void nrm2(ssize_t n, XP x, ptrdiff_t incx, RP r){auto rr = BLAS(dnrm2) (n, (d const*)static_cast<X*>(x), incx); std::memcpy((s*)static_cast<R*>(r), &rr, sizeof(d));} // NOLINT(google-readability-casting)

template<class XP, class X = typename std::pointer_traits<XP>::element_type, class RP, class R = typename std::pointer_traits<RP>::element_type, enable_if_t<is_c<X>{} and is_s<R>{} and std::is_assignable<R&, decltype(std::norm(X{}))>{}, int> =0> void nrm2(ssize_t n, XP x, ptrdiff_t incx, RP r){auto rr = BLAS(scnrm2)(n, (c const*)static_cast<X*>(x), incx); std::memcpy((s*)static_cast<R*>(r), &rr, sizeof(s));} // NOLINT(google-readability-casting)
template<class XP, class X = typename std::pointer_traits<XP>::element_type, class RP, class R = typename std::pointer_traits<RP>::element_type, enable_if_t<is_z<X>{} and is_d<R>{} and std::is_assignable<R&, decltype(std::norm(X{}))>{}, int> =0> void nrm2(ssize_t n, XP x, ptrdiff_t incx, RP r){auto rr = BLAS(dznrm2)(n, (z const*)static_cast<X*>(x), incx); std::memcpy((s*)static_cast<R*>(r), &rr, sizeof(d));} // NOLINT(google-readability-casting)

//	template<class S>    v nrm2 (S n, typename add_const_ptr<std::complex<double>>::type x, S incx, d* r){*r = BLAS(dznrm2  )(BC(n), x, BC(incx));}
	ixamax(s)      ixamax(d)       ixamax(c)       ixamax(z)
} // end namespace core

#undef xnrm2
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

namespace core{

//xgemv(s) xgemv(d) xgemv(c) xgemv(z)
//xger(s)   xger(d)
//                  xgeru(c) xgeru(z)
//                  xgerc(c) xgerc(z)


using std::enable_if_t;
using std::is_assignable;

template<class A, class M, class X, class B, class Y, enable_if_t<is_s<M>{} and is_s<X>{} and is_s<Y>{} and is_assignable<Y&, decltype(A{}*M{}*X{}+B{}*Y{})>{}, int> =0> void gemv(char trans, size_t m, size_t n, A const& a, M* ma, size_t lda, X* x, size_t incx, B b, Y* y, size_t incy){BLAS(sgemv)(trans, m, n, a, (s const*)ma, lda, (s const*)x, incx, b, (s*)y, incy);} // NOLINT(google-readability-casting)
template<class A, class M, class X, class B, class Y, enable_if_t<is_d<M>{} and is_d<X>{} and is_d<Y>{} and is_assignable<Y&, decltype(A{}*M{}*X{}+B{}*Y{})>{}, int> =0> void gemv(char trans, size_t m, size_t n, A const& a, M* ma, size_t lda, X* x, size_t incx, B b, Y* y, size_t incy){BLAS(dgemv)(trans, m, n, a, (d const*)ma, lda, (d const*)x, incx, b, (d*)y, incy);} // NOLINT(google-readability-casting)
template<class A, class M, class X, class B, class Y, enable_if_t<is_c<M>{} and is_c<X>{} and is_c<Y>{} and is_assignable<Y&, decltype(A{}*M{}*X{}+B{}*Y{})>{}, int> =0> void gemv(char trans, size_t m, size_t n, A const& a, M* ma, size_t lda, X* x, size_t incx, B b, Y* y, size_t incy){BLAS(cgemv)(trans, m, n, a, (c const*)ma, lda, (c const*)x, incx, b, (c*)y, incy);} // NOLINT(google-readability-casting)
template<class A, class M, class X, class B, class Y, enable_if_t<is_z<M>{} and is_z<X>{} and is_z<Y>{} and is_assignable<Y&, decltype(A{}*M{}*X{}+B{}*Y{})>{}, int> =0> void gemv(char trans, size_t m, size_t n, A const& a, M* ma, size_t lda, X* x, size_t incx, B b, Y* y, size_t incy){BLAS(zgemv)(trans, m, n, a, (z const*)ma, lda, (z const*)x, incx, b, (z*)y, incy);} // NOLINT(google-readability-casting)

//template<class SX, class SY, enable_if_t<is_s<SX>{} and is_s<SY>{} and is_assignable<SY&, SX&>{},int> =0> void copy(size_t n, SX* x, size_t incx, SY* y, size_t incy){BLAS(scopy)(n, (             float   const*)(x), incx, (             float  *)(y), incy);}
//template<class DX, class DY, enable_if_t<is_d<DX>{} and is_d<DY>{} and is_assignable<DY&, DX&>{},int> =0> void copy(size_t n, DX* x, size_t incx, DY* y, size_t incy){BLAS(dcopy)(n, (             double  const*)(x), incx, (             double *)(y), incy);}
//template<class CX, class CY, enable_if_t<is_c<CX>{} and is_c<CY>{} and is_assignable<CY&, CX&>{},int> =0> void copy(size_t n, CX* x, size_t incx, CY* y, size_t incy){BLAS(ccopy)(n, (std::complex<float > const*)(x), incx, (std::complex<float >*)(y), incy);}
//template<class ZX, class ZY, enable_if_t<is_z<ZX>{} and is_z<ZY>{} and is_assignable<ZY&, ZX&>{},int> =0> void copy(size_t n, ZX* x, size_t incx, ZY* y, size_t incy){BLAS(zcopy)(n, (std::complex<double> const*)(x), incx, (std::complex<double>*)(y), incy);}


} // end namespace core

template<class T> 
struct blas2{
//	template<class S>
//	static v trsv(char ulA, char transA, char di, S m, T const* A, S lda, T* X, S incx) = delete;
};

template<> struct blas2<s>{template<class... As> static v    trsv(As... as)                              {BLAS(strsv)(as...);}};
template<> struct blas2<d>{template<class... As> static v    trsv(As... as)                              {BLAS(dtrsv)(as...);}};
template<> struct blas2<c>{template<class... As> static v    trsv(As... as)                              {BLAS(ctrsv)(as...);}};
template<> struct blas2<z>{template<class... As> static auto trsv(As... as)->decltype(BLAS(ztrsv)(as...)){BLAS(ztrsv)(as...);}};

namespace core{
	template<typename TconstP, typename TP, typename S=std::size_t, typename C=char> v trsv(C ulA, C transA, C diA, S n, TconstP A, S lda, TP X, S incx){blas2<std::decay_t<typename std::pointer_traits<TP>::element_type>>::trsv(ulA, transA, diA, n, A, lda, X, incx);}
} // end namespace core

//#undef xgemv
#undef xger
#undef xgeru
#undef xgerc

///////////////////////////////////////////////////////////////////////////////
// LEVEL 3

#if 0
#define xsyrk(T) \
template<class UL, class C, class S>             v syrk(        UL ul, C transA,             S n, S k, T    alpha, T const* A, S lda,             T    beta, T* CC, S ldc){ \
	MULTI_MARK_SCOPE("cpu_syrk"); BLAS(T##syrk)(      ul, transA,            BC(n), BC(k), alpha, A, BC(lda),        beta, CC, BC(ldc));}
#endif

namespace core{

using std::is_convertible_v;
using std::pointer_traits;
using std::enable_if_t;
using std::max;

#define xsyrk(T) \
template<class UL, class C, class S, class ALPHA, class AAP, class AA = typename pointer_traits<AAP>::element_type, class BETA, class CCP, class CC = typename pointer_traits<CCP>::element_type, \
enable_if_t< \
	is_##T<AA>{} and is_##T<CC>{} and is_assignable<CC&, decltype(ALPHA{}*AA{}*AA{})>{} and \
	is_convertible_v<AAP, AA*> and is_convertible_v<CCP, CC*> \
, int> =0> \
v syrk(        UL ul, C transA,             S n, S k, ALPHA const* alpha, AAP aa, S lda,             BETA const* beta, CCP cc, S ldc) \
/*=delete;*/ \
{ \
	if(transA == 'N' or transA == 'n') MULTI_ASSERT1( lda >= max(1L, n) ); else MULTI_ASSERT1( lda >= max(1L, k) ); \
	MULTI_ASSERT1( ldc >= max(1L, n) ); \
	MULTI_MARK_SCOPE("cpu_herk"); \
	BLAS(T##syrk)(      ul, transA,            BC(n), BC(k), *(T const*)alpha, aa, BC(lda),        *(T const*)beta, cc, BC(ldc)); \
}

#define xherk(T) \
template<class UL, class C, class S, class ALPHA, class AAP, class AA = typename pointer_traits<AAP>::element_type, class BETA, class CCP, class CC = typename pointer_traits<CCP>::element_type, class Real = typename T::value_type,\
enable_if_t< \
	is_##T<AA>{} and is_##T<CC>{} and is_assignable<CC&, decltype(ALPHA{}*AA{}*AA{})>{} and \
	is_convertible_v<AAP, AA*> and is_convertible_v<CCP, CC*> \
, int> =0> \
v herk(        UL ul, C transA,             S n, S k, ALPHA const* alpha, AAP aa, S lda,             BETA const* beta, CCP cc, S ldc) \
/*=delete;*/ \
{ \
	if(transA == 'N' or transA == 'n'){MULTI_ASSERT1( lda >= max(1L, n) );}else{MULTI_ASSERT1( lda >= max(1L, k) );} /* NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay)*/ \
	MULTI_ASSERT1( ldc >= max(1L, n) ); /* NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay)*/ \
	MULTI_MARK_SCOPE("cpu_herk"); \
	BLAS(T##herk)(      ul, transA,            BC(n), BC(k), *(Real const*)alpha, aa, BC(lda),        *(Real const*)beta, cc, BC(ldc)); \
}

#define xgemm(T) \
template<class ALPHA, class AAP, class AA = typename pointer_traits<AAP>::element_type, class BBP, class BB = typename pointer_traits<BBP>::element_type, class BETA, class CCP, class CC = typename pointer_traits<CCP>::element_type, \
enable_if_t< \
	is_##T<AA>{} and is_##T<BB>{} and is_##T<CC>{} and is_assignable<CC&, decltype(ALPHA{}*AA{}*BB{})>{} and \
	is_convertible_v<AAP, AA*> and is_convertible_v<BBP, BB*> and is_convertible_v<CCP, CC*> \
, int> =0 > \
v gemm(char transA, char transB, ssize_t m, ssize_t n, ssize_t k, ALPHA const* alpha, AAP aa, ssize_t lda, BBP bb, ssize_t ldb, BETA const* beta, CCP cc, ssize_t ldc) \
{ \
	MULTI_MARK_SCOPE("cpu_gemm");			                                                                                                                            \
	using std::max;                                                                                                                                                     \
	if(transA =='N'){MULTI_ASSERT1(lda >= max(1L, m));}else{MULTI_ASSERT1(lda >= max(1L, k));}                                                                           \
	if(transB =='N'){MULTI_ASSERT1(ldb >= max(1L, k));}else{MULTI_ASSERT1(ldb >= max(1L, n));}                                                                           \
	MULTI_ASSERT1( aa != cc );                                                                                                                                                 \
	MULTI_ASSERT1( bb != cc );                                                                                                                                                 \
	MULTI_ASSERT1(ldc >= max(ssize_t{1}, m));                                                                                                                                          \
	if(*beta != 0.) MULTI_ASSERT1((is_assignable<CC&, decltype(ALPHA{}*AA{}*BB{} + BETA{}*CC{})>{}));                                                                          \
	BLAS(T##gemm)(transA, transB, BC(m), BC(n), BC(k), *(T const*)alpha, (T const*)static_cast<AA*>(aa), BC(lda), (T const*)static_cast<BB*>(bb), BC(ldb), *(T const*)beta, (T*)static_cast<CC*>(cc), BC(ldc));               \
}

xgemm(s) xgemm(d) xgemm(c) xgemm(z) // NOLINT(readability-function-cognitive-complexity) : 36 of 25
#undef xgemm

#define xtrsm(T) \
template<class ALPHA, class AAP, class AA = typename pointer_traits<AAP>::element_type, class BBP, class BB = typename pointer_traits<BBP>::element_type, \
enable_if_t< \
	is_##T<AA>{} and is_##T<BB>{} and is_assignable<BB&, decltype(AA{}*BB{}/ALPHA{})>{} and is_assignable<BB&, decltype(ALPHA{}*BB{}/AA{})>{} and \
	is_convertible_v<AAP, AA*> and is_convertible_v<BBP, BB*> \
,int> =0> \
v trsm(char side, char ul, char transA, char diag, ssize_t m, ssize_t n, ALPHA alpha, AAP aa, ssize_t lda, BBP bb, ssize_t ldb){ \
	MULTI_MARK_SCOPE("cpu_trsm");											\
	assert( side   == 'L' or side   == 'R' );                  /* NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay)*/  \
	assert( ul     == 'U' or ul     == 'L' );                  /* NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay)*/  \
	assert( transA == 'N' or transA == 'T' or transA == 'C' ); /* NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay)*/  \
	assert( diag     == 'U' or diag     == 'N' );              /* NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay)*/  \
	MULTI_ASSERT1( m >= 0 and n >= 0 ); \
	using std::max; \
	if(side == 'L'){MULTI_ASSERT1(lda >= max(ssize_t{1}, m));}else if(side == 'R'){assert( lda >= max(ssize_t{1}, n) );}   /* NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay,hicpp-no-array-decay)*/  \
	MULTI_ASSERT1( ldb >= max(ssize_t{1}, m) ); \
	BLAS(T##trsm)(side, ul, transA, diag, BC(m), BC(n), alpha, (T const*)static_cast<AA*>(aa), BC(lda), (T*)static_cast<BB*>(bb), BC(ldb)); \
}
xtrsm(s) xtrsm(d) xtrsm(c) xtrsm(z) // NOLINT(readability-function-cognitive-complexity) : 29 of 25
#undef xtrsm

xsyrk(s) xsyrk(d) xsyrk(c) xsyrk(z)
	              xherk(c) xherk(z)

} // end namespace core

#undef xsyrk
#undef xherk
#undef xtrsm

#undef BC

struct context{ // stateless (and thread safe)
	template<class... As>
	static auto axpy(As... as)
	->decltype(core::axpy(as...)){
		return core::axpy(as...);}

	template<class... As>
	static auto gemv(As... as)
	->decltype(core::gemv(as...)){
		return core::gemv(as...);}

	template<class... As>
	static auto gemm(As&&... as)
	->decltype(core::gemm(std::forward<As>(as)...)){
		return core::gemm(std::forward<As>(as)...);}

	template<class... As>
	static auto dot(As&&... as)
	->decltype(core::dot(std::forward<As>(as)...)){
		return core::dot(std::forward<As>(as)...);}

	template<class... As>
	static auto dotc(As&&... as)
	->decltype(core::dotc(std::forward<As>(as)...)){
		return core::dotc(std::forward<As>(as)...);}

	template<class... As>
	static auto dotu(As&&... as)
	->decltype(core::dotu(std::forward<As>(as)...)){
		return core::dotu(std::forward<As>(as)...);}

	template<class... As>
	static auto trsm(As&&... as)
	->decltype(core::trsm(std::forward<As>(as)...)){
		return core::trsm(std::forward<As>(as)...);}

	template<class... As>
	static auto herk(As&&... as)
	->decltype(core::herk(std::forward<As>(as)...)){
		return core::herk(std::forward<As>(as)...);}
};

template<class Context> struct is_context : std::false_type{};
template<> struct is_context<context> : std::true_type{};
template<> struct is_context<context&&> : std::true_type{};
template<> struct is_context<context&> : std::true_type{};
template<> struct is_context<context const&> : std::true_type{};

template<> struct is_context<void*&> : std::true_type{};

namespace core{
template<class Context, class... As>
auto copy(Context&& /*unused*/, As... as)
->decltype(core::copy(as...)){
	return core::copy(as...);}
} // end namespace core

template<class TPtr, std::enable_if_t<std::is_convertible<TPtr, typename std::pointer_traits<TPtr>::element_type*>{}, int> =0> 
auto default_context_of(TPtr const& /*unused*/) -> blas::context*{return {};}

} // end namespace blas

} // end namespace multi
} // end namespace boost

#endif


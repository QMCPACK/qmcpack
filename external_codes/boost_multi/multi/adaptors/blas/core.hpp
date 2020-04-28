#ifdef COMPILATION_INSTRUCTIONS
(echo '#include"'$0'"'>$0.cpp)&&$CXX -Wall -Wextra -Wpedantic -Wfatal-errors -D_TEST_MULTI_ADAPTORS_BLAS_CORE $0.cpp -o $0x `pkg-config --libs blas`&&$0x&&(rm $0x $0.cpp; for a in `find tests/ -name '*.cpp'`; do sh $a || break; done); exit
#endif
// https://software.intel.com/en-us/articles/intel-mkl-link-line-advisor
// © Alfredo A. Correa 2019

#ifndef MULTI_ADAPTORS_BLAS_CORE_HPP
#define MULTI_ADAPTORS_BLAS_CORE_HPP

//#include <cblas/cblas.h>

#include<iostream> // debug
#include<cassert>
#include<complex>
#include<cstdint> // int64_t
#include<limits> // numeric_limits

#ifdef CBLAS_H
#define BLAS(NamE) cblas_##NamE
#else
#define BLAS(NamE) NamE##_
extern "C"{

#ifndef _BLAS_INT
#define _BLAS_INT __INTPTR_WIDTH__
#endif

#define s float
#define d double
#define c std::complex<s>
#define z std::complex<d>
#define v void 
#define C _Complex s
#define Z _Complex d
#if(_BLAS_INT==32)
#define INT std::int32_t
#endif
#if(_BLAS_INT==64)
#define INT std::int64_t
#endif
#define INTEGER INT const&
#define N INTEGER n
#define INCX INTEGER incx
#define INCY INTEGER incy

static_assert(sizeof(INT)==32/8 or sizeof(INT)==64/8, "please set _BLAS_INT to int32_t or int64_t");

#define xROTG(T1, T2)     v BLAS(   T1##rotg)(T1 const*, T1 const*, T2*, T1*)
#define xROTMG(T)         v BLAS(   T##rotmg)(T*, T*, T*, T const&, T(&param)[5])
#define xROT(TT, T, S)    v BLAS(  TT##rot  )(N,              T       *x, INCX, T       *y, INCY, S const&, S const&)
#define xROTM(T)          v BLAS(   T##rotm )(N, T* x, INCX, T* y, INCY, T const(&p)[5])
#define xSWAP(T)          v BLAS(   T##swap )(N,              T       *x, INCX, T       *y, INCY)
#define xSCAL(TT, TA, TX) v BLAS(  TT##scal )(N, TA const& a, TX      *x, INCX                  )
#define xCOPY(T)          v BLAS(   T##copy )(N,              T const *x, INCX, T       *y, INCY) 
#define xAXPY(T)          v BLAS(   T##axpy )(N,  T const& a, T const *x, INCX, T       *y, INCY)
#define xDOT(R, TT, T)    R BLAS(  TT##dot  )(N,              T const *x, INCX, T const *y, INCY)
#define xDOTU(R, T)       R BLAS(   T##dotu )(N,              T const *x, INCX, T const *y, INCY)
#define xDOTC(R, T)       R BLAS(   T##dotc )(N,              T const *x, INCX, T const *y, INCY)
#define xxDOT(TT, T)      T BLAS(  TT##dot  )(N,  T const& a, T const *x, INCX, T const *y, INCY)
#define xNRM2(R, TT, T)   R BLAS(  TT##nrm2 )(N,              T const *x, INCX                  )   
#define xASUM(R, TT, T)   R BLAS(  TT##asum )(N,              T const *x, INCX                  )
#define IxAMAX(T)       INT BLAS(i##T##amax )(N,              T const* x, INCX                  )

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

#define xGEMV(T) void BLAS(T##gemv)(TRANS, NR, NC, T const& a, T const* A, LDA, T const* X, INCX, T const& beta, T*       Y, INCY)
#define xGER(T)  void BLAS(T##ger )(       NR, NC, T const& a,                  T const* X, INCX,                T const* Y, INCY, T* A, LDA)
#define xGERU(T) void BLAS(T##geru)(       NR, NC, T const& a,                  T const* X, INCX,                T const* Y, INCY, T* A, LDA)
#define xGERC(T) void BLAS(T##gerc)(       NR, NC, T const& a,                  T const* X, INCX,                T const* Y, INCY, T* A, LDA)
#define xTRSV(T) void BLAS(T##trsv)(UPLO, TRANS, DIAG, N, T const* A, LDA, T* X, INCX)

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

#define xGEMM(T)     void BLAS(T##gemm)(TRANSA, TRANSB, NR, NC, NK, T const& a, T const* A, LDA, T const* B, LDB, T const& b, T const* CC, LDC)
#define xSYRK(T)     void BLAS(T##syrk)(UPLO, TRANSA, NR, NK, T const& a, T const* A, LDA, T const& b, T* CC, LDC) 
#define xHERK(TT, T) void BLAS(T##herk)(UPLO, TRANSA, NR, NK, TT const& a, T const* A, LDA, TT const& b, T* CC, LDC) 
#define xTRSM(T) void BLAS(T##trsm)(SIDE, UPLO, TRANSA, DIAG, NR, NK, T const& a, T const* A, LDA, T const* B, LDB) 
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

using s = float;
using d = double;
using c = std::complex<s>;
using z = std::complex<d>;
using v = void;

#define BC(x) [](auto xx){assert(xx>=std::numeric_limits<INT>::min() and xx<std::numeric_limits<INT>::max()); return xx;}(x)

#define xrotg(T1, T2)                       v   rotg (T1 const& a, T1 const& b, T2& cc, T1& ss                                   ){BLAS(T1##rotg )(const_cast<T1*>(&a), const_cast<T1*>(&b), &cc, &ss);}
#define xrotmg(T)                           v   rotmg(T& d1, T& d2, T& A, T const& B, T(&p)[5]                                   ){BLAS(T##rotmg )(&d1, &d2, &A, B, p);}
#define xrot(T, TT, CS)   template<class S> v   rot  (S n,       T       *x, S incx, T       *y, S incy, CS const& c, CS const& s){BLAS(TT##rot )(BC(n),    x, BC(incx), y, BC(incy), c, s);}
#define xrotm(T)          template<class S> v   rotm (S n,       T       *x, S incx, T       *y, S incy, T const(&p)[5]          ){BLAS( T##rotm)(BC(n),    x, BC(incx), y, BC(incy), p);              }
#define xswap(T)          template<class S> v   swap (S n,       T       *x, S incx, T       *y, S incy                          ){BLAS( T##swap)(BC(n),    x, BC(incx), y, BC(incy));                 }
#define xscal(XX, TA, TX) template<class S> TX* scal (S n, TA* a, TX      *x, S incx                                              ){BLAS(XX##scal)(BC(n), *a, x, BC(incx)             ); return x+n*incx;}
#define xcopy(T)          template<class S> v   copy (S n,       T const *x, S incx, T       *y, S incy                          ){BLAS( T##copy)(BC(n),    x, BC(incx), y, BC(incy));                 }
#define xaxpy(T)          template<class S> T*  axpy (S n, T  a, T const *x, S incx, T       *y, S incy                          ){BLAS( T##axpy)(BC(n), a, x, BC(incx), y, BC(incy)); return y+n*incy;}
#define xdot(R, TT, T)    template<class S> v   dot  (S n,       T const *x, S incx, T const *y, S incy, R* r                    ){*r = BLAS(TT##dot )(BC(n),    x, BC(incx), y, BC(incy));                 }

xrotg(s, s)    xrotg(d, d) //MKL extension xrotg(c, s); xrotg(z, d);
xrotmg(s)      xrotmg(d)
xrot(s, s, s)  xrot(d, d, d)  xrot(c, cs, s) xrot(z, zd, d)
xrotm(s)       xrotm(d)
xswap(s)       xswap(d)       xswap(c)       xswap(z)

namespace core{
xscal(s, s, s) xscal(d, d, d) xscal(c, c, c) xscal(z, z, z) xscal(zd, d, z) xscal(cs, s, c)
xcopy(s)       xcopy(d)       xcopy(c)       xcopy(z)

xdot(s, s, s)  xdot(d, d, d)                                xdot(d, ds, s)
}
xaxpy(s)       xaxpy(d)       xaxpy(c)       xaxpy(z)

template<class R, class S, class T> R dot(S n, T const* x, S incx, T const* y, S incy){
	R ret;
	dot(n, x, incx, y, incy, &ret);
	return ret;
}

template<class S, class T> T dot(S n, T const* x, S incx, T const* y, S incy){
	return dot<T, S, T>(n, x, incx, y, incy);
}

#undef xrotg
#undef xrot
#undef xswap
#undef xscal
#undef xcopy
#undef xaxpy
#undef xdot

#ifndef CBLAS_H
//#define xdotu(T) template<class S> v dotu(S n, T const* x, S incx, T const* y, S incy, T* r){*r = BLAS(T##dotu)(BC(n), x, BC(incx), y, BC(incy));}
#define xdotu(T) template<class S> v dotu(S n, T const* x, S incx, T const* y, S incy, T* r){*r = (T)(BLAS(T##dotu)(BC(n), x, BC(incx), y, BC(incy)));}
#define xdotc(T) template<class S> v dotc(S n, T const* x, S incx, T const* y, S incy, T* r){*r = (T)(BLAS(T##dotc)(BC(n), x, BC(incx), y, BC(incy)));}
namespace core{
xdotu(c) xdotu(z)
xdotc(c) xdotc(z)
}
//                 template<class S> z dot(S n,  c const *x, S incx, c const *y, S incy){return dotc(n, x, incx, y, incy);}
//                 template<class S> z dot(S n,  z const *x, S incx, z const *y, S incy){return dotc(n, x, incx, y, incy);}

#undef xdotu
#undef xdotc
#else
#define xdotu(T) template<class S> v dotu(S n, T const* x, S incx, T const* y, S incy, T* r){BLAS(T##dotu_sub)(BC(n), x, BC(incx), y, BC(incy), r);}
#define xdotc(T) template<class S> v dotc(S n, T const* x, S incx, T const* y, S incy, T* r){BLAS(T##dotc_sub)(BC(n), x, BC(incx), y, BC(incy), r);}
namespace core{
xdotu(c) xdotu(z)
xdotc(c) xdotc(z)
}
#undef xdotu
#undef xdotc
#endif

namespace core{
template<class S> s dot(S n, s const& b, s const* x, S incx, s const* y, S incy){return BLAS(sdsdot)(BC(n), b, x, BC(incx), y, BC(incy));}

//template<class S> void dot(S n, s const& b, s const* x, S incx, s const* y, S incy, s* result){*result = BLAS(sdsdot)(BC(n), b, x, BC(incx), y, BC(incy));}

}

#define xnrm2(R, T, TT) template<class S>    v nrm2 (S n, T const* x, S incx, R* r){*r = BLAS(TT##nrm2  )(BC(n), x, BC(incx));}
#define xasum(T, TT)    template<class S> auto asum (S n, T const* x, S incx){return BLAS(TT##asum  )(BC(n), x, BC(incx));}
#define ixamax(T)       template<class S> auto iamax(S n, T const* x, S incx){return BLAS(i##T##amax)(BC(n), x, BC(incx)) - 1;}
xasum(s, s)    xasum(d, d)                        xasum (c, sc)                  xasum(z, dz)
namespace core{
	xnrm2(s, s, s) xnrm2(d, d, d)  xnrm2(s, c, sc) xnrm2(d, z, dz)
	ixamax(s)      ixamax(d)       ixamax(c)       ixamax(z)
}
#undef xnrm2
#undef xasum
#undef ixamax


///////////////////////////////////////////////////////////////////////////////
// LEVEL2
#define xgemv(T) template<class C, class S> v gemv(C trans, S m, S n, T const& a, T const* A, S lda, T const* X, S incx, T beta, T*       Y, S incy             ){BLAS(T##gemv)(trans, BC(m), BC(n), a, A, BC(lda), X, BC(incx), beta, Y, BC(incy)            );}
#define xger(T)  template<         class S> v ger(          S m, S n, T const& a,                    T const* X, S incx,         T const* Y, S incy, T* A, S lda){BLAS(T##ger )(       BC(m), BC(n), a,             X, BC(incx),       Y, BC(incy), A, BC(lda));}
                 template<         class S> v ger(          S m, S n, c const& a,                    c const* X, S incx,         c const* Y, S incy, c* A, S lda){BLAS(cgeru  )(       BC(m), BC(n), a,             X, BC(incx),       Y, BC(incy), A, BC(lda));}
                 template<         class S> v ger(          S m, S n, z const& a,                    z const* X, S incx,         z const* Y, S incy, z* A, S lda){BLAS(zgeru  )(       BC(m), BC(n), a,             X, BC(incx),       Y, BC(incy), A, BC(lda));}
#define xgeru(T) template<         class S> v geru(         S m, S n, T const& a,                    T const* X, S incx,         T const* Y, S incy, T* A, S lda){BLAS(T##geru)(       BC(m), BC(n), a,             X, BC(incx),       Y, BC(incy), A, BC(lda));}
#define xgerc(T) template<         class S> v gerc(         S m, S n, T const& a,                    T const* X, S incx,         T const* Y, S incy, T* A, S lda){BLAS(T##gerc)(       BC(m), BC(n), a,             X, BC(incx),       Y, BC(incy), A, BC(lda));}
xgemv(s) xgemv(d) xgemv(c) xgemv(z)
xger(s)   xger(d)
                  xgeru(c) xgeru(z)
                  xgerc(c) xgerc(z)

template<class T> 
struct blas2{
//	template<class S>
//	static v trsv(char ulA, char transA, char di, S m, T const* A, S lda, T* X, S incx) = delete;
};

template<> struct blas2<s>{template<class... As> static v trsv(As... as){BLAS(strsv)(as...);}};
template<> struct blas2<d>{template<class... As> static v trsv(As... as){BLAS(dtrsv)(as...);}};
template<> struct blas2<c>{template<class... As> static v trsv(As... as){BLAS(ctrsv)(as...);}};
template<> struct blas2<z>{template<class... As> static auto trsv(As... as)->decltype(BLAS(ztrsv)(as...)){BLAS(ztrsv)(as...);}};

namespace core{
template<typename TconstP, typename TP, typename S=std::size_t, typename C=char> v trsv(C ulA, C transA, C diA, S n, TconstP A, S lda, TP X, S incx){blas2<std::decay_t<typename std::pointer_traits<TP>::element_type>>::trsv(ulA, transA, diA, n, A, lda, X, incx);}
}

#undef xgemv
#undef xger
#undef xgeru
#undef xgerc

///////////////////////////////////////////////////////////////////////////////
// LEVEL 3
#define xgemm(T) \
template<class C, class S> v gemm(C transA, C transB, S m, S n, S k, T const* a, T const* A, S lda, T const* B, S ldb, T const* beta, T* CC, S ldc){ \
	if(transA == 'N' or transA == 'n') {assert(lda >= std::max(1l, m));} else {assert(lda >= std::max(1l, k));} \
	assert(ldb >= std::max(1l, transB=='N'?k:n)); \
	/*if(transB == 'N' or transB == 'n') {std::cerr<<"ldb,k,n,m ="<< ldb <<','<<k<<','<<n<<','<<m<<std::endl;*/ \
		/*if(ldb==1 and n==1) return gemv(transA, m, k, a, A, lda, B, S{1}, beta, CC, S{1});*/ \
		/*assert(ldb >= std::max(1l, k));} else {std::cerr<<"ldb ="<< ldb <<std::endl; assert(ldb >= std::max(1l, n));}*/ \
	assert(ldc >= std::max(1l, m)); \
	BLAS(T##gemm)(transA, transB, BC(m), BC(n), BC(k), *a, A, BC(lda), B, BC(ldb), *beta, CC, BC(ldc));\
}
#define xsyrk(T) template<class UL, class C, class S> v syrk(UL ul, C transA, S n, S k, T alpha, T const* A, S lda, T beta, T* CC, S ldc){BLAS(T##syrk)(ul, transA, BC(n), BC(k), alpha, A, BC(lda), beta, CC, BC(ldc));}
#define xherk(T) template<class UL, class C, class S, class Real> v herk(UL ul, C transA, S n, S k, Real alpha, T const* A, S lda, Real beta, T* CC, S ldc){BLAS(T##herk)(ul, transA, BC(n), BC(k), alpha, A, BC(lda), beta, CC, BC(ldc));}
#define xtrsm(T) template<class C, class UL, class Di, class S> v trsm(C side, UL ul, C transA, Di di, S m, S n, T alpha, T const* A, S lda, T* B, S ldb){BLAS(T##trsm)(side, ul, transA, di, BC(m), BC(n), alpha, A, lda, B, ldb);}

namespace core{
xgemm(s) xgemm(d) xgemm(c) xgemm(z)
}

namespace core{
xsyrk(s) xsyrk(d) xsyrk(c) xsyrk(z)
                  xherk(c) xherk(z)
xtrsm(s) xtrsm(d) xtrsm(c) xtrsm(z)
}

#undef xgemm
#undef xsyrk
#undef xherk
#undef xtrsm

#undef BC

}}}
///////////////////////////////////////////////////////////////////////////////

#if _TEST_MULTI_ADAPTORS_BLAS_CORE

#include "../../array.hpp"
#include "../../utility.hpp"

#include<complex>
#include<cassert>
#include<iostream>
#include<numeric>
#include<algorithm>

using std::cout;
namespace multi = boost::multi;

int main(){}

#endif
#endif


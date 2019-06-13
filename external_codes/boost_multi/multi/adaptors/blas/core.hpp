#ifdef COMPILATION_INSTRUCTIONS
(echo "#include\""$0"\"" > $0x.cpp) && clang++ `#-DNDEBUG` -O3 -std=c++14 -Wall -Wextra -Wpedantic -Wfatal-errors -D_TEST_MULTI_ADAPTORS_BLAS_CORE -DADD_ $0x.cpp -o $0x.x -lblas && time $0x.x $@ && rm -f $0x.x $0x.cpp; exit
#endif
// Alfredo A. Correa 2019 ©

#ifndef MULTI_ADAPTORS_BLAS_CORE_HPP
#define MULTI_ADAPTORS_BLAS_CORE_HPP

#include<iostream>
//#include <cblas/cblas.h>

#include<cassert>
#include<complex>
#include "../../utility.hpp"
#include "../../array.hpp" // allocating multi::arrays for output
#include<numeric> // inner_product
//#include<algorithm>
#include<iostream>

#ifdef CBLAS_H
#define BLAS(NamE) cblas_##NamE
#else
#define BLAS(NamE) NamE##_
extern "C"{

#define s float
#define d double
#define c std::complex<s>
#define z std::complex<d>
#define v void 
#define C _Complex s
#define Z _Complex d
#define INT int
#define INTEGER INT const&
#define N INTEGER n
#define INCX INTEGER incx
#define INCY INTEGER incy

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
#define NR const int& nr
#define NC const int& nc
#define LDA const int& lda

#define xGEMV(T) void BLAS(T##gemv)(TRANS, NR, NC, T const& a, T const* A, LDA, T const* X, INCX, T const& beta, T*       Y, INCY)
#define xGER(T)  void BLAS(T##ger )(       NR, NC, T const& a,                  T const* X, INCX,                T const* Y, INCY, T* A, LDA)
#define xGERU(T) void BLAS(T##geru)(       NR, NC, T const& a,                  T const* X, INCX,                T const* Y, INCY, T* A, LDA)
#define xGERC(T) void BLAS(T##gerc)(       NR, NC, T const& a,                  T const* X, INCX,                T const* Y, INCY, T* A, LDA)

xGEMV(s); xGEMV(d); xGEMV(c); xGEMV(z);
xGER(s); xGER(d);
xGERU(c); xGERU(z);
xGERC(c); xGERC(z);

#define TRANSA const char& transa
#define TRANSB const char& transb
#define NK const int& nk
#define LDB const int& ldb
#define LDC const int& ldc

#define xGEMM(T) void BLAS(T##gemm)(TRANSA, TRANSB, NR, NC, NK, T const& a, T const* A, LDA, T const* B, LDB, T const& b, T* CC, LDC)
xGEMM(s); xGEMM(d); xGEMM(c); xGEMM(z);

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

#define xrotg(T1, T2)                       v   rotg (T1 const& a, T1 const& b, T2& cc, T1& ss){BLAS(T1##rotg )(const_cast<T1*>(&a), const_cast<T1*>(&b), &cc, &ss);}
#define xrotmg(T)                           v   rotmg(T& d1, T& d2, T& A, T const& B, T(&p)[5]){BLAS(T##rotmg )(&d1, &d2, &A, B, p);}
#define xrot(T, TT, CS)   template<class S> v   rot  (S n,       T       *x, S incx, T       *y, S incy, CS const& c, CS const& s){       BLAS(TT##rot )(n,    x, incx, y, incy, c, s);}
#define xrotm(T)          template<class S> v   rotm (S n,       T       *x, S incx, T       *y, S incy, T const(&p)[5]){BLAS(T##rotm)(n, x, incx, y, incy, p);}
#define xswap(T)          template<class S> v   swap (S n,       T       *x, S incx, T       *y, S incy){       BLAS( T##swap)(n,    x, incx, y, incy);} 
#define xscal(XX, TA, TX) template<class S> TX* scal (S n, TA a, TX      *x, S incx                    ){       BLAS(XX##scal)(n, a, x, incx         ); return x+n*incx;}
#define xcopy(T)          template<class S> v   copy (S n,       T const *x, S incx, T       *y, S incy){       BLAS( T##copy)(n,    x, incx, y, incy);} 
#define xaxpy(T)          template<class S> T*  axpy (S n, T  a, T const *x, S incx, T       *y, S incy){       BLAS( T##axpy)(n, a, x, incx, y, incy); return y+n*incy;}
#define xdot(R, TT, T)    template<class S> v   dot  (S n,       T const *x, S incx, T const *y, S incy, R& r){r = BLAS(TT##dot)(n,    x, incx, y, incy);               }

xrotg(s, s)    xrotg(d, d) //MKL extension xrotg(c, s); xrotg(z, d);
xrotmg(s)      xrotmg(d)
xrot(s, s, s)  xrot(d, d, d)  xrot(c, cs, s) xrot(z, zd, d)
xrotm(s)       xrotm(d)
xswap(s)       xswap(d)       xswap(c)       xswap(z)
xscal(s, s, s) xscal(d, d, d) xscal(c, c, c) xscal(z, z, z)                xscal(zd, d, z) xscal(cs, s, c)
xcopy(s)       xcopy(d)       xcopy(c)       xcopy(z)
xaxpy(s)       xaxpy(d)       xaxpy(c)       xaxpy(z)
xdot(s, s, s)  xdot(d, d, d)                                xdot(d, ds, s)

template<class R, class S, class T> R dot(S n, T const* x, S incx, T const* y, S incy){
	R ret;
	dot(n, x, incx, y, incy, ret);
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
#define xdotu(T) template<class S> T dotu(S n, T const* x, S incx, T const* y, S incy){return BLAS(T##dotu)(n, x, incx, y, incy);}
#define xdotc(T) template<class S> T dotc(S n, T const* x, S incx, T const* y, S incy){return BLAS(T##dotc)(n, x, incx, y, incy);}
xdotu(c) xdotu(z)
xdotc(c) xdotc(z)

                 template<class S> z dot(S n,  c const *x, S incx, c const *y, S incy){return dotc(n, x, incx, y, incy);}
                 template<class S> z dot(S n,  z const *x, S incx, z const *y, S incy){return dotc(n, x, incx, y, incy);}

#undef xdotu
#undef xdotc
#else
#define xdotu(T) template<class S> T dotu(S n, T const* x, S incx, T const* y, S incy){T ret; BLAS(T##dotu_sub)(n, x, incx, y, incy, &ret); return ret;}
#define xdotc(T) template<class S> T dotc(S n, T const* x, S incx, T const* y, S incy){T ret; BLAS(T##dotc_sub)(n, x, incx, y, incy, &ret); return ret;}
xdotu(c) xdotu(z)
xdotc(c) xdotc(z)
#undef xdotu
#undef xdotc
#endif

template<class S> s dot(S n, s const& b, s const* x, S incx, s const* y, S incy){return BLAS(sdsdot)(n, b, x, incx, y, incy);}

#define xnrm2(R, T, TT) template<class S>    R nrm2 (S n, T const* x, S incx){return BLAS(TT##nrm2  )(n, x, incx);}
#define xasum(T, TT)    template<class S> auto asum (S n, T const* x, S incx){return BLAS(TT##asum  )(n, x, incx);}
#define ixamax(T)       template<class S> auto iamax(S n, T const* x, S incx){return BLAS(i##T##amax)(n, x, incx);}
 xnrm2(s, s, s) xnrm2(d, d, d)                                    xnrm2(s, c, sc)              xnrm2(d, z, dz)
 xasum(s, s)    xasum(d, d)                        xasum (c, sc)                  xasum(z, dz)
ixamax(s)      ixamax(d)       ixamax(c) ixamax(z)
#undef xnrm2
#undef xasum
#undef ixamax

///////////////////////////////////////////////////////////////////////////////
// LEVEL2
#define xgemv(T) template<class C, class S> v gemv(C trans, S m, S n, T const& a, T const* A, S lda, T const* X, S incx, T beta, T*       Y, S incy             ){BLAS(T##gemv)(trans, m, n, a, A, lda, X, incx, beta, Y, incy        );}
#define xger(T)  template<         class S> v ger(          S m, S n, T const& a,                    T const* X, S incx,         T const* Y, S incy, T* A, S lda){BLAS(T##ger )(       m, n, a,         X, incx,       Y, incy, A, lda);}
                 template<         class S> v ger(          S m, S n, c const& a,                    c const* X, S incx,         c const* Y, S incy, c* A, S lda){BLAS(cgeru  )(       m, n, a,         X, incx,       Y, incy, A, lda);}
                 template<         class S> v ger(          S m, S n, z const& a,                    z const* X, S incx,         z const* Y, S incy, z* A, S lda){BLAS(zgeru  )(       m, n, a,         X, incx,       Y, incy, A, lda);}
#define xgeru(T) template<         class S> v geru(         S m, S n, T const& a,                    T const* X, S incx,         T const* Y, S incy, T* A, S lda){BLAS(T##geru)(       m, n, a,         X, incx,       Y, incy, A, lda);}
#define xgerc(T) template<         class S> v gerc(         S m, S n, T const& a,                    T const* X, S incx,         T const* Y, S incy, T* A, S lda){BLAS(T##gerc)(       m, n, a,         X, incx,       Y, incy, A, lda);}
xgemv(s) xgemv(d) xgemv(c) xgemv(z)
xger(s)   xger(d)
                  xgeru(c) xgeru(z)
                  xgerc(c) xgerc(z)
#undef xgemv
#undef xger
#undef xgeru
#undef xgerc

///////////////////////////////////////////////////////////////////////////////
// LEVEL 3
#define xgemm(T) template<class C, class S> v gemm(C transA, C transB, S m, S n, S k, T const& a, T const* A, S lda, T const* B, S ldb, T const& beta, T* CC, S ldc){BLAS(T##gemm)(transA, transB, m, n, k, a, A, lda, B, ldb, beta, CC, ldc);}
xgemm(s) xgemm(d) xgemm(c) xgemm(z)

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


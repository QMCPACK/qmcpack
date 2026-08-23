// Copyright 2019-2024 Alfredo A. Correa
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_MULTI_ADAPTORS_LAPACK_CORE_HPP
#define BOOST_MULTI_ADAPTORS_LAPACK_CORE_HPP

//#include<iostream>
#include<cassert>
#include<complex>

// #include <cblas/cblas.h>
// #include<lapacke.h>
// #include<lapack.h>

#define s float
#define d double
#define c std::complex<s>
#define z std::complex<d>
#define v void 

#define INT int
#define INTEGER INT const&

//#define N INTEGER n
#define CHARACTER char const&
#define UPLO CHARACTER
#define JOBZ CHARACTER
#define LAPACK(NamE) NamE##_
#define LWORK INTEGER lwork
#define LIWORK INTEGER liwork
#define IWORK int*

// cppcheck-suppress [preprocessorErrorDirective] bug in cppcheck 2.11
#define xPOTRF(T)     v LAPACK(T##potrf)(UPLO, int const& N, T*, int const& LDA, int& INFO)  // NOLINT(bugprone-macro-parentheses,readability-identifier-length)
#define xSYEV(T)      v LAPACK(T##syev) (JOBZ, UPLO, int const& N, T*, int const& LDA, T*, T*, LWORK, int& INFO)  // NOLINT(bugprone-macro-parentheses)
#define xSYEVD(T)     v LAPACK(T##syevd)(JOBZ, UPLO, int const& N, T*, int const& LDA, T*, T*, LWORK, IWORK, LIWORK, int& INFO)  // NOLINT(bugprone-macro-parentheses)
#define xHEEV(T)      v LAPACK(T##heev) (JOBZ, UPLO, int const& N, T*, int const& LDA, T*, T*, LWORK, int& INFO)  // NOLINT(bugprone-macro-parentheses)

#define subroutine  void
#define integer     int const&
#define integer_out int&
#define integer_ptr int*
#define integer_cptr int const*
#define character   char const&

// http://www.netlib.org/lapack/explore-html/dd/d9a/group__double_g_ecomputational_ga0019443faea08275ca60a734d0593e60.html
#define xGETRF(T)     \
/*NOLINTBEGIN(readability-identifier-length,bugprone-macro-parentheses)*/ \
void T##getrf_( \
	integer     M,    /*The number of rows of the matrix A.  M >= 0.*/           \
	integer     N,    /*The number of columns of the matrix A.  N >= 0.*/        \
	T*          A,    /*On entry, the M-by-N matrix to be factored.*/      \
                      /*On exit, the factors L and U from the factorization*/    \
	integer     LDA,  /*The leading dimension of the array A.  LDA >= max(1,M).*/\
	integer_ptr IPIV, /*The pivot indices; for 1 <= i <= min(M,N), row i of the matrix was interchanged with row IPIV(i).*/\
	integer_out INFO  /*= 0:  successful exit*/\
                      /*< 0:  if INFO = -i, the i-th argument had an illegal value*/\
                      /*> 0:  if INFO = i, U(i,i) is exactly zero. The factorization has been completed, but the factor U is exactly singular, and division by zero will occur if it is used to solve a system of equations.*/\
) \
/*NOLINTEND(readability-identifier-length,bugprone-macro-parentheses)*/

// NOLINTBEGIN(bugprone-macro-parentheses)
// NOLINTBEGIN(readability-identifier-length)
// http://www.netlib.org/lapack/explore-html/d8/ddc/group__real_g_ecomputational_gaa00bcf4d83a118cb6f0b6619d6ffaa24.html
#define xGETRS(T)     \
void T##getrs_( \
	character    TRANS,/*Specifies the form of the system of equations:             */\
	                   /* = 'N':  A * X = B  (No transpose)                         */\
	                   /* = 'T':  A**T* X = B  (Transpose)                          */\
	                   /* = 'C':  A**T* X = B  (Conjugate transpose = Transpose)    */\
	integer      N,    /*The order of the matrix A.  N >= 0.                        */\
	integer      NRHS, /*The number of right hand sides, i.e., the number of columns*/\
	                   /*of the matrix B.  NRHS >= 0.                               */\
	T const*     A,    /* The factors L and U from the factorization A = P*L*U      */\
	                   /*as computed by SGETRF.                                     */\
	integer      LDA,  /*The leading dimension of the array A.  LDA >= max(1,N).    */\
	integer_cptr IPIV, /*The pivot indices from SGETRF; for 1<=i<=N, row i of the   */\
	                   /*matrix was interchanged with row IPIV(i).                  */\
	T*           B,    /*On entry, the right hand side matrix B.                    */\
	                   /*On exit, the solution matrix X.                            */\
	integer      LDB,  /*The leading dimension of the array B.  LDB >= max(1,N).    */\
	integer_out  INFO  /*= 0:  successful exit                                      */\
	                   /*< 0:  if INFO = -i, the i-th argument had an illegal value */\
)
// NOLINTEND(readability-identifier-length)
// NOLINTEND(bugprone-macro-parentheses)

// TODO(correaa) // http://www.netlib.org/lapack/explore-html/d7/d3b/group__double_g_esolve_ga5ee879032a8365897c3ba91e3dc8d512.html

extern "C"{
xGETRF(s)   ; xGETRF(d)   ; xGETRF(c)   ; xGETRF(z)   ;
xGETRS(s)   ; xGETRS(d)   ; xGETRS(c)   ; xGETRS(z)   ;
}

namespace core {
// http://www.netlib.org/lapack/explore-html/da/d30/a18643_ga5b625680e6251feb29e386193914981c.html

using lapack_int = int;

// TODO(correaa) make into a template, then remove inline
inline auto getrf(lapack_int m, lapack_int n, double* A, lapack_int lda, int* ipiv) -> int {  // NOLINT(readability-identifier-length) lapack conventional name
	assert( m >= 0 );
	assert( n >= 0 );
	assert( lda >= std::max(lapack_int{1}, m) );
	int info;  // NOLINT(cppcoreguidelines-init-variables) delayed initialization
	dgetrf_(m, n, A, lda, ipiv, info);
	assert(info >= 0);
	return info;
}

// TODO(correaa) make into a template, then remove inline
inline void getrs(char trans, lapack_int const n, lapack_int const nrhs, double const* A, lapack_int const lda, int const* ipiv, double* B, lapack_int const ldb) {  // NOLINT(readability-identifier-length) lapack conventional name
	assert( trans == 'T' || trans == 'N' || trans == 'C' );
	assert( n >= 0 );
	assert( nrhs >= 0 );
	assert( lda >= std::max(1, n) );
	int info;  // NOLINT(cppcoreguidelines-init-variables) delayed initialization
	dgetrs_(trans, n, nrhs, A, lda, ipiv, B, ldb, info);

	switch(info){
		case -1: throw std::logic_error{"transa != 'N', 'T', or 'C'"};
		case -2: throw std::logic_error{"n < 0"                    };
		case -3: throw std::logic_error{"nrhs < 0"                 };
		case -4: throw std::logic_error{"n > lda"                  };
		case -5: throw std::logic_error{"lda ≤ 0"                  };
		case -6: throw std::logic_error{"n > ldb"                  };
		case -7: throw std::logic_error{"ldb ≤ 0"                  };
		case -8: throw std::logic_error{"error!"                  };
		default: assert(info == 0);
	}
}

}  // end namespace core

namespace lapack {

struct context{
	template<class... Args> static auto getrf(Args&&... args)->decltype(core::getrf(args...)){return core::getrf(std::forward<Args>(args)...);}
	template<class... Args> static auto getrs(Args&&... args)->decltype(core::getrs(args...)){return core::getrs(std::forward<Args>(args)...);}
};

}  // end namespace lapack

extern "C" {
xPOTRF(s)   ; xPOTRF(d)    ;
xPOTRF(c)   ; xPOTRF(z)    ;

//xSYEV(s)    ; xSYEV(d)     ;
//xSYEVD(s)   ; xSYEVD(d)    ;
//xHEEV(c)    ; xHEEV(z)     ;
}

#undef subroutine
#undef integer
#undef character

#undef JOBZ
#undef UPLO
#undef INFO
#undef CHARACTER
#undef N
#undef LDA

#undef INTEGER
#undef INT


#define xpotrf(T) template<class S> v potrf(char uplo, S n, T *x, S incx, int& info){LAPACK(T##potrf)(uplo, n, x, incx, info);}  // NOLINT(bugprone-macro-parentheses,readability-identifier-length)

namespace core {
xpotrf(s) xpotrf(d)
xpotrf(c) xpotrf(z)
}  //end namespace core

// NOLINTBEGIN(bugprone-macro-parentheses)

// http://www.netlib.org/lapack/explore-html/d2/d8a/group__double_s_yeigen_ga442c43fca5493590f8f26cf42fed4044.html
#define xsyev(T) template<class S> v syev(char jobz, char uplo, S n, T* a, S lda, T* w, T* work, S lwork, int& info){LAPACK(T##syev)(jobz, uplo, n, a, lda, w, work, lwork, info);}

// http://www.netlib.org/lapack/explore-html/d2/d8a/group__double_s_yeigen_ga77dfa610458b6c9bd7db52533bfd53a1.html
#define xsyevd(T) template<class S> v syevd(char jobz, char uplo, S n, T* a, S lda, T* w, T* work, S lwork, int* iwork, S liwork, int& info){ \
	if(n <= 1               ){assert(lwork >= 1              ); assert(liwork >=1       );} \
	if(jobz == 'N' and n > 1){assert(lwork >= 2*n+1          ); assert(liwork >= 1      );} \
	if(jobz == 'V' and n > 1){assert(lwork >= 1 + 6*n + 2*n*n); assert(liwork >= 3 + 5*n);} \
	LAPACK(T##syevd)(jobz, uplo, n, a, lda, w, work, lwork, iwork, liwork, info); \
}

#define xheev(T) template<class S> v heev(char jobz, char uplo, S n, T* a, S lda, T* w, T* work, S lwork, int& info){LAPACK(T##heev)(jobz, uplo, n, a, lda, w, work, lwork, info);}
// NOLINTEND(bugprone-macro-parentheses)

// namespace core{
// // xsyev (s) xsyev (d)
// // xsyevd(s) xsyevd(d)
// //                     xheev(c) xheev(z)
// }

#undef s
#undef d
#undef c
#undef z
#undef v

#define TRANS const char& trans

#endif

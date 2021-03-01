#ifdef COMPILATION_INSTRUCTIONS
(echo "#include\""$0"\"" > $0x.cpp) && clang++ `#-DNDEBUG` -O3 -std=c++14 -Wall -Wextra -Wpedantic -Wfatal-errors -D_TEST_MULTI_ADAPTORS_LAPACK_CORE -DADD_ $0x.cpp -o $0x.x -lblas -llapack && time $0x.x $@ && rm -f $0x.x $0x.cpp; exit
#endif
// Alfredo A. Correa 2019 ©

#ifndef MULTI_ADAPTORS_LAPACK_CORE_HPP
#define MULTI_ADAPTORS_LAPACK_CORE_HPP

//#include<iostream>
#include<cassert>
#include<complex>

//#include <cblas/cblas.h>

#define s float
#define d double
#define c std::complex<s>
#define z std::complex<d>
#define v void 

#define INT int
#define INTEGER INT const&

#define LDA const int& lda
#define N INTEGER n
#define CHARACTER char const&
#define INFO int& info
#define UPLO CHARACTER
#define JOBZ CHARACTER
#define LAPACK(NamE) NamE##_
#define LWORK INTEGER lwork
#define LIWORK INTEGER liwork
#define IWORK int*

#define xPOTRF(T)     v LAPACK(T##potrf)(UPLO, N, T*, LDA, INFO)
#define xSYEV(T)      v LAPACK(T##syev)(JOBZ, UPLO, N, T*, LDA, T*, T*, LWORK, INFO)
#define xSYEVD(T)     v LAPACK(T##syevd)(JOBZ, UPLO, N, T*, LDA, T*, T*, LWORK, IWORK, LIWORK, INFO)
#define xHEEV(T)      v LAPACK(T##heev)(JOBZ, UPLO, N, T*, LDA, T*, T*, LWORK, INFO)

extern "C"{
xPOTRF(s)   ; xPOTRF(d)    ;
xPOTRF(c)   ; xPOTRF(z)    ;

xSYEV(s)    ; xSYEV(d)     ;
xSYEVD(s)   ; xSYEVD(d)    ;
                             xHEEV(c)    ; xHEEV(z)     ;
}

#undef JOBZ
#undef UPLO
#undef INFO
#undef CHARACTER
#undef N
#undef LDA

#undef INTEGER
#undef INT


#define xpotrf(T) template<class S> v potrf(char uplo, S n, T *x, S incx, int& info){LAPACK(T##potrf)(uplo, n, x, incx, info);}

namespace core{
xpotrf(s) xpotrf(d)
xpotrf(c) xpotrf(z)
}

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

namespace core{
xsyev (s) xsyev (d)
xsyevd(s) xsyevd(d)
                    xheev(c) xheev(z)
}

#undef s
#undef d
#undef c
#undef z
#undef v

#define TRANS const char& trans

///////////////////////////////////////////////////////////////////////////////

#if _TEST_MULTI_ADAPTORS_LAPACK_CORE

#include "../../array.hpp"
#include "../../utility.hpp"

#include<iostream>
#include<numeric>
#include<vector>

namespace multi = boost::multi;
using std::cout; 

int main(){
	using core::potrf;

	std::vector<double> v = {
		2., 1.,
		1., 2.
	};
	cout 
		<< v[0] <<'\t'<< v[1] <<'\n'
		<< v[2] <<'\t'<< v[3] <<'\n' << std::endl
	;
	int info;
	potrf('U', 2, v.data(), 2, info);
	cout << "error " << info << std::endl;
	cout 
		<< v[0] <<'\t'<< v[1] <<'\n'
		<< v[2] <<'\t'<< v[3] <<'\n'
	;
	cout << std::endl;
}

#endif
#endif


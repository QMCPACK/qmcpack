#ifdef COMPILATION_INSTRUCTIONS
(echo "#include\""$0"\"" > $0x.cpp) && clang++ `#-DNDEBUG` -O3 -std=c++14 -Wall -Wextra -Wpedantic -Wfatal-errors -D_TEST_MULTI_ADAPTORS_BLAS_GEMV -DADD_ $0x.cpp -o $0x.x -lblas && time $0x.x $@ && rm -f $0x.x $0x.cpp; exit
#endif
// Alfredo A. Correa 2019 ©

#ifndef MULTI_ADAPTORS_BLAS_GEMV_HPP
#define MULTI_ADAPTORS_BLAS_GEMV_HPP

#include "../blas/core.hpp"

namespace boost{
namespace multi{
namespace blas{

struct op{enum : char{N='N', T='T', C='C'};};

struct conj{template<class T> auto operator()(T const& t) const{using std::conj; return conj(t);}};

template<class Op, class AA, class BB, class It1, class Size1, class It2, class Size2, class Out>
auto gemm_n(
	Op opA, Op opB, AA const& a, 
	It1 Af, Size1 An, It2 Bf, Size2 Bn, BB const& b, 
	Out Cf
){
	assert( Af->stride() == 1 );
	assert( Bf->stride() == 1 );
	assert( Cf->stride() == 1 );
	switch(opA){
	case op::N: assert( Cf->size() == Af->size() );
		switch(opB){
		case op::N:
			gemm(opA, opB, Cf->size(), Bn, An, a, base(Af), stride(Af), base(Bf), stride(Bf), b, base(Cf), stride(Cf));
			return Cf + Bn;
		case op::T: case op::C:
			gemm(opA, opB, Cf->size(), Bn, An, a, base(Af), stride(Af), base(Bf), stride(Bf), b, base(Cf), stride(Cf));
			return Cf + Bf->size();
		}
	case op::T: case op::C: assert( Cf->size() == An );
		switch(opB){
		case op::N:
			gemm(opA, opB, Bn, An, Cf->size(), a, base(Af), stride(Af), base(Bf), stride(Bf), b, base(Cf), stride(Cf));
			return Cf + Bn;
		case op::T: case op::C:
			gemm(opA, opB, Cf->size(), An, Bn, a, base(Af), stride(Af), base(Bf), stride(Bf), b, base(Cf), stride(Cf));
			return Cf + Bf->size();
		}
	}
	assert(0);
	return Cf;
}

template<class Op, class AA, class BB, class It1, class It2, class Out>
auto gemm(Op opA, Op opB, AA a, It1 Af, It1 Al, It2 Bf, It2 Bl, BB b, Out Cf){
	assert( stride(Af)==stride(Al) and Af->size()==Al->size() );
	assert( stride(Bf)==stride(Bl) and Bf->size()==Bl->size() );
	return gemm_n(opA, opB, a, Af, Al - Af, Bf, Bl - Bf, b, Cf);
}

template<class Op, class AA, class BB, class A2D, class B2D, class C2D>
C2D&& gemm(Op opA, Op opB, AA a, A2D const& A, B2D const& B, BB b, C2D&& C){
	auto e = gemm(opA, opB, a, begin(A), end(A), begin(B), end(B), b, begin(C));
	assert( end(C) == e );
	return std::forward<C2D>(C);
}

//template<class Op, class A2D, class B2D, class C2D>
//C2D&& gemm(Op TA, Op TB, A2D const& A, B2D const& B, C2D&& C){
//	return gemm(TA, TB, 1., A, B, 0., std::forward<C2D>(C));
//}

template<class AA, class BB, class A2D, class B2D, class C2D>
C2D&& gemm(AA a, A2D const& A, B2D const& B, BB b, C2D&& C){
	switch(stride(C)){
		case  1: gemm(a, rotated(B), rotated(A), b, rotated(C)); break;
		default: switch(stride(A)){
			case  1: switch(stride(B)){
				case  1: gemm('T', 'T', a, rotated(B), rotated(A), b, C); break;
				default: gemm('N', 'T', a, B, rotated(A), b,          C);
			}; break;
			default: switch(stride(B)){
				case  1: gemm('T', 'N', a, rotated(B), A, b, C); break;
				default: gemm('N', 'N', a, B         , A, b, C);
			}
		}
	}
	return std::forward<C2D>(C);
#if 0
	switch(stride(A)){
		case 1: switch(stride(B)){
			case  1: switch(stride(C)){
				case  1: gemm('N', 'N', a, rotated(B), rotated(A), b, rotated(C)); break; // C = A.B, C^T = (A.B)^T, C^T = B^T.A^T
				default: gemm('T', 'T', a, rotated(B), rotated(A), b,         C );
			} break;
			default: switch(stride(C)){
				case  1: gemm('T', 'T', a, A, rotated(A), b, rotated(C)); break;
				default: gemm('N', 'T', a, B, rotated(A), b,         C );
			}
		}; break;
		default: switch(stride(B)){
			case 1: switch(stride(C)){
				case 1:  gemm(a, rotated(B), rotated(A), b, rotated(C)); break;
				default: gemm('T', 'N', a, rotated(B), A, b,         C );
			}; break;
			default: switch(stride(C)){
				case 1:  gemm(a, rotated(B), rotated(A), b, rotated(C) ); break; // C^T = (A*B)^T
				default: gemm('N', 'N', a, B, A, b,         C ); // C = A.B
			}
		}
	}
	return std::forward<C2D>(C);
#endif
}


}}}

#if _TEST_MULTI_ADAPTORS_BLAS_GEMV

#include "../../array.hpp"
#include "../../utility.hpp"

#include<complex>
#include<cassert>
#include<iostream>
#include<numeric>
#include<algorithm>

using std::cout;
namespace multi = boost::multi;

template<class M> void print(M const& C){
	for(int i = 0; i != size(C); ++i){
		for(int j = 0; j != size(C[i]); ++j)
			std::cout << C[i][j] << ' ';
		std::cout << std::endl;
	}
	std::cout << std::endl;
}

int main(){
	{
		multi::array<double, 2> const A = {
			{ 1., 1., 1.},
			{ 0., 0., 0.},
			{ 0., 0., 0.}
		};
		multi::array<double, 2> const B = {
			{ 1., 0., 0.},
			{ 1., 0., 0.},
			{ 1., 0., 0.},
		};
		using multi::blas::gemm;

		multi::array<double, 2> C1({3, 3});
		gemm('N', 'N', 1., B, A, 0., C1); // C = A*B , C^T = (B^T).(A^T) , if A, B, C are c-ordering
		print(C1);

		multi::array<double, 2> C2({3, 3});
		gemm('N', 'N', 1., A, B, 0., C2); // C = B*A , C^T = (A^T).(B^T) , if A, B, C are c-ordering
		print(C2);
	}
	{
		multi::array<double, 2> const A = {
			{ 1., 1., 1.},
			{ 0., 0., 0.},
			{ 0., 0., 0.},
			{ 0., 0., 0.}
		};
		multi::array<double, 2> const B = {
			{ 1., 0., 0., 0., 0.},
			{ 1., 0., 0., 0., 0.},
			{ 1., 0., 0., 0., 0,}
		};
		multi::array<double, 2> C({4, 5});

		using multi::blas::gemm;
		gemm('N', 'N', 1., B, A, 0., C); // C = A*B , C^T = (B^T).(A^T) , if A, B, C are c-ordering
		print(C);
	}
	{
		multi::array<double, 2> const A = {
			{ 0., 0., 0.},
			{ 1., 1., 1.}
		};
		multi::array<double, 2> const B = {
			{ 0., 1., 0.},
			{ 0., 1., 0.},
			{ 0., 1., 0.}
		};
		multi::array<double, 2> C({2, 3});

		using multi::blas::gemm;
		gemm('N', 'N', 1., B, A, 0., C); // C = A*B , C^T = (B^T).(A^T) , if A, B, C are c-ordering
		print(C);
	}
	{
		multi::array<double, 2> const A = {
			{ 0., 0., 0.},
			{ 0., 0., 1.},
			{ 0., 0., 0.}
		};
		multi::array<double, 2> const B = {
			{ 0., 1., 0.},
			{ 0., 0., 0.},
			{ 0., 0., 0.}
		};
		multi::array<double, 2> C({3, 3});

		using multi::blas::gemm; 
		gemm('T', 'T', 1., B, A, 0., C); // C = (B*A)^T , C^T = A*B, C=A^T*B^T , if A, B, C are c-ordering
		print(C);

	}
	{
		multi::array<double, 2> const A = {
			{ 0., 0., 0., 0.},
			{ 0., 1., 0., 0.},
			{ 0., 0., 0., 0.},
			{ 0., 0., 0., 0.},
			{ 0., 0., 0., 0.}
		};
		multi::array<double, 2> const B = {
			{ 0., 0., 0., 0., 0.},
			{ 0., 1., 0., 0., 0.},
			{ 0., 0., 0., 0., 0.}
		};
		multi::array<double, 2> C({4, 3});

		using multi::blas::gemm;
		gemm('T', 'T', 1., B, A, 0., C); //C = (B*A)^T, C^T = A*B, C=A^T*B^T, if A, B, C are c-ordering
		print(C);
	}
	{
		multi::array<double, 2> const A = {
			{ 0., 0., 0.},
			{ 0., 1., 0.},
			{ 0., 0., 0.},
			{ 0., 0., 0.}
		};
		multi::array<double, 2> const B = {
			{ 0., 0., 0., 0.},
			{ 0., 1., 0., 0.},
		};
		multi::array<double, 2> C({3, 2});

		using multi::blas::gemm;
		gemm('T', 'T', 1., B, A, 0., C); //C = (B*A)^T, C^T = A*B, C=A^T*B^T, if A, B, C are c-ordering
		print(C);
	}
	{
		multi::array<double, 2> const A = {
			{ 0., 1.},
			{ 0., 1.}
		};
		multi::array<double, 2> const B = {
			{ 1., 1.},
			{ 0., 0.}
		};
		multi::array<double, 2> C({2, 2});
		using multi::blas::gemm; 
		gemm('T', 'N', 1., B, A, 0., C); // C = A*(B^T) , C^T = B*(A^T) , if A, B, C are c-ordering
		print(C);
	}
	{
		multi::array<double, 2> const A = {
			{ 0., 1., 0., 0.},
			{ 0., 1., 0., 0.}
		};
		multi::array<double, 2> const B = {
			{ 1., 1., 1., 1.},
			{ 0., 0., 0., 0.},
			{ 0., 0., 0., 0.}
		};
		multi::array<double, 2> C({2, 3});
		using multi::blas::gemm;
		gemm('T', 'N', 1., B, A, 0., C); // C = A*(B^T) , C^T = B*(A^T) , if A, B, C are c-ordering
		print(C);
	}
	{
		multi::array<double, 2> const A = {
			{ 0., 1.},
			{ 0., 1.}
		};
		multi::array<double, 2> const B = {
			{ 1., 1.},
			{ 0., 0.}
		};
		multi::array<double, 2> C({2, 2});
		using multi::blas::gemm;
		gemm('N', 'T', 1., B, A, 0., C); // C = ((B^T)*A)^T , C^T = B^T*A, C = (A^T)*B , if A, B, C are c-ordering
		print(C);
	}
	{
		multi::array<double, 2> const A = {
			{ 0., 1., 0., 0.},
			{ 0., 1., 0., 0.},
		};
		multi::array<double, 2> const B = {
			{ 1., 1., 1.},
			{ 0., 0., 0.}
		};
		multi::array<double, 2> C({4, 3});
		using multi::blas::gemm;
		gemm('N', 'T', 1., B, A, 0., C); // C = ((B^T)*A)^T , C^T = B^T*A, C = (A^T)*B , if A, B, C are c-ordering
		print(C);
	}
	{
		multi::array<double, 2> const A = {
			{ 1., 2.},
			{ 3., 4.}
		};
		multi::array<double, 2> const B = {
			{ 5., 6.},
			{ 7., 8.}
		};
		multi::array<double, 2> C({2, 2});
		using multi::blas::gemm;
		gemm(1., A, B, 0., C); // C = A.B
		print(C);
	}
	{
		multi::array<double, 2> const A = {
			{ 1., 2.},
			{ 3., 4.}
		};
		multi::array<double, 2> const B = {
			{ 5., 6.},
			{ 7., 8.}
		};
		multi::array<double, 2> C({2, 2});
		using multi::blas::gemm;
		gemm(1., rotated(A), B, 0., C); // C = A^T.B
		print(C);
	}
	{
		multi::array<double, 2> const A = {
			{ 1., 2.},
			{ 3., 4.}
		};
		multi::array<double, 2> const B = {
			{ 5., 6.},
			{ 7., 8.}
		};
		multi::array<double, 2> C({2, 2});
		using multi::blas::gemm;
		gemm(1., A, rotated(B), 0., C); // C = A.B^T
		print(C);
	}
	{
		multi::array<double, 2> const A = {
			{ 1., 2.},
			{ 3., 4.}
		};
		multi::array<double, 2> const B = {
			{ 5., 6.},
			{ 7., 8.}
		};
		multi::array<double, 2> C({2, 2});
		using multi::blas::gemm;
		gemm(1., rotated(A), rotated(B), 0., C); // C = A^T.B^T
		print(C);
	}
	{
		multi::array<double, 2> const A = {
			{ 1., 2.},
			{ 3., 4.}
		};
		multi::array<double, 2> const B = {
			{ 5., 6.},
			{ 7., 8.}
		};
		multi::array<double, 2> C({2, 2});
		using multi::blas::gemm;
		gemm(1., A, B, 0., rotated(C)); // C^T = A.B, C = B^T.A^T
		print(C);
	}
	{
		multi::array<double, 2> const A = {
			{ 1., 2.},
			{ 3., 4.}
		};
		multi::array<double, 2> const B = {
			{ 5., 6.},
			{ 7., 8.}
		};
		multi::array<double, 2> C({2, 2});
		using multi::blas::gemm;
		gemm(1., rotated(A), B, 0., rotated(C)); // C^T = A^T.B, C = B^T.A
		print(C);
	}
	{
		multi::array<double, 2> const A = {
			{ 1., 2.},
			{ 3., 4.}
		};
		multi::array<double, 2> const B = {
			{ 5., 6.},
			{ 7., 8.}
		};
		multi::array<double, 2> C({2, 2});
		using multi::blas::gemm;
		gemm(1., A, rotated(B), 0., rotated(C)); // C^T = A.B^T, C = B.A^T
		print(C);
	}
	{
		multi::array<double, 2> const A = {
			{ 1., 2.},
			{ 3., 4.}
		};
		multi::array<double, 2> const B = {
			{ 5., 6.},
			{ 7., 8.}
		};
		multi::array<double, 2> C({2, 2});
		using multi::blas::gemm;
		gemm(1., rotated(A), rotated(B), 0., rotated(C)); // C^T = A^T.B^T, C = B.A
		print(C);
	}
}

#endif
#endif


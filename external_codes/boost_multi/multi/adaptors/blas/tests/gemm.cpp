#ifdef COMPILATION_INSTRUCTIONS
c++ -std=c++14 -Wall -Wextra -Wpedantic -lblas -DADD_ $0 -o $0x.x  && time $0x.x $@ && rm -f $0x.x; exit
#endif

#include "../../blas.hpp"

#include "../../../array.hpp"

#include<complex>
#include<cassert>

using std::cout;
namespace multi = boost::multi;

int main(){
	{
		multi::array<double, 2> const A = {
			{1.,  2.,  3.},
			{5.,  6.,  7.},
			{9., 10., 11.}
		};
		multi::array<double, 2> const B = {
			{10.,  12.,  31.},
			{51.,  61.,  74.},
			{91., 101., 111.}
		};
		multi::array<double, 2> C({size(A), std::get<1>(sizes(B))});
		using multi::blas::gemm;

// Miguel proposal:
//		assert(A.stride(1) == 1 and B.stride(1) == 1 and C.stride(1) == 1)
		gemm('N', 'N', A         , B         , C         ); // C = B*A
//		gemm('N', 'N', rotated(A), rotated(B), rotated(C)); // fail [T(C) = T(A)*T(B),  C = B*A]
//		gemm('T', 'T', A         , B         , C         ); // C = T(A)*T(B)   T(C) = B*A
//		gemm('N', 'T', ...)                              // C = A*T(B)
//		gemm('T', 'N', ...)                              // C = T(A)*B
//		gemm('N', 'H'  ...)                              // C = A*H(B)
//		A2 = rotated(A);

//		gemm('N', 'N', A2       , B        , C        );

//		gemm('T', 'T', rotate(A), rotate(B), C        ); // C = B*A = T(T(A)*T(B)) or T(C) = T(A)*T(B)
//		gemm('T', 'T', A,         B        , rotate(C)); // C = A*B = T(T(B)*T(A))



//
		gemm('N', 'N', A        , B        , C        ); // C = B*A = T(T(A)*T(B)) or T(C) = T(A)*T(B)
/*		gemm('N', 'N', rotate(A), rotate(B), rotate(C)); // C = A*B

		gemm('T', 'T', rotate(A), rotate(B), C        ); // C = B*A = T(T(A)*T(B)) or T(C) = T(A)*T(B)
		gemm('T', 'T', A,         B        , rotate(C)); // C = A*B = T(T(B)*T(A))
*/
//        auto C = gemm(          A, B);
		
		
		
//		gemm('T', 'T', A, B, C); // 
//		for(auto i = 0; i != 3; ++i){
//			for(auto j = 0; j != 3; ++j)
//				cout << C[i][j] <<' ';
//			cout << std::endl;
//		}
	}
}


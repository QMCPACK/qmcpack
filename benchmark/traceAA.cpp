#if COMPILATION_INSTRUCTIONS
c++ -std=c++17 -O2 -Ofast -DNDEBUG -Wall -Wextra -Wpedantic $0 -o $0x.x && $0x.x 16384 128 128 $@ && rm -f $0x.x; exit
#endif

#include<iostream>
#include<cstdlib>
#include<sys/time.h>
#include<numeric>

#include "../include/multi/array.hpp"

using std::cout;
using std::cerr;
using std::endl;

double getTime(){
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return double(tv.tv_sec)+double(tv.tv_usec)/1000000.0;
}

struct timer{
	double ct_;
	timer(){ct_ = getTime();}
	template<class T> 
	auto
	operator()(T&& t) const{
		return std::make_pair(getTime()-ct_, std::forward<T>(t));
	}
};

double trace_AA_miguel(double const* A, int const N, int const ni, int const nj)
{
  double s(0.0);
  int nbi = (N + ni - 1) / ni;
  int nbj = (N + nj - 1) / nj;
  double* ss = new double[N];
  for(int bi=0; bi<nbi; ++bi) {
    int i0 = bi*ni;
    int iN = std::min((bi+1)*ni,N);
    for(int bj=0; bj<nbj; ++bj) {
      int j0 = bj*nj;
      int jN = std::min((bj+1)*nj,N);
	  std::fill_n(ss, N, 0.);
      for(int i=i0; i<iN; ++i) {
        for(int j=j0; j<jN; ++j)
          ss[i] += A[i*N+j]*A[j*N+i];
      }
      for(int i = i0; i < iN; ++i) s += ss[i];
    }
  }
  delete[] ss;
  return s;
}

double trace_AA(double const* A, std::ptrdiff_t const N, std::ptrdiff_t const ni, std::ptrdiff_t const nj){
	double s = 0.;
	std::ptrdiff_t const nbi = (N + ni - 1) / ni;
	std::ptrdiff_t const nbj = (N + nj - 1) / nj;
	for(std::ptrdiff_t bi = 0; bi != nbi; ++bi){
	    std::ptrdiff_t i0 = bi*ni;
	    std::ptrdiff_t iN = std::min((bi+1)*ni,N);
	    for(std::ptrdiff_t bj = 0; bj != nbj; ++bj) {
			std::ptrdiff_t j0 = bj*nj;
			std::ptrdiff_t jN = std::min((bj+1)*nj,N);
			for(std::ptrdiff_t i = i0; i != iN; ++i){
				for(std::ptrdiff_t j = j0; j != jN; ++j){
					s += A[i*N+j]*A[j*N+i];
				}
     		}
	    }
	}
	return s;
}

template<class Matrix, class index = typename Matrix::index>
double trace_AA(Matrix const& A, index const ni, index const nj){
	index const N = A.size();
	double s = 0.;
	index const nbi = (N + ni - 1) / ni;
	index const nbj = (N + nj - 1) / nj;
//	auto const& At = A.rotated();
	for(index bi = 0; bi != nbi; ++bi){
		index const i0 = bi*ni;
		index const iN = std::min((bi+1)*ni, N);
		for(index bj=0; bj != nbj; ++bj) {
			index const j0 = bj*nj;
			index const jN = std::min((bj+1)*nj, N);
			for(index i=i0; i != iN; ++i){
			//	auto const& Ai = A[i];
			//	auto const& Ati = At[i];
				for(int j=j0; j != jN; ++j){
					s += A[i][j]*A[j][i]; //Ai[j]*Ati[j]; // s += A[i*N+j]*A[j*N+i];
				}
			}
		}
	}
	return s;
}

template<class Matrix, class index = typename Matrix::index>
double trace_AA2(Matrix const& A, index const ni, index const nj){
	index const N = A.size();
	double s = 0.;
	index const nbi = (N + ni - 1) / ni;
	index const nbj = (N + nj - 1) / nj;
	for(index bi = 0; bi != nbi; ++bi){
		index const i0 = bi*ni;
		index const iN = std::min((bi+1)*ni, N);
		auto const& rowt = A.sliced(i0, iN).rotated();
		auto const& col = A.rotated().sliced(i0, iN).unrotated();
		for(index bj = 0; bj != nbj; ++bj) {
			index const j0 = bj*nj;
			index const jN = std::min((bj+1)*nj, N);
			auto const& B1t = rowt.sliced(j0, jN);
			auto const& B2  = col.sliced(j0, jN);
			index const B1ts0 = B1t.size();
			index const B1ts1 = B2.size();
			for(index j = 0; j != B1ts0; ++j){
				auto const& B1tj = B1t[j];
				auto const& B2j = B2[j];
				for(index i = 0; i != B1ts1; ++i){
					s += B1tj[i]*B2j[i];
				}
			}
		}
	}
	return s;
}

template<class Matrix1, class Matrix2, class index = typename Matrix1::index>
double trace_AA4(Matrix1 const& A, Matrix2 const& B, index const ni, index const nj){
	index const AN = A.size();
	index const BM = B.size();
	double s = 0.;
	index const nbi = (AN + ni - 1) / ni;
	index const nbj = (BM + nj - 1) / nj;
	for(index bi = 0; bi != nbi; ++bi){
		index const i0 = bi*ni;
		index const iN = std::min((bi+1)*ni, AN);
		auto const& Arowt = A.sliced(i0, iN).rotated();
		auto const& Bcol = B.sliced(i0, iN).unrotated();
		for(index bj = 0; bj != nbj; ++bj) {
			index const j0 = bj*nj;
			index const jN = std::min((bj+1)*nj, BM);
			auto const& B1t = Arowt.sliced(j0, jN);
			auto const& B2  = Bcol.sliced(j0, jN);
			if(ni >= 32 and nj >= 32){
				s += trace_AA4(B1t, B2, ni/4, nj/4);
			}else
			{
				index const B1ts0 = B1t.size();
				index const B1ts1 = B2.size();
				for(index j = 0; j != B1ts0; ++j){
					auto const& B1tj = B1t[j];
					auto const& B2j = B2[j];
					for(index i = 0; i != B1ts1; ++i){
						s += B1tj[i]*B2j[i];
					}
				}
			}
		}
	}
	return s;
}

int main(int, char*[]){
	cout << sizeof(double);
//	if(argc < 4) {
//		cerr<<" test N ni nj\n";
//		exit(1);
//	}
	int const N = 32768;// atoi(argv[1]);
	cout << "N = " << N << '\n';
	int const ni = 64; // atoi(argv[2]);
	int const nj = 64; //atoi(argv[3]);
	{
		cout << "Miguel traceAA\n";
		double* A = new double[N*N];
		std::iota(A, A + N*N, 1.11);
		double ss = trace_AA_miguel(A, N, ni, nj);
		std::cout << ss << std::endl;
		std::iota(A, A + N*N, 1.23);
		double t1 = getTime();
		double s = trace_AA_miguel(A, N, ni, nj);
		t1 = getTime()-t1;
		cout<<" S: " << s <<endl;	
		cout<<" Time: " << t1 <<endl;
		delete [] A;
	}
	{
		cout << "Miguel with ptrdiff_t traceAA\n";
		double* A = new double[N*N];
		std::iota(A, A + N*N, 1.11);
		double ss = trace_AA(A, N, ni, nj);
		std::cout << ss <<'\n';
		std::iota(A, A + N*N, 1.23);
		double t1 = getTime();
		double s = trace_AA(A, N, ni, nj);
		t1 = getTime()-t1;
		cout<<" S: " <<s <<endl;	
		cout<<" Time: " << t1 <<endl;
		delete [] A;
	}
	{
		cout << "Multi manual blocking traceAA\n";
		boost::multi::array<double, 2> A({N, N});
		std::iota(A.elements().begin(), A.elements().end(), 1.11);
		double ss = trace_AA(A, ni, nj);
		std::cout << ss <<'\n';
		std::iota(A.elements().begin(), A.elements().end(), 1.23);
		double t1 = getTime();
		double s = trace_AA(A, ni, nj);
		t1 = getTime()-t1;
		cout<<" S: " <<s <<endl;
		cout<<" Time: " <<t1 <<endl;
	}
	{
		cout << "Multi subarray blocking traceAA\n";
		boost::multi::array<double, 2> A({N, N});
		std::iota(A.elements().begin(), A.elements().end(), 1.11);
		double ss = trace_AA2(A, ni, nj);
		std::cout << ss <<'\n';
		std::iota(A.elements().begin(), A.elements().end(), 1.23);
		double t1 = getTime();
		double s = trace_AA2(A, ni, nj);
		t1 = getTime()-t1;
		cout<<" S: " << s <<endl;
		cout<<" Time: " << t1 <<endl;
	}
	{
		cout << "Multi RECURSIVE subarray blocking traceAA\n";
		boost::multi::array<double, 2> A({N, N});
		std::iota(A.elements().begin(), A.elements().end(), 1.11);
		double ss = trace_AA4(A, A.rotated(), A.size()/4, A.size()/4);
		std::cout << ss <<'\n';

		std::iota(A.elements().begin(), A.elements().end(), 1.23);
		double t1 = getTime();
		double s = trace_AA4(A, A.rotated(), A.size()/4, A.size()/4);
		t1 = getTime()-t1;
		cout<<" S: " <<s <<endl;
		cout<<" Time: "<< t1 <<endl;
	}

	return 0;
}

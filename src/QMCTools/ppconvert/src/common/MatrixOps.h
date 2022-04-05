//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Paul R. C. Kent, kentpr@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Alfredo A. Correa, correaa@llnl.gov, Lawrence Livermore National Laboratory
//
// File created by: Paul R. C. Kent, kentpr@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


//           http://pathintegrals.info                     //
/////////////////////////////////////////////////////////////

#ifndef MATRIX_OPS_H
#define MATRIX_OPS_H
#include "Blitz.h"

#include<iostream>  // for std::cerr

void LinearLeastSquares(Array<double, 2>& A, Array<double, 1>& x, Array<double, 1>& b);
void MatVecProd(Array<double, 2>& A, Array<double, 1>& x, Array<double, 1>& Ax);

//void LUdecomp(Array<double, 2>& A, Array<int, 1>& perm, double& sign);

//void LUsolve(Array<double, 2>& LU, Array<int, 1>& perm, Array<double, 1>& b);


struct lup{ // LU method for decomposition and solution

// translated from  https://en.wikipedia.org/wiki/LU_decomposition#C_code_example
template<class Matrix, class Permutation>
static auto decompose(Matrix&& A, Permutation&& P, double tol = std::numeric_limits<double>::epsilon()){
	throw;
	std::iota(begin(P), end(P), typename std::decay_t<Permutation>::value_type{0});
	auto const N = std::min(size(A), size(~A));
	assert( P.size() >= N );

	auto&& ret = A({0, N}, {0, N});
	for(auto i : extension(ret)){
		if(lup::permute_max_diagonal(A, P, i) < tol) return A({0, i}, {0, i});

		for(auto&& row : A({i + 1, N})){
			auto&& urow = row({i + 1, N});
			std::transform(
				cbegin(urow), cend(urow), cbegin(A[i]({i + 1, N})), begin(urow), 
				[f = row[i] /= A[i][i]](auto const& a, auto const& b){return a - f*b;}
			);
		}
	}
	return std::move(ret);
}

template<class Matrix, class Permutation, class VectorSol>
static auto solve(Matrix const& LU, Permutation const& P, VectorSol&& x) -> VectorSol&&{
	return upper_solve(LU, lower_solve(LU, permute(P, x)));
}

private:

template<class Matrix, class Permutation, class Index>
static auto permute_max_diagonal(Matrix&& LU, Permutation&& P, Index i){
	auto mi = std::max_element(begin(LU) + i, end(LU), [i](auto const& a, auto const& b){return std::abs(a[i]) < std::abs(b[i]);}) - begin(LU);
	     swap(LU[i], LU[mi]); 
	std::swap(P [i], P [mi]);
	return std::abs(LU[i][i]);
}

template<class Permutation, class Vector>
static auto permute(Permutation const& p, Vector&& data) -> Vector&&{
	assert(size(p) <= size(data));
	using index = typename Permutation::size_type;
	for(index i = 0; i != size(p); ++i){
		index k = p[i];
		for( ; k > i; k = p[k]){}
		index pk = p[k];
		if(k >=i and pk != i){
			auto const t = data[i];
			for( ; pk != i; k = pk, pk = p[k]){
				data[k] = data[pk];
			};
			data[k] = t;
		}
	}
	return std::forward<Vector>(data);
}

template<class LUMatrix, class Vector>
static auto lower_solve(LUMatrix const& LU, Vector&& x) -> Vector&&{
	assert(size(LU) <= size(x));
	auto const N = size(LU);
	for(typename LUMatrix::size_type i = 0; i != N; ++i){
		auto const& Lrowi = LU[i]({0, i});
		x[i] -= std::inner_product(begin(Lrowi), end(Lrowi), cbegin(x), 0.);
	}
	return std::forward<Vector>(x);
}

template<class LUMatrix, class Vector>
static auto upper_solve(LUMatrix const& LU, Vector&& x) -> Vector&&{
	assert(size(LU) <= size(x));
	auto const N = size(LU);
	for(typename LUMatrix::size_type i = N - 1; i >= 0; --i){
		auto const& Urowi = LU[i]({i + 1, N});
		(x[i] -= std::inner_product(begin(Urowi), end(Urowi), cbegin(x) + i + 1, 0.)) /= LU[i][i];
	}
	return std::forward<Vector>(x);
}

};

void SVdecomp(Array<double, 2>& A, Array<double, 2>& U, Array<double, 1>& S, Array<double, 2>& V);

void SVdecomp(Array<std::complex<double>, 2>& A,
              Array<std::complex<double>, 2>& U,
              Array<double, 1>& S,
              Array<std::complex<double>, 2>& V);


void SymmEigenPairs(const Array<double, 2>& A, int NumPairs, Array<double, 1>& Vals, Array<double, 2>& Vectors);

void SymmEigenPairs(const Array<std::complex<double>, 2>& A,
                    int NumPairs,
                    Array<double, 1>& Vals,
                    Array<std::complex<double>, 2>& Vectors);

// Orthogonalizes matrix A with the polar decomposition.  This returns
// the matrix closest to A which is orthogonal.  A must be square.
void PolarOrthogonalize(Array<std::complex<double>, 2>& A);

const Array<double, 2> operator*(const Array<double, 2>& A, const Array<double, 2>& B);

double InnerProduct(const Array<double, 1>& A, const Array<double, 1>& B);

void OuterProduct(const Array<double, 1>& A, const Array<double, 1>& B, Array<double, 2>& AB);


void MatMult(const Array<double, 2>& A, const Array<double, 2>& B, Array<double, 2>& C);
void MatMult(const Array<std::complex<double>, 2>& A,
             const Array<std::complex<double>, 2>& B,
             Array<std::complex<double>, 2>& C);

double Determinant(const Array<double, 2>& A);
std::complex<double> Determinant(const Array<std::complex<double>, 2>& A);

/// This function returns the determinant of A and replaces A with its
/// cofactors.
double DetCofactors(Array<double, 2>& A, Array<double, 1>& work);

/// This function returns the worksize needed by the previous function.
int DetCofactorsWorksize(int N);


/// Complex versions of the functions above
std::complex<double> ComplexDetCofactors(Array<std::complex<double>, 2>& A, Array<std::complex<double>, 1>& work);
int ComplexDetCofactorsWorksize(int N);

double GJInverse(Array<double, 2>& A);
double GJInversePartial(Array<double, 2>& A);

Array<double, 2> Inverse(Array<double, 2>& A);

inline void OutOfPlaceTranspose(Array<double, 2>& A)
{
  int m = A.rows();
  int n = A.cols();
  Array<double, 2> Atrans(n, m);
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      Atrans(j, i) = A(i, j);
  A.resize(n, m);
  A = Atrans;
}

inline void OutOfPlaceTranspose(Array<std::complex<double>, 2>& A)
{
  int m = A.rows();
  int n = A.cols();
  Array<std::complex<double>, 2> Atrans(n, m);
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      Atrans(j, i) = A(i, j);
  A.resize(n, m);
  A = Atrans;
}

inline void Transpose(Array<double, 2>& A)
{
  int m = A.rows();
  int n = A.cols();
  if (m != n)
    OutOfPlaceTranspose(A);
  else
  {
    for (int i = 0; i < m; i++)
      for (int j = i + 1; j < m; j++)
        std::swap(A(i, j), A(j, i));
  }
}

inline void Transpose(Array<std::complex<double>, 2>& A)
{
  int m = A.rows();
  int n = A.cols();
  if (m != n)
    OutOfPlaceTranspose(A);
  else
  {
    for (int i = 0; i < m; i++)
      for (int j = i + 1; j < m; j++)
        std::swap(A(i, j), A(j, i));
  }
}


inline void Print(const Mat3& A)
{
  fprintf(stderr, "%14.6e %14.6e %14.6e\n", A(0, 0), A(0, 1), A(0, 2));
  fprintf(stderr, "%14.6e %14.6e %14.6e\n", A(1, 0), A(1, 1), A(1, 2));
  fprintf(stderr, "%14.6e %14.6e %14.6e\n", A(2, 0), A(2, 1), A(2, 2));
}

inline void CubicFormula(double a, double b, double c, double d, double& x1, double& x2, double& x3)
{
  double A = b / a;
  double B = c / a;
  double C = d / a;
  double Q = (A * A - 3.0 * B) / 9.0;
  double R = (2.0 * A * A * A - 9.0 * A * B + 27.0 * C) / 54.0;
  //cerr << "Q = " << Q << " R = " << R << "\n";
  if ((R * R) < (Q * Q * Q))
  {
    double theta    = acos(R / sqrt(Q * Q * Q));
    double twosqrtQ = 2.0 * sqrt(Q);
    double third    = 1.0 / 3.0;
    double thirdA   = third * A;
    x1              = -twosqrtQ * cos(third * theta) - thirdA;
    x2              = -twosqrtQ * cos(third * (theta + 2.0 * M_PI)) - thirdA;
    x3              = -twosqrtQ * cos(third * (theta - 2.0 * M_PI)) - thirdA;
  }
  else
  {
    std::cerr << "Complex roots detected in CubicFormula.\n";
    exit(1);
  }
}


inline void TestCubicFormula()
{
  double x1 = 2.3;
  double x2 = 1.1;
  double x3 = 9.2;
  double a  = 1.0;
  double b  = -(x1 + x2 + x3);
  double c  = (x1 * x2 + x1 * x3 + x2 * x3);
  double d  = -x1 * x2 * x3;

  double y1, y2, y3;
  CubicFormula(a, b, c, d, y1, y2, y3);
  std::cerr << "y1 = " << y1 << "\n";
  std::cerr << "y2 = " << y2 << "\n";
  std::cerr << "y3 = " << y3 << "\n";
}


inline void EigVals(const Mat3& M, double& lambda1, double& lambda2, double& lambda3)
{
  double a = -1.0;
  double b = M(0, 0) + M(1, 1) + M(2, 2);
  double c = (-M(0, 0) * M(1, 1) - M(0, 0) * M(2, 2) - M(1, 1) * M(2, 2) + M(0, 1) * M(1, 0) + M(0, 2) * M(2, 0) +
              M(1, 2) * M(2, 1));
  double d = (M(0, 0) * M(1, 1) * M(2, 2) + M(0, 1) * M(1, 2) * M(2, 0) + M(0, 2) * M(1, 0) * M(2, 1) -
              M(0, 1) * M(1, 0) * M(2, 2) - M(0, 2) * M(1, 1) * M(2, 0) - M(0, 0) * M(1, 2) * M(2, 1));
  CubicFormula(a, b, c, d, lambda1, lambda2, lambda3);
}


inline void Eig(const Mat3& M, Mat3& U, Vec3& Lambda)
{
  EigVals(M, Lambda(0), Lambda(1), Lambda(2));
  for (int i = 0; i < 3; i++)
  {
    double L    = Lambda(i);
    U(0, i)     = 1.0;
    U(1, i)     = (M(1, 2) / M(0, 2) * (M(0, 0) - L) - M(1, 0)) / (M(1, 1) - M(1, 2) / M(0, 2) * M(0, 1) - L);
    U(2, i)     = ((L - M(1, 1)) * U(1, i) - M(1, 0)) / M(1, 2);
    double norm = 1.0 / sqrt(U(0, i) * U(0, i) + U(1, i) * U(1, i) + U(2, i) * U(2, i));
    U(0, i) *= norm;
    U(1, i) *= norm;
    U(2, i) *= norm;
  }
}


inline double det(const Mat2& C) { return (C(0, 0) * C(1, 1) - C(0, 1) * C(1, 0)); }

// inline double det (const Mat3 &C)
// {
//   double d0, d1, d2;
//   d0 = C(1,1)*C(2,2) - C(1,2)*C(2,1);
//   d1 = C(1,0)*C(2,2) - C(2,0)*C(1,2);
//   d2 = C(1,0)*C(2,1) - C(2,0)*C(1,1);
//   return (C(0,0)*d0 - C(0,1)*d1 + C(0,2)*d2);
// }

// inline Mat3 Inverse (const Mat3 &A)
// {
//   Mat3 CoFacts, Inv;
//   CoFacts(0,0) = A(1,1)*A(2,2) - A(1,2)*A(2,1);
//   CoFacts(0,1) = -(A(1,0)*A(2,2) - A(2,0)*A(1,2));
//   CoFacts(0,2) = A(1,0)*A(2,1) - A(2,0)*A(1,1);
//   CoFacts(1,0) = -(A(0,1)*A(2,2) - A(0,2)*A(2,1));
//   CoFacts(1,1) = A(0,0)*A(2,2) - A(0,2)*A(2,0);
//   CoFacts(1,2) = -(A(0,0)*A(2,1)-A(0,1)*A(2,0));
//   CoFacts(2,0) = A(0,1)*A(1,2) - A(0,2)*A(1,1);
//   CoFacts(2,1) = -(A(0,0)*A(1,2) - A(1,0)*A(0,2));
//   CoFacts(2,2) = A(0,0)*A(1,1) - A(0,1)*A(1,0);

//   double det = (A(0,0) * CoFacts(0,0) +
// 		A(0,1) * CoFacts(0,1) +
// 		A(0,2) * CoFacts(0,2));
//   double detinv = 1.0/det;
//   Inv(0,0)=CoFacts(0,0)*detinv;
//   Inv(0,1)=CoFacts(1,0)*detinv;
//   Inv(0,2)=CoFacts(2,0)*detinv;

//   Inv(1,0)=CoFacts(0,1)*detinv;
//   Inv(1,1)=CoFacts(1,1)*detinv;
//   Inv(1,2)=CoFacts(2,1)*detinv;

//   Inv(2,0)=CoFacts(0,2)*detinv;
//   Inv(2,1)=CoFacts(1,2)*detinv;
//   Inv(2,2)=CoFacts(2,2)*detinv;
//   return (Inv);
// }


// inline Mat3 Transpose(const Mat3 &A)
// {
//   Mat3 B;
//   B(0,0) = A(0,0);
//   B(0,1) = A(1,0);
//   B(0,2) = A(2,0);
//   B(1,0) = A(0,1);
//   B(1,1) = A(1,1);
//   B(1,2) = A(2,1);
//   B(2,0) = A(0,2);
//   B(2,1) = A(1,2);
//   B(2,2) = A(2,2);
//   return (B);
// }

void CholeskyBig(Array<double, 2>& A);


inline Mat3 Cholesky(Mat3 C)
{
  Mat3 L;
  double L00Inv;
  L       = 0.0;
  L(0, 0) = sqrt(C(0, 0));
  L00Inv  = 1.0 / L(0, 0);
  L(1, 0) = C(1, 0) * L00Inv;
  L(2, 0) = C(2, 0) * L00Inv;
  L(1, 1) = sqrt(C(1, 1) - L(1, 0) * L(1, 0));
  L(2, 1) = (C(2, 1) - L(2, 0) * L(1, 0)) / L(1, 1);
  L(2, 2) = sqrt(C(2, 2) - (L(2, 0) * L(2, 0) + L(2, 1) * L(2, 1)));
  return (L);
}


inline void TestCholesky()
{
  Mat3 V, D, C, L;
  V(0, 0) = 1.0;
  V(0, 1) = 2.0;
  V(0, 2) = 3.0;
  V(1, 0) = 6.0;
  V(1, 1) = 1.0;
  V(1, 2) = 9.0;
  V(2, 0) = 7.0;
  V(2, 1) = 5.0;
  V(2, 2) = 7.0;

  D(0, 0) = 1.0;
  D(0, 1) = 0.0;
  D(0, 2) = 0.0;
  D(1, 0) = 0.0;
  D(1, 1) = 2.0;
  D(1, 2) = 0.0;
  D(2, 0) = 0.0;
  D(2, 1) = 0.0;
  D(2, 2) = 3.0;

  C = V * D * Transpose(V);
  std::cerr << "C = " << C << "\n";
  std::cerr << "V = " << V << "\n";
  L = Cholesky(C);
  std::cerr << "L = " << L << "\n";
  std::cerr << "L L^T = " << L * Transpose(L) << "\n";

  Mat3 E, U;
  Vec3 Lambda;
  for (int i = 0; i < 10000000; i++)
  {
    Eig(C, U, Lambda);
    //E = Inverse(C);
    //L = Cholesky (E);
  }
}

#endif

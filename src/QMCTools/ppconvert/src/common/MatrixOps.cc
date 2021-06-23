//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Paul R. C. Kent, kentpr@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Paul R. C. Kent, kentpr@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
    
    
//           http://pathintegrals.info                     //
/////////////////////////////////////////////////////////////

#include "MatrixOps.h"

#include "../config.h"

//#ifdef NOUNDERSCORE 
//#define FORT(name) name
//#else
#define FORT(name) name ## _
//#endif 

#define F77_DGELS  F77_FUNC(dgels,DGELS)
#define F77_ZGESVD F77_FUNC(zgesvd,ZGESVD)
#define F77_DGESVD F77_FUNC(dgesvd,DGESVD)
#define F77_DGETRF F77_FUNC(dgetrf,DGETRF)
#define F77_ZGETRF F77_FUNC(zgetrf,ZGETRF)
#define F77_DGETRI F77_FUNC(dgetri,DGETRI)
#define F77_ZGETRI F77_FUNC(zgetri,ZGETRI)
#define F77_DGEMM  F77_FUNC(dgemm,DGEMM)
#define F77_ZGEMM  F77_FUNC(zgemm,ZGEMM)
#define F77_DSYEVR F77_FUNC(dsyevr,DSYEVR)
#define F77_ZHEEVR F77_FUNC(zheevr,ZHEEVR)
#define F77_DGEMV  F77_FUNC(dgemv,DGEMV)
#define F77_DPOTRF  F77_FUNC(dpotrf,DPOTRF)


extern "C" void 
F77_DGELS (char *transa, int *m, int *n, int *nrhs, 
	   double *A, int *LDA, double *B,int *LDB,
	   double *work, int *LDWORK,int *INFO);

extern "C" void 
F77_DGESVD (char *JOBU, char* JOBVT, int *M, int *N,
	    double *A, int *LDA, double *S, double *U,
	    int *LDU, double *VT, int *LDVT, double *work,
	    int *LWORK, int *INFO);

extern "C" void 
F77_ZGESVD (char *JOBU, char* JOBVT, int *M, int *N,
	    std::complex<double> *A, int *LDA, double *S, std::complex<double> *U,
	    int *LDU, std::complex<double> *VT, int *LDVT, std::complex<double> *work,
	    int *LWORK, double *work2, int *INFO);


extern "C" void 
F77_DGETRF(int *m, int *n, double A[], int *lda, int ipiv[], int *info);

extern "C" void 
F77_DPOTRF(char *UPLO, int *n, double A[], int *lda, int *info);


extern "C" void 
F77_ZGETRF(int *m, int *n, std::complex<double> A[], 
	   int *lda, int ipiv[], int *info);

extern "C" void 
F77_DGETRI (int *N, double A[], int *lda, int ipiv[], double work[], 
	    int *lwork, int *info);

extern "C" void 
F77_ZGETRI (int *N, std::complex<double> A[], int *lda, int ipiv[], 
	    std::complex<double> work[], int *lwork, int *info);

extern "C" void 
F77_DGEMM (char *transA, char *transB, int *m, int *n, int *k,
	   double *alpha, const double *A, int *lda, const double *B, int *ldb,
	   double *beta,  double *C, int *ldc);

extern "C" void 
F77_ZGEMM (char *transA, char *transB, int *m, int *n, int *k,
	   std::complex<double> *alpha, const std::complex<double> *A, int *lda, const std::complex<double> *B, 
	   int *ldb, std::complex<double> *beta,  std::complex<double> *C, int *ldc);

extern "C" void 
F77_DSYEVR (char *JobType, char *Range, char *UpperLower, 
	    int *N, double *Amat, int *LDA,
	    double *VL, double *VU,
	    int *IL, int *IU, 
	    double *AbsTolerance, int *M,
	    double *EigVals, 
	    double *EigVecs, int *LDEigVecs, int *ISuppZ,
	    double *Work, int *Lwork, 
	    int *IWorkSpace, int *LIwork,
	    int *Info);

extern "C" void 
F77_ZHEEVR (char *JobType, char *Range, char *UpperLower, 
	    int *N, std::complex<double> *Amat, int *LDA,
	    double *VL, double *VU,
	    int *IL, int *IU, 
	    double *AbsTolerance, int *M,
	    double *EigVals, 
	    std::complex<double> *EigVecs, int *LDEigVecs, 
	    int *ISuppZ,
	    std::complex<double> *Work, int *Lwork, 
	    double *Rwork, int *LRwork,
	    int *IWorkSpace, int *LIwork,
	    int *Info);


const Array<double,2> operator*(const Array<double,2> &A,
				const Array<double,2> &B)
{
  int m = A.rows();
  int n = B.cols();
  int k = A.cols();
  assert (B.rows() == k);
  // We use "transpose" operation because we have C ordering, which fortran
  // thinks is transposed.
  char transA = 'T';
  char transB = 'T';
  double alpha = 1.0;
  double beta = 0.0;
//  blitz::GeneralArrayStorage<2> colMajor;
//  colMajor.ordering() = blitz::firstDim, blitz::secondDim;
//  Array<double,2> C(m,n,colMajor);
	Array<double, 2> C({m, n}, 0.);
  F77_DGEMM (&transA, &transB, &m, &n, &k, &alpha, A.data(), &k, 
           B.data(), &n, &beta, C.data(), &m);
  return C;
}

const Array<std::complex<double>,2> operator*(const Array<std::complex<double>,2> &A,
					 const Array<std::complex<double>,2> &B)
{
  int m = A.rows();
  int n = B.cols();
  int k = A.cols();
  assert (B.rows() == k);
  // We use "transpose" operation because we have C ordering, which fortran
  // thinks is transposed.
  char transA = 'T';
  char transB = 'T';
  std::complex<double> alpha(1.0, 0.0);
  std::complex<double> beta(0.0, 0.0);
//  blitz::GeneralArrayStorage<2> colMajor;
//  colMajor.ordering() = blitz::firstDim, blitz::secondDim;
//  Array<std::complex<double>,2> C(m,n,colMajor);
	Array<std::complex<double>, 2> C({m, n}, std::complex<double>{});
  F77_ZGEMM (&transA, &transB, &m, &n, &k, &alpha, A.data_elements(), &k, 
           B.data(), &n, &beta, C.data(), &m);
  return C;
}


//void MatMult (const Array<double,2> &A, const Array<double,2> &B,
//	      Array<double,2> &C)
//{
//  int m = A.rows();
//  int n = B.cols();
//  int k = A.cols();
//  assert (B.rows() == k);
//  // We use "transpose" operation because we have C ordering, which fortran
//  // thinks is transposed.
//  char transA = 'T';
//  char transB = 'T';
//  double alpha = 1.0;
//  double beta = 0.0;
//  blitz::GeneralArrayStorage<2> colMajor;
//  colMajor.ordering() = blitz::firstDim, blitz::secondDim;
//  F77_DGEMM (&transA, &transB, &m, &n, &k, &alpha, A.data(), &k, 
//	     B.data(), &n, &beta, C.data(), &m);
//}

//double 
//InnerProduct(const Array<double,1> &A,
//	     const Array<double,1> &B)
//{
//  assert(A.size()==B.size());
//  double total=0.0;
//  for (int counter=0;counter<A.size();counter++)
//    total+=A(counter)*B(counter);
//  return total;

//}


//void
//OuterProduct(const Array<double,1> &A,
//	     const Array<double,1> &B,
//	     Array<double,2> &AB)
//{
//  double total=0.0;
//  AB.resize(A.size(),B.size());
//  for (int i=0;i<A.size();i++)
//    for (int j=0;j<B.size();j++)
//      AB(i,j)=A(i)*B(j);
//}

//Note that A gets corrupted in this process and b gets returned
//as the answer
//void
//LinearLeastSquares(Array<double,2> &A, Array<double,1> &x,
//		   Array<double,1> &b)
//{
//  char trans;
//  trans='T';
//  ///These n and m are "transposed" because of Fortran ordering
//  int n=A.rows();
//  int m=A.cols();
//  ///

//  int ldb=b.size();
//  //  cerr<<"The value of ldb is "<<ldb<<endl;
//  Array<double,1> work(1);
//  int ldwork=-1;
//  int info=0;
//  int nrhs=1;
//  F77_DGELS(&trans,&m,&n,&nrhs,A.data(),&m,b.data(),&ldb,work.data(),
//	    &ldwork,&info);
//  work.resize((int)work(0));
//  ldwork=work.size();
//  F77_DGELS(&trans,&m,&n,&nrhs,A.data(),&m,b.data(),&ldb,work.data(),
//	    &ldwork,&info);

//}
//void MatMult (const Array<std::complex<double>,2> &A, const Array<std::complex<double>,2> &B,
//	      Array<std::complex<double>,2> &C)
//{
//  int m = A.rows();
//  int n = B.cols();
//  int k = A.cols();
//  assert (B.rows() == k);
//  // We use "transpose" operation because we have C ordering, which fortran
//  // thinks is transposed.
//  char transA = 'T';
//  char transB = 'T';
//  std::complex<double> alpha = 1.0;
//  std::complex<double> beta = 0.0;
//  blitz::GeneralArrayStorage<2> colMajor;
//  colMajor.ordering() = blitz::firstDim, blitz::secondDim;
//  F77_ZGEMM (&transA, &transB, &m, &n, &k, &alpha, A.data(), &k, 
//	     B.data(), &n, &beta, C.data(), &m);
//  Transpose(C);
//}


//double Determinant (const Array<double,2> &A)
//{
//  int m = A.rows();
//  int n = A.cols();
//  assert (m == n);  // Cannot take a determinant of a non-square
//		    // matrix
//  if (A.rows() == 1)
//    return (A(0,0));
//  if (A.rows() == 2) 
//    return (A(0,0)*A(1,1)-A(0,1)*A(1,0));
//  else {
//    Array<double,2> LU(m,m);
//    Array<int,1> ipiv(m);
//    int info;
//    LU = A;
//    // Do LU factorization
//    F77_DGETRF (&m, &n, LU.data(), &m, ipiv.data(), &info);
//    double det = 1.0;
//    int numPerm = 0;
//    for (int i=0; i<m; i++) {
//      det *= LU(i,i);
//      numPerm += (ipiv(i) != (i+1));
//    }
//    if (numPerm & 1)
//      det *= -1.0;
//    
//    return det;
//  }
//}

//std::complex<double> 
//Determinant (const Array<std::complex<double>,2> &A)
//{
//  int m = A.rows();
//  int n = A.cols();
//  assert (m == n);  // Cannot take a determinant of a non-square
//		    // matrix
//  if (A.rows() == 1)
//    return (A(0,0));
//  if (A.rows() == 2) 
//    return (A(0,0)*A(1,1)-A(0,1)*A(1,0));
//  else {
//    Array<std::complex<double>,2> LU(m,m);
//    Array<int,1> ipiv(m);
//    int info;
//    LU = A;
//    // Do LU factorization
//    F77_ZGETRF (&m, &n, LU.data(), &m, ipiv.data(), &info);
//    std::complex<double> det = 1.0;
//    int numPerm = 0;
//    for (int i=0; i<m; i++) {
//      det *= LU(i,i);
//      numPerm += (ipiv(i) != (i+1));
//    }
//    if (numPerm & 1)
//      det *= -1.0;
//    
//    return det;
//  }
//}


// The cofactors of A are given by 
// cof(A) = det(A) transpose(A^{-1})
//void Cofactors (const Array<double,2> &A, 
//		Array<double,2> &cof,
//		Array<double,2> &scratch)
//{
//  const int maxSize = 2000;
//  int m = A.rows();
//  int n = A.cols();
//  assert (m == n);  // Cannot take cofactors of a non-square matrix
//  assert (A.cols() < maxSize);
//  int  ipiv[maxSize];
//  int info;
//  // Copy and transpose for FORTRAN ordering
//  for (int i=0; i<m; i++)
//    for (int j=0; j<m; j++)
//      scratch(i,j) = A(j,i);
//  // Do LU decomposition
//  F77_DGETRF (&m, &n, scratch.data(), &m, ipiv, &info);
//  // Now scratch contains LU matrix in fortran ordering with pivots in ipiv
//  // Put identity matrix in cof
//  cof = 0.0;
//  for (int i=0; i<m; i++)
//    cof(i,i) = 1.0;
//  int numPerm = 0;
//  // Now apply permutation matrix to cof
//  for (int row=0; row<m; row++) {
//    int ip = ipiv[row]-1;
//    if (ip != row) { 
//      SwapRow(cof, row, ip);
//      numPerm++;
//    }
//  }
//}


//void SVdecomp (Array<double,2> &A,
//	       Array<double,2> &U, Array<double,1> &S,
//	       Array<double,2> &V)
//{
//  int M = A.rows();
//  int N = A.cols();
//  Array<double,2> Atrans(M,N);
//  // U will be Utrans after lapack call
//  U.resize(std::min(M,N),M);
//  V.resize(N,std::min(M,N));
//  
//  S.resize(std::min(N,M));
//  Atrans = A;

//  // Transpose U for FORTRAN ordering
//  Transpose(Atrans);
//  char JOBU  = 'S'; // return min (M,N) columns of U
//  char JOBVT = 'S'; // return min (M,N) columns of V
//  int LDA = M;
//  int LDU = M;
//  int LDVT = std::min(M,N);
//  int LWORK = 10 * std::max(3*std::min(M,N)+std::max(M,N),5*std::min(M,N));
//  Array<double,1> WORK(LWORK);
//  int INFO;

//  F77_DGESVD (&JOBU, &JOBVT, &M, &N, Atrans.data(), &LDA,
//	      S.data(), U.data(), &LDU, V.data(), &LDVT,
//	      WORK.data(), &LWORK, &INFO);
//  assert (INFO == 0);
//  // Transpose U to get back to C ordering
//  // V was really Vtrans so we don't need to transpose
//  Transpose(U);
//}


void SVdecomp (Array<std::complex<double>,2> &A, Array<std::complex<double>,2> &U, 
	       Array<double,1> &S, Array<std::complex<double>,2> &V)
{
  int M = A.rows();
  int N = A.cols();
  Array<std::complex<double>,2> Atrans(M,N);
  // U will be Utrans after lapack call
  U.resize(std::min(M,N),M);
  V.resize(N,std::min(M,N));
  
  S.resize(std::min(N,M));
  Atrans = A;

  // Transpose U for FORTRAN ordering
  Transpose(Atrans);
  char JOBU  = 'S'; // return min (M,N) columns of U
  char JOBVT = 'S'; // return min (M,N) columns of V
  int LDA = M;
  int LDU = M;
  int LDVT = std::min(M,N);
  int LWORK = 10 * std::max(3*std::min(M,N)+std::max(M,N),5*std::min(M,N));
  Array<std::complex<double>,1> WORK(LWORK);
  Array<double,1> WORK2(5*std::min(M,N));
  int INFO;

  F77_ZGESVD (&JOBU, &JOBVT, &M, &N, Atrans.data(), &LDA,
	      S.data(), U.data(), &LDU, V.data(), &LDVT,
	      WORK.data(), &LWORK, WORK2.data(), &INFO);
  assert (INFO == 0);
  // Transpose U to get back to C ordering
  // V was really Vtrans so we don't need to transpose
  Transpose(U);
}


//void PolarOrthogonalize (Array<std::complex<double>,2> &A)
//{
//  int M = A.rows();
//  int N = A.cols();
//  if (M != N) {
//    std::cerr << "Error:  nonsquare matrix in PolarOrthogonalize. Aborting.\n";
//    abort();
//  }
//  Array<std::complex<double>,2> U, V;
//  Array<double,1> S;

//  SVdecomp (A, U, S, V);
//  Transpose(V);
//  for (int i=0; i<V.rows(); i++)
//    for (int j=0; j<V.cols(); j++)
//      V(i,j) = conj(V(i,j));
//  A = U * V;
//}

//void SymmEigenPairs (const Array<scalar,2> &A, int NumPairs,
//		     Array<scalar,1> &Vals,
//		     Array<scalar,2> &Vectors)
//{
//  char JobType = 'V';    // Find eigenvectors and eignevalues
//  char Range   = 'I';    // Find eigenpairs in a range of indices
//  char UpperLower = 'U'; // Use upper triagle of A

//  int N   = A.rows();
//  double *Amat = new double[N*N];
//  int LDA = N;
//  double VL = 0.0;
//  double VU = 0.0;
//  int IL = 1;
//  int IU = NumPairs;
//  double AbsTolerance = -1.0;
//  int NumComputed;
//  double *EigVals = new double[N];
//  double *EigVecs = new double[N*NumPairs];
//  int LDEigVecs = N;
//  int *ISuppZ = new int[2*NumPairs];
//  int Info;

//  // First do workspace query
//  int Lwork = -1;
//  int LIwork = -1;
//  double WorkSize;
//  int IWorkSize;
//  
//  
//   F77_DSYEVR (&JobType, &Range, &UpperLower, &N, Amat, &LDA, &VL, &VU,
//	       &IL, &IU, &AbsTolerance, &NumComputed, EigVals, EigVecs,
//	       &LDEigVecs, ISuppZ, &WorkSize, &Lwork, &IWorkSize, &LIwork, 
//	       &Info);

//   // Now allocate WorkSpace;
//   Lwork = (int) floor(WorkSize+0.5);
//   LIwork = IWorkSize;
//   double *WorkSpace = new double[Lwork];
//   int *IWorkSpace = new int[LIwork];
//   
//   // Copy A int Amat
//   for (int row=0; row<N; row++)
//     for (int col=0; col<N; col++)
//       *(Amat+(col*N)+row) = A(row,col);
//  
//   F77_DSYEVR (&JobType, &Range, &UpperLower, &N, Amat, &LDA, &VL, &VU,
//	       &IL, &IU, &AbsTolerance, &NumComputed, EigVals, EigVecs,
//	       &LDEigVecs, ISuppZ, WorkSpace, &Lwork, IWorkSpace, &LIwork, 
//	       &Info);

//   if (Info !=0) 
//     std::cerr << "Lapack error in DSYEVR: " << Info << std::endl;

//   // Now copy over output of Vectors and Vals
//   Vals.resize(NumPairs);
//   Vectors.resize(NumPairs,N);
//   
//   for (int i=0; i<NumPairs; i++)
//     {
//       Vals(i) = *(EigVals+i);
//       for (int j=0; j<N; j++)
//	 Vectors(i,j) = *(EigVecs+(i*N)+j);
//     }

//   // Now free allocate memory
//   delete[] Amat;
//   delete[] EigVals;
//   delete[] EigVecs; 
//   delete[] WorkSpace; 
//   delete[] IWorkSpace; 
//   delete[] ISuppZ;
//}




//void SymmEigenPairs (const Array<std::complex<double>,2> &A, int NumPairs,
//		     Array<scalar,1> &Vals,
//		     Array<std::complex<double>,2> &Vectors)
//{
//  char JobType = 'V';    // Find eigenvectors and eignevalues
//  char Range   = 'I';    // Find eigenpairs in a range of indices
//  char UpperLower = 'U'; // Use upper triagle of A

//  int N   = A.rows();
//  std::complex<double> *Amat = new std::complex<double>[N*N];
//  int LDA = N;
//  double VL = 0.0;
//  double VU = 0.0;
//  int IL = 1;
//  int IU = NumPairs;
//  double AbsTolerance = 0.0;
//  int NumComputed;
//  double *EigVals = new double[N];
//  std::complex<double> *EigVecs = new std::complex<double>[N*NumPairs];
//  int LDEigVecs = N;
//  int *ISuppZ = new int[2*NumPairs];
//  int Info;

//  // First do workspace query
//  int Lwork = -1;
//  int LIwork = -1;
//  int LRwork = -1;
//  std::complex<double> WorkSize;
//  double RWorkSize;
//  int IWorkSize;
//  
//  
//   F77_ZHEEVR(&JobType, &Range, &UpperLower, &N, Amat, &LDA, &VL, &VU,
//	      &IL, &IU, &AbsTolerance, &NumComputed, EigVals, EigVecs,
//	      &LDEigVecs, ISuppZ, &WorkSize, &Lwork, &RWorkSize, &LRwork, 
//	      &IWorkSize, &LIwork, &Info);
////    fprintf (stderr, "WorkSize  = %1.8f\n", WorkSize.real());
////    fprintf (stderr, "RWorkSize = %1.8f\n", RWorkSize);
////    fprintf (stderr, "IWorkSize = %d\n", IWorkSize);

//   // Now allocate WorkSpace;
//   Lwork = (int) floor(WorkSize.real()+0.5);
//   LIwork = IWorkSize;
//   LRwork = (int) floor (RWorkSize+0.5);
//   std::complex<double> *WorkSpace = new std::complex<double>[Lwork];
//   double * RWorkSpace = new double[LRwork];
//   int *IWorkSpace = new int[LIwork];
//   
//   // Copy A int Amat
//   for (int row=0; row<N; row++)
//     for (int col=0; col<N; col++)
//       *(Amat+(col*N)+row) = A(row,col);
//   
//   F77_ZHEEVR(&JobType, &Range, &UpperLower, &N, Amat, &LDA, &VL, &VU,
//	      &IL, &IU, &AbsTolerance, &NumComputed, EigVals, EigVecs,
//	      &LDEigVecs, ISuppZ, WorkSpace, &Lwork, RWorkSpace, &LRwork, 
//	      IWorkSpace, &LIwork, &Info);

//   if (Info !=0) {
//     fprintf (stderr, "Lapack error in zheevr.  Exitting.\n");
//     exit(-1);
//   }

//   // Now copy over output of Vectors and Vals
//   Vals.resize(NumPairs);
//   Vectors.resize(NumPairs,N);
//   
//   for (int i=0; i<NumPairs; i++) {
//     Vals(i) = *(EigVals+i);
//     for (int j=0; j<N; j++)
//       Vectors(i,j) = *(EigVecs+(i*N)+j);
//   }

//   // Now free allocate memory
//   delete[] Amat;
//   delete[] EigVals;
//   delete[] EigVecs;
//   delete[] WorkSpace;
//   delete[] IWorkSpace;
//   delete[] ISuppZ;
//   

//}


///// This function returns the determinant of A and replaces A with its
///// cofactors.
//double 
//DetCofactors (Array<double,2> &A, Array<double,1> &work)
//{
//  const int maxN = 2000;
//  int ipiv[maxN];
//  int N = A.rows();
//  int M = A.cols();
//  assert (N == M);
//  assert (N <= maxN);
//  // First, transpose A for fortran ordering
////   for (int i=0; i<N; i++)
////     for (int j=0; j<i; j++) {
////       double tmp = A(i,j);
////       A(i,j) = A(j,i);
////       A(j,i) = tmp;
////     }
//  Transpose(A);
//  
//  int info;
//  // Do LU factorization
//  F77_DGETRF (&N, &M, A.data(), &N, ipiv, &info);
//  double det = 1.0;
//  int numPerm = 0;
//  for (int i=0; i<N; i++) {
//    det *= A(i,i);
//    numPerm += (ipiv[i] != (i+1));
//  }
//  if (numPerm & 1)
//    det *= -1.0;
//  
//  int lwork = work.size();
//  // Now, do inverse
//  F77_DGETRI (&N, A.data(), &N, ipiv, work.data(), &lwork, &info);

//  // Now, we have the transpose of Ainv.  Now, just multiply by det:
//  A = det * A;
//  // And we're done!
//  return det;
//}

//int 
//DetCofactorsWorksize(int N)
//{
//  double work;
//  double dummy;
//  int info;
//  int ipiv;
//  int lwork = -1;
//  
//  F77_DGETRI(&N, &dummy, &N, &ipiv, &work, &lwork, &info);

//  return ((int)ceil(work));
//}



/// This function returns the determinant of A and replaces A with its
/// cofactors.
//std::complex<double>
//ComplexDetCofactors (Array<std::complex<double>,2> &A, 
//		     Array<std::complex<double>,1> &work)
//{
//  const int maxN = 2000;
//  int ipiv[maxN];
//  int N = A.rows();
//  int M = A.cols();
//  assert (N == M);
//  assert (N <= maxN);
//  // First, transpose A for fortran ordering
////   for (int i=0; i<N; i++)
////     for (int j=0; j<i; j++) {
////       double tmp = A(i,j);
////       A(i,j) = A(j,i);
////       A(j,i) = tmp;
////     }
//  Transpose(A);
//  
//  int info;
//  // Do LU factorization
//  F77_ZGETRF (&N, &M, A.data(), &N, ipiv, &info);
//  std::complex<double> det = 1.0;
//  int numPerm = 0;
//  for (int i=0; i<N; i++) {
//    det *= A(i,i);
//    numPerm += (ipiv[i] != (i+1));
//  }
//  if (numPerm & 1)
//    det = -det;
//  
//  int lwork = work.size();
//  // Now, do inverse
//  F77_ZGETRI (&N, A.data(), &N, ipiv, work.data(), &lwork, &info);

//  // Now, we have the transpose of Ainv.  Now, just multiply by det:
//  A = det * A;
//  // And we're done!
//  return det;
//}

//int 
//ComplexDetCofactorsWorksize(int N)
//{
//  std::complex<double> work;
//  std::complex<double> dummy;
//  int info;
//  int ipiv;
//  int lwork = -1;
//  
//  F77_ZGETRI(&N, &dummy, &N, &ipiv, &work, &lwork, &info);

//  return ((int)ceil(work.real()));
//}

extern "C" void   F77_DGEMV (char *TRANS, const int *M, const int *N, 
			     double *alpha, const void *A, 
			     const int *LDA, const void *X, 
			     const int *INCX, double *beta, 
			     const void *Y, const int *INCY);


//void
//MatVecProd (Array<double,2> &A, Array<double,1> &x, Array<double,1> &Ax)
//{
//  assert (A.cols() == x.size());
//  assert (A.rows() == Ax.size());

//  double zero(0.0);
//  double one(1.0);
//  char trans = 'T';

//  int n = A.rows();
//  int m = A.cols();
//  int inc = 1;

//  F77_DGEMV(&trans, &m, &n, &one, A.data(), &m,
//	    x.data(), &inc, &zero, Ax.data(), &inc);
//}



//void CholeskyBig (Array<double,2> &A)
//{
//  int n=A.extent(0);
//  int lda=A.extent(1);
//  int info;
//  char upper='L';
//  F77_DPOTRF(&upper,&n,A.data(),&lda,&info);
//  for (int i=0;i<A.extent(0);i++)
//    for (int j=0;j<i;j++)
//      A(i,j)=0.0;
//  
//}

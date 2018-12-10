
#include<iterator>

#include"Configuration.h"
//#include <AFQMC/Matrix/DenseMatrix.h>
#include "AFQMC/config.h"
#include "Numerics/OhmmsBlas.h"
#include "Numerics/Blasf.h"
#include "AFQMC/Numerics/DenseMatrixOperations.h"

namespace qmcplusplus
{

namespace DenseMatrixOperators
{

bool symEigenSysAll(int N, double* A, int LDA, double* eigVal, double* eigVec, int LDV)
{

  std::vector<double> A0(N*N);
  for(int i=0; i<N; i++)
  for(int j=0; j<N; j++) A0[i*N+j] = A[i*LDA+j];

  char JOBZ('V');
  char RANGE('A');
  char UPLO('U');
  double VL=0;
  double VU=0;
  int IL=0;
  int IU=0;
  double ABSTOL=1e-8;//DLAMCH( 'Safe minimum' );
  int M; // output: total number of eigenvalues found
  std::vector<int> ISUPPZ(2*N);
  std::vector<double> WORK(1); // set with workspace query
  int LWORK=-1;
  std::vector<int> IWORK(1);    
  int LIWORK=-1;
  int INFO;
  const double one=1.0;
  const double zero=0.0;

  dsyevr(JOBZ, RANGE, UPLO, N, &(A0[0]), N, VL, VU, IL, IU, ABSTOL, M, eigVal, eigVec, LDV, &(ISUPPZ[0]), &(WORK[0]), LWORK, &(IWORK[0]), LIWORK, INFO);

  LWORK = int(WORK[0]);  
  WORK.resize(LWORK);
  LIWORK = int(IWORK[0]);  
  IWORK.resize(LIWORK);

  // remember that Z comes out transposed!!!
  dsyevr (JOBZ, RANGE, UPLO, N, A0.data(), N, VL, VU, IL, IU, ABSTOL, M, eigVal, eigVec, LDV, &(ISUPPZ[0]), &(WORK[0]), LWORK, &(IWORK[0]), LIWORK, INFO);

  if(INFO != 0) {
    app_error()<<" Problems with eigenvalue/eigenvector calculation during eigenSysAll; INFO: " <<INFO <<std::endl;  
    return false; 
  } 
  if(M != N) {
    app_error()<<" Problems with eigenvalue/eigenvector calculation during eigenSysAll. Found too few eigenvalues. M: " <<M <<std::endl;  
    return false; 
  } 
  // transpose Z 
  transpose<double>(N,eigVec,LDV); 
  return true;
}

bool symEigenSysAll(int N, std::complex<double>* A, int LDA, double* eigVal, std::complex<double>* eigVec, int LDV)
{

  std::vector<std::complex<double> > A0(N*N);
  // transposing matrix since lapack expects fortran ordering
  // not a problem since routine destroys input matrix and we want to preserve it anyway  
  for(int i=0; i<N; i++)
  for(int j=0; j<N; j++) A0[i*N+j] = A[j*LDA+i];

  char JOBZ('V');
  char RANGE('A');
  char UPLO('U');
  double VL=0;
  double VU=0;
  int IL=0;
  int IU=0;
  double ABSTOL=1e-8;//DLAMCH( 'Safe minimum' );
  int M; // output: total number of eigenvalues found
  std::vector<int> ISUPPZ(2*N);
  std::vector<std::complex<double> > WORK(1); // set with workspace query
  int LWORK=-1;
  std::vector<double > RWORK(1); // set with workspace query
  int LRWORK=-1;
  std::vector<int> IWORK(1);    
  int LIWORK=-1;
  int INFO;
  const std::complex<double> one=1.0;
  const std::complex<double> zero=0.0;


  zheevr (JOBZ, RANGE, UPLO, N, &(A0[0]), N, VL, VU, IL, IU, ABSTOL, M, eigVal, eigVec, LDV, &(ISUPPZ[0]), &(WORK[0]), LWORK, &(RWORK[0]), LRWORK, &(IWORK[0]), LIWORK, INFO);

  LWORK = int(WORK[0].real());  
  WORK.resize(LWORK);
  LRWORK = int(RWORK[0]);  
  RWORK.resize(LRWORK);
  LIWORK = int(IWORK[0]);  
  IWORK.resize(LIWORK);

  // remember that Z comes out transposed!!!
  zheevr (JOBZ, RANGE, UPLO, N, A0.data(), N, VL, VU, IL, IU, ABSTOL, M, eigVal, eigVec, LDV, &(ISUPPZ[0]), &(WORK[0]), LWORK, &(RWORK[0]), LRWORK, &(IWORK[0]), LIWORK, INFO);

  if(INFO != 0) {
    app_error()<<" Problems with eigenvalue/eigenvector calculation during eigenSysAll; INFO: " <<INFO <<std::endl;  
    return false; 
  } 
  if(M != N) {
    app_error()<<" Problems with eigenvalue/eigenvector calculation during eigenSysAll. Found too few eigenvalues. M: " <<M <<std::endl;  
    return false; 
  } 
  // transpose Z 
  transpose<std::complex<double> >(N,eigVec,LDV); 
  return true;

}

bool symEigenSysSelect(int N, std::complex<double>* A, int LDA, int neig, double* eigVal, bool getEigV, std::complex<double>* eigVec, int LDV)
{

  std::vector<std::complex<double> > A0(N*N);
  // transposing matrix since lapack expects fortran ordering
  // not a problem since routine destroys input matrix and we want to preserve it anyway  
  for(int i=0; i<N; i++)
  for(int j=0; j<N; j++) A0[i*N+j] = A[j*LDA+i];

  char JOBZ('V');
  if(!getEigV) JOBZ = 'N';
  char RANGE('I');
  char UPLO('U');
  double VL=0;
  double VU=0;
  int IL=1;
  int IU=neig;
  double ABSTOL=1e-8;//DLAMCH( 'Safe minimum' );
  int M; // output: total number of eigenvalues found
  std::vector<int> ISUPPZ(2*N);
  std::vector<std::complex<double> > WORK(1); // set with workspace query
  int LWORK=-1;
  std::vector<double > RWORK(1); // set with workspace query
  int LRWORK=-1;
  std::vector<int> IWORK(1);    
  int LIWORK=-1;
  int INFO;
  const std::complex<double> one=1.0;
  const std::complex<double> zero=0.0;

  zheevr (JOBZ, RANGE, UPLO, N, &(A0[0]), N, VL, VU, IL, IU, ABSTOL, M, eigVal, eigVec, LDV, &(ISUPPZ[0]), &(WORK[0]), LWORK, &(RWORK[0]), LRWORK, &(IWORK[0]), LIWORK, INFO);

  LWORK = int(WORK[0].real());  
  WORK.resize(LWORK);
  LRWORK = int(RWORK[0]);  
  RWORK.resize(LRWORK);
  LIWORK = int(IWORK[0]);  
  IWORK.resize(LIWORK);

  // remember that Z comes out transposed!!!
  zheevr (JOBZ, RANGE, UPLO, N, A0.data(), N, VL, VU, IL, IU, ABSTOL, M, eigVal, eigVec, LDV, &(ISUPPZ[0]), &(WORK[0]), LWORK, &(RWORK[0]), LRWORK, &(IWORK[0]), LIWORK, INFO);

  if(INFO != 0) {
    app_error()<<" Problems with eigenvalue/eigenvector calculation during eigenSysAll; INFO: " <<INFO <<std::endl;  
    return false; 
  } 
  if(M != neig) {
    app_error()<<" Problems with eigenvalue/eigenvector calculation during eigenSysAll. Found too few eigenvalues. M: " <<M <<std::endl;  
    return false; 
  } 
  return true;
}

bool symEigenSysSelect(int N, double* A, int LDA, int neig, double* eigVal, bool getEigV, double* eigVec, int LDV) 
{

  std::vector<double> A0(N*N);
  for(int i=0; i<N; i++)
  for(int j=0; j<N; j++) A0[i*N+j] = A[i*LDA+j];

  char JOBZ('V');
  if(!getEigV) JOBZ = 'N';
  char RANGE('I');
  char UPLO('U');
  double VL=0;
  double VU=0;
  int IL=1;
  int IU=neig;
  double ABSTOL=1e-8;//DLAMCH( 'Safe minimum' );
  int M; // output: total number of eigenvalues found
  std::vector<int> ISUPPZ(2*N);
  std::vector<double> WORK(1); // set with workspace query
  int LWORK=-1;
  std::vector<int> IWORK(1);    
  int LIWORK=-1;
  int INFO;
  const double one=1.0;
  const double zero=0.0;

  dsyevr(JOBZ, RANGE, UPLO, N, &(A0[0]), N, VL, VU, IL, IU, ABSTOL, M, eigVal, eigVec, LDV, &(ISUPPZ[0]), &(WORK[0]), LWORK, &(IWORK[0]), LIWORK, INFO);

  LWORK = int(WORK[0]);  
  WORK.resize(LWORK);
  LIWORK = int(IWORK[0]);  
  IWORK.resize(LIWORK);

  // remember that Z comes out transposed!!!
  dsyevr (JOBZ, RANGE, UPLO, N, A0.data(), N, VL, VU, IL, IU, ABSTOL, M, eigVal, eigVec, LDV, &(ISUPPZ[0]), &(WORK[0]), LWORK, &(IWORK[0]), LIWORK, INFO);

  if(INFO != 0) {
    app_error()<<" Problems with eigenvalue/eigenvector calculation during eigenSysAll; INFO: " <<INFO <<std::endl;  
    return false; 
  } 
  if(M != neig) {
    app_error()<<" Problems with eigenvalue/eigenvector calculation during eigenSysAll. Found too few eigenvalues. M: " <<M <<std::endl;  
    return false; 
  } 
  return true;

}

bool genHermitianEigenSysSelect(int N, std::complex<double>* A, int LDA, std::complex<double>* B, int LDB, int neig, double* eigVal, bool getEigV, std::complex<double>* eigVec, int LDV, int* IFAIL)
{

  std::vector<std::complex<double> > A0(N*N);
  std::vector<std::complex<double> > B0(N*N);
  // transposing matrix since lapack expects fortran ordering
  // not a problem since routine destroys input matrix and we want to preserve it anyway  
  for(int i=0; i<N; i++)
  for(int j=0; j<N; j++) A0[i*N+j] = A[j*LDA+i];
  for(int i=0; i<N; i++)
  for(int j=0; j<N; j++) B0[i*N+j] = B[j*LDB+i];

  int ITYPE = 1;
  char JOBZ('N');
  char RANGE('I');
  char UPLO('U');
  double VL=0;
  double VU=0;
  int IL=1;
  int IU=neig;
  double ABSTOL=1e-8;//DLAMCH( 'Safe minimum' );
  int M; // output: total number of eigenvalues found
  std::vector<int> ISUPPZ(2*N);
  std::vector<std::complex<double> > WORK(1); // set with workspace query
  int LWORK=-1;
  std::vector<double > RWORK(7*N); // set with workspace query
  std::vector<int> IWORK(5*N);    
  int INFO;
  const std::complex<double> one=1.0;
  const std::complex<double> zero=0.0;

  if(getEigV) JOBZ = 'V'; 

  zhegvx (ITYPE, JOBZ, RANGE, UPLO, N, &(A0[0]), N, &(B0[0]), N, VL, VU, IL, IU, ABSTOL, M, eigVal, eigVec, LDV, &(WORK[0]), LWORK, &(RWORK[0]), &(IWORK[0]), IFAIL, INFO);

  LWORK = int(WORK[0].real());  
  WORK.resize(LWORK);

  // remember that Z comes out transposed!!!
  zhegvx (ITYPE, JOBZ, RANGE, UPLO, N, &(A0[0]), N, &(B0[0]), N, VL, VU, IL, IU, ABSTOL, M, eigVal, eigVec, LDV, &(WORK[0]), LWORK, &(RWORK[0]), &(IWORK[0]), IFAIL, INFO);

  if(INFO != 0) {
    app_error()<<" Problems with generalized eigenvalue/eigenvector calculation during genHermitianEigenSysSelect; INFO: " <<INFO <<std::endl;  
    return false; 
  } 
  return true;
}

bool exponentiateHermitianMatrix(int N, std::complex<double>* A, int LDA, std::complex<double>* expA, int LDEXPA) 
{

  std::vector<std::complex<double> > A0(N*N);  // temporary storage for later
  std::vector<double> W(N);  // computed eigenvalues in ascending order
  std::vector<std::complex<double> > Z(N*N); // computed eigenvectors
  const std::complex<double> one=std::complex<double>(1.0,0.0);
  const std::complex<double> zero=std::complex<double>(0.0,0.0);

  if(!symEigenSysAll(N,A,LDA,W.data(),Z.data(),N)) {
    app_error()<<" Problems in call to eigSysAll in exponentiateHermitianMatrix. \n" <<std::endl;
    return false; 
  } 

  // always do this test
  for(int i=0; i<N; i++) 
  for(int j=0; j<N; j++) 
    expA[i*LDEXPA+j] = std::complex<double>(0.0);
  for(int i=0; i<N; i++) expA[i*LDEXPA+i] = W[i]; 
  // A0 = V*W
  product(N,N,N,one,Z.data(),N,expA,LDEXPA,zero,A0.data(),N); 
  // expA = A0*V^* = V*expA*V^*
  BLAS::gemm('C','N', N, N, N,
          one, &(Z[0]), N, A0.data(), N, 
          zero, expA, LDEXPA);
  // expA should be equal to A

  RealType s=0.0;
  for(int i=0; i<N; i++)
  for(int j=0; j<N; j++)  
    s += std::abs(A[i*LDA+j]-expA[i*LDEXPA+j]);
  if( std::abs(s) > 1e-8) {
    std::cerr<<std::endl <<std::endl <<" Error in reconstruction of A: " <<s <<std::endl <<std::endl; 
    return false; 
  }

  // now exp(A) = Z*exp(M)*Z^*, where A = Z*M*Z^* and M is the diagonal matrix of eigenvalues 
  
  for(int i=0; i<N; i++) 
  for(int j=0; j<N; j++) 
    expA[i*LDEXPA+j] = std::complex<double>(0.0);
  for(int i=0; i<N; i++) expA[i*LDEXPA+i] = std::exp(W[i]); 
  // A0 = V*expA
  product(N,N,N,one,Z.data(),N,expA,LDEXPA,zero,A0.data(),N);
  // expA = A0*V^* = V*expA*V^*
  BLAS::gemm('C','N', N, N, N,
          one, &(Z[0]), N, A0.data(), N, 
          zero, expA, LDEXPA);
  // expA should be equal to A

  return true;
}

bool exponentiateHermitianMatrix(int N, double* A, int LDA, double* expA, int LDEXPA) 
{

  std::vector<double> A0(N*N);  // temporary storage for later
  std::vector<double> W(N);  // computed eigenvalues in ascending order
  std::vector<double> Z(N*N); // computed eigenvectors
  const double one=1.0;
  const double zero=0.0;

  if(!symEigenSysAll(N,A,LDA,W.data(),Z.data(),N)) {
    app_error()<<" Problems in call to eigSysAll in exponentiateHermitianMatrix. \n" <<std::endl;
    return false; 
  } 

  // always do this test
  for(int i=0; i<N; i++) 
  for(int j=0; j<N; j++) 
    expA[i*LDEXPA+j] = 0.0;
  for(int i=0; i<N; i++) expA[i*LDEXPA+i] = W[i]; 
  // A0 = V*W
  product(N,N,N,one,Z.data(),N,expA,LDEXPA,zero,A0.data(),N);
  // expA = A0*V' = V*expA*'V
  BLAS::gemm('T','N', N, N, N,
          one, &(Z[0]), N, A0.data(), N, 
          zero, expA, LDEXPA);
  // expA should be equal to A

  RealType s=0.0;
  for(int i=0; i<N; i++)
  for(int j=0; j<N; j++)  
    s += std::abs(A[i*LDA+j]-expA[i*LDEXPA+j]);
  if( std::abs(s) > 1e-8) {
    std::cerr<<std::endl <<std::endl <<" Error in reconstruction of A: " <<s <<std::endl <<std::endl; 
    return false; 
  }

  // now exp(A) = Z*exp(M)*Z', where A = Z*M*Z' and M is the diagonal matrix of eigenvalues 
  
  for(int i=0; i<N; i++) 
  for(int j=0; j<N; j++) 
    expA[i*LDEXPA+j] = 0.0;
  for(int i=0; i<N; i++) expA[i*LDEXPA+i] = std::exp(W[i]); 
  // A0 = V*expA
  product(N,N,N,one,Z.data(),N,expA,LDEXPA,zero,A0.data(),N);
  // expA = A0*V' = V*expA*'V
  BLAS::gemm('T','N', N, N, N,
          one, &(Z[0]), N, A0.data(), N, 
          zero, expA, LDEXPA);
  // expA should be equal to A

  return true;
}

void GeneralizedGramSchmidt(std::complex<double>* A, int LDA, int nR, int nC)
{
  //  void zgeqrf( const int *M, const int *N, std::complex<double> *A, const int *LDA, std::complex<double> *TAU, std::complex<double> *WORK, const int *LWORK, int *INFO );
  //  void zungqr( const int *M, const int *N, const int *K, std::complex<double> *A, const int *LDA, std::complex<double> *TAU, std::complex<double> *WORK, const int *LWORK, int *INFO );
  //
 
  // temporary
  std::vector<std::complex<double> > AT(nR*nC);
  for(int i=0; i<nR; i++)
   for(int j=0; j<nC; j++)
    AT[ j*nR+i ] = A[ i*LDA+j ];  
  
  int K = std::min(nR,nC);
  std::vector<std::complex<double> > TAU(K),WORK(1);
  int info,lwork=-1; 

  zgeqrf( nR, nC, AT.data(), nR, TAU.data(), WORK.data(), lwork, info);

  lwork = int(WORK[0].real());
  WORK.resize(lwork);

  zgeqrf( nR, nC, AT.data(), nR, TAU.data(), WORK.data(), lwork, info);

  if(info != 0) {
    app_error()<<" Problems with QR decomposition; INFO: " <<info <<std::endl;
    APP_ABORT("Problems with QR decomposition. \n");
  }
  
  zungqr( nR, nC, K, AT.data(), nR, TAU.data(), WORK.data(), lwork, info);

  if(info != 0) {
    app_error()<<" Problems with QR decomposition (zungqr); INFO: " <<info <<std::endl;
    APP_ABORT("Problems with QR decomposition (zungqr). \n");
  }

  for(int i=0; i<nR; i++)
   for(int j=0; j<nC; j++)
    A[ i*LDA+j ] = AT[ j*nR+i ];  

}

} // namespace DenseMatrixOperators

} // namespace qmcplusplus



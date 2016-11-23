
#include<iterator>
#include<tuple>
#include<cassert>

#include "AFQMC/config.h"
#include "sys/sysinfo.h"

#if defined(HAVE_MKL)
 #include "mkl.h"
 #include "mkl_service.h"
 #include "mkl_solvers_ee.h"
#elif defined(HAVE_ESSL)

#endif

namespace qmcplusplus
{

namespace SparseMatrixOperators
{

// Performs a product between a sparse matrix stored in format s2D and a dense 
// matrix stored in c format
//   N: number of rows in B/C 
//   M: number of columns in B/C
//   LDB: leading dimension of B
//
//   For eack term in A, aik 
//   C(i,:) += aik * B(k.:) 
void product_SD(const IndexType K,
             const s2D<ComplexType>* A, const int nterms,
             ComplexType* B, const IndexType LDB,  
             ComplexType* C, IndexType LDC )
{
  ComplexType aik=0;
  IndexType ii=0,kk=0; 
  ComplexType* lit;
  ComplexType* rit;
  for(int cnt1=0; cnt1<nterms; cnt1++) {
    std::tie(ii,kk,aik) = *(A++); 
//    if(!(ii>=0 && ii<28 && kk>=0 && kk <28))
//      std::cout<<ii <<" " <<kk <<std::endl;
    assert(ii>=0 && ii<28 && kk>=0 && kk <28);
//std::cout<<" SD: " <<cnt1 <<" " <<ii <<" " <<kk <<" " <<aik <<std::endl;
    lit=C+ii*LDC;
    rit=B+kk*LDB;  
    for(int cnt2=0; cnt2<K; cnt2++)
      *(lit++) += aik*(*(rit++));  
  }

}

template<class T>
void product_SpMatV(int nrows,
             const T& A,
             const ComplexType* B,
             ComplexType* C )
{
#if defined(HAVE_MKL)
    char trans = 'N';
    mkl_cspblas_zcsrgemv (&trans, &nrows, A.values() , A.row_index(), A.column_data(),  B, C);
//#elif defined(HAVE_ESSL)

#else
  ComplexType zero = ComplexType(0,0);
  const ComplexType* val = A.values();
  const int* cols = A.column_data();
  int disp = (A.zero_base())?0:-1; 
  if( A.format() == 0) {  // CSR
    const int* rows = A.row_index();
    for(int nr=0; nr<nrows; nr++,C++,rows++) {
      *C=zero;
      for(int i=*rows; i<*(rows+1); i++,val++,cols++)
        *C += (*val) * ( *( B + (*cols) + disp) ); 
    }
  } else  {  // EESL: Compressed Matrix 
     
  }

#endif
}

template<class T>
void product_SpMatV(const int M, const int K,
             const ComplexType& alpha,
             const T& A,
             const ComplexType* B,
             const ComplexType& beta, 
             ComplexType* C )
{
#if defined(HAVE_MKL)
  char trans = 'N';
  char matdes[6];
  matdes[0] = 'G';
  matdes[3] = 'C';
  mkl_zcsrmv( &trans, &M, &K, &alpha, matdes, A.values() , A.column_data(),  A.row_index() ,  A.row_index()+1, B, &beta, C );
//#elif defined(HAVE_ESSL)

#else

  ComplexType zero = ComplexType(0,0);
  const ComplexType* val = A.values();
  const int* cols = A.column_data();
  int disp = (A.zero_base())?0:-1;
  if( A.format() == 0) {  // CSR
    const int* rows = A.row_index();
    for(int nr=0; nr<M; nr++,C++,rows++) {
      *C*=beta;
      for(int i=*rows; i<*(rows+1); i++,val++,cols++)
        *C += alpha * (*val) * ( *( B + (*cols) + disp) );
    }
  } else  {  // EESL: Compressed Matrix 

  }

#endif
}

template<class T>
void product_SpMatV(const int M, const int K,
             const RealType& alpha,
             const T& A,
             const RealType* B,
             const RealType& beta, 
             RealType* C )
{
#if defined(HAVE_MKL)
  char trans = 'N';
  char matdes[6];
  matdes[0] = 'G';
  matdes[3] = 'C';
  mkl_dcsrmv( &trans, &M, &K, &alpha, matdes, A.values() , A.column_data(),  A.row_index() ,  A.row_index()+1, B, &beta, C );
//#elif defined(HAVE_ESSL)

#else

  RealType zero = RealType(0);
  const RealType* val = A.values();
  const int* cols = A.column_data();
  int disp = (A.zero_base())?0:-1;
  if( A.format() == 0) {  // CSR
    const int* rows = A.row_index();
    for(int nr=0; nr<M; nr++,C++,rows++) {
      *C*=beta;
      for(int i=*rows; i<*(rows+1); i++,val++,cols++)
        *C += alpha * (*val) * ( *( B + (*cols) + disp) );
    }
  } else  {  // EESL: Compressed Matrix 

  }

#endif
}

void product_SpMatV(const int M, const int K,
             const RealType& alpha,
             const RealType* val,
             const int* col,
             const int* row,
             const RealType* B,
             const RealType& beta,
             RealType* C )
{
#if defined(HAVE_MKL)
  char trans = 'N';
  char matdes[6];
  matdes[0] = 'G';
  matdes[3] = 'C';
  mkl_dcsrmv( &trans, &M, &K, &alpha, matdes, val , col,  row ,  row+1, B, &beta, C );
#else
APP_ABORT("ERROR: product_SpMatV only implemented with MKL. \n");
#endif
}

void product_SpMatV(const int M, const int K,
             const ComplexType& alpha,
             const ComplexType* val,
             const int* col,
             const int* row,
             const ComplexType* B,
             const ComplexType& beta,
             ComplexType* C )
{
#if defined(HAVE_MKL)
  char trans = 'N';
  char matdes[6];
  matdes[0] = 'G';
  matdes[3] = 'C';
  mkl_zcsrmv( &trans, &M, &K, &alpha, matdes, val , col,  row ,  row+1, B, &beta, C );
#else
APP_ABORT("ERROR: product_SpMatV only implemented with MKL. \n");
#endif
}

void product_SpMatTV(const int M, const int K,
             const ComplexType& alpha,
             const ComplexType* val,
             const int* col,
             const int* row,
             const ComplexType* B,
             const ComplexType& beta,
             ComplexType* C )
{
#if defined(HAVE_MKL)
  char trans = 'T';
  char matdes[6];
  matdes[0] = 'G';
  matdes[3] = 'C';
  mkl_zcsrmv( &trans, &M, &K, &alpha, matdes, val , col,  row ,  row+1, B, &beta, C );
#else
APP_ABORT("ERROR: product_SpMatV only implemented with MKL. \n");
#endif
}

template<class T>
void product_SpMatTV(const int M, const int K,
             const ComplexType& alpha,
             const T& A,
             const ComplexType* B,
             const ComplexType& beta, 
             ComplexType* C )
{
#if defined(HAVE_MKL)
  char trans = 'T';
  char matdes[6];
  matdes[0] = 'G';
  matdes[3] = 'C';
  mkl_zcsrmv( &trans, &M, &K, &alpha, matdes, A.values() , A.column_data(),  A.row_index() ,  A.row_index()+1, B, &beta, C );
//#elif defined(HAVE_ESSL)

#else

APP_ABORT("ERROR: product_SpMatTV only implemented with MKL. \n");
  ComplexType zero = ComplexType(0,0);
  const ComplexType* val = A.values();
  const int* cols = A.column_data();
  int disp = (A.zero_base())?0:-1;
  if( A.format() == 0) {  // CSR
    const int* rows = A.row_index();
    for(int nr=0; nr<M; nr++,C++,rows++) {
      *C*=beta;
      for(int i=*rows; i<*(rows+1); i++,val++,cols++)
        *C += alpha * (*val) * ( *( B + (*cols) + disp) );
    }
  } else  {  // EESL: Compressed Matrix 

  }

#endif
}

template<class T>
void product_SpMatM(const int M, const int N, const int K,
             const ComplexType& alpha,
             const T& A,
             const ComplexType* B, const int ldb,
             const ComplexType& beta, 
             ComplexType* C, int ldc )
{
#if defined(HAVE_MKL)
  char trans = 'N';
  char matdes[6];
  matdes[0] = 'G';
  matdes[3] = 'C';
  mkl_zcsrmm( &trans, &M, &N, &K, &alpha, matdes, A.values() , A.column_data(),  A.row_index() ,  A.row_index()+1, B, &ldb, &beta, C, &ldc );

#else

APP_ABORT("ERROR: product_SpMatM only implemented with MKL. \n");

#endif
}

template<class T>
void product_SpMatM(const int M, const int N, const int K,
             const RealType& alpha,
             const T& A,
             const RealType* B, const int ldb,
             const RealType& beta,
             RealType* C, int ldc )
{
#if defined(HAVE_MKL)
  char trans = 'N';
  char matdes[6];
  matdes[0] = 'G';
  matdes[3] = 'C';
  mkl_dcsrmm( &trans, &M, &N, &K, &alpha, matdes, A.values() , A.column_data(),  A.row_index() ,  A.row_index()+1, B, &ldb, &beta, C, &ldc );

#else

APP_ABORT("ERROR: product_SpMatM only implemented with MKL. \n");

#endif
}

template<class T>
void product_SpMatM(const int M, const int N, const int K,
             const float& alpha,
             const T& A,
             const float* B, const int ldb,
             const float& beta,
             float* C, int ldc )
{
#if defined(HAVE_MKL)
  char trans = 'N';
  char matdes[6];
  matdes[0] = 'G';
  matdes[3] = 'C';
  mkl_scsrmm( &trans, &M, &N, &K, &alpha, matdes, A.values() , A.column_data(),  A.row_index() ,  A.row_index()+1, B, &ldb, &beta, C, &ldc );

#else

APP_ABORT("ERROR: product_SpMatM only implemented with MKL. \n");

#endif
}

// Performs a product between a sparse matrix stored in format s2D and a dense 
// matrix stored in c format
//   N: number of rows in B/C 
//   M: number of columns in B/C
//   LDB: leading dimension of B
//
//   For eack term in A, aik 
//   C(i,:) += aik * B(k.:) 
void product_SD(const IndexType K,
             const s2D<RealType>* A, const int nterms,
             ComplexType* B, const IndexType LDB,                                            
             ComplexType* C, IndexType LDC )
{ 
 
  register RealType aik;
  register IndexType ii,kk; 
  int cnt1,cnt2;
  ComplexType* lit;
  ComplexType* rit;
  for(cnt1=0; cnt1<nterms; cnt1++,A++) {
    std::tie(ii,kk,aik) = *A;   
    for(cnt2=0,lit=C+ii*LDC,rit=B+kk*LDB; cnt2<K; cnt2++,lit++,rit++)
      *lit += aik*(*rit);
  }

}

// Dot product between 2 sparse vectors
template<class T>
T product_SpVSpV(const int n1,const  int* indx1, const T* A1, const int n2, const int* indx2, const T* A2) {
  T res=T(0);
  int i=0, j=0;
  while( i<n1 && j<n2 ) {
    if( *(indx1+i) < *(indx2+j)   )
      ++i;
    else if( *(indx2+j) < *(indx1+i) ) 
      ++j;
    else {
      res += *(A1+i) * (*(A2+j)); 
      ++i;++j;
    }
  }
  return res;
}

bool sparseEigenSystem(RealSpMat &A, int& m0, RealType *eigval, RealType* eigVec, double Emin )
{

  if(A.cols() != A.rows()) {
    std::cerr<<"Problems in sparseEigenSystem: A matrix not squared. \n" <<std::endl;
    return false;
  }

  if(!A.isCompressed()) {
    std::cerr<<"Problems in sparseEigenSystem: A matrix not compressed. \n" <<std::endl;
    return false;
  }

#if defined(HAVE_MKL)

  char UPLO('F');
  int N = A.rows();

  /* Declaration of FEAST variables */
  MKL_INT       fpm[128];      /* Array to pass parameters to Intel MKL Extended Eigensolvers */

  //double        Emin, Emax;    /* Lower/upper bound of search interval [Emin,Emax] */
  double        Emax = 100;    /* Lower/upper bound of search interval [Emin,Emax] */

  double        epsout;        /* Relative error of the trace */
  MKL_INT       loop;          /* Number of refinement loop */
  //MKL_INT       M0 = m0_;            /* Initial guess for subspace dimension to be used */
  MKL_INT       M;             /* Total number of eigenvalues found in the interval */

  /* Declaration of local variables */
  MKL_INT       info;          /* Errors */

  std::vector<double> res(m0);  

  /* Step 1. Call  FEASTINIT to define the default values for the input FEAST parameters */
  feastinit(
      fpm /* OUT: Array is used to pass parameters to Intel MKL Extended Eigensolvers */
      );

  m0   = 10;
  M    = 10;
  loop = 0;
  info = 0;
  epsout = 0.0;
  Emin = 1.0;
  Emax = 2.0;

/*

   N = 11;
    MKL_INT       rows[12] = { 1, 5, 10, 16, 23, 30, 37, 44, 51, 57, 62, 66 };
    MKL_INT       cols[65] = {    1,   2,   3,   4,
                                  1,   2,   3,   4,   5,
                                  1,   2,   3,   4,   5,   6,
                                  1,   2,   3,   4,   5,   6,   7,
                                       2,   3,   4,   5,   6,   7,   8,
                                            3,   4,   5,   6,   7,   8,   9,
                                                 4,   5,   6,   7,   8,   9,  10,
                                                      5,   6,   7,   8,   9,  10,  11,
                                                           6,   7,   8,   9,  10,  11,
                                                                7,   8,   9,  10,  11,
                                                                     8,   9,  10,  11
                            };

    double        val[65] = {   5.0, 2.0, 1.0, 1.0,
                                2.0, 6.0, 3.0, 1.0, 1.0,
                                1.0, 3.0, 6.0, 3.0, 1.0, 1.0,
                                1.0, 1.0, 3.0, 6.0, 3.0, 1.0, 1.0,
                                     1.0, 1.0, 3.0, 6.0, 3.0, 1.0, 1.0,
                                          1.0, 1.0, 3.0, 6.0, 3.0, 1.0, 1.0,
                                               1.0, 1.0, 3.0, 6.0, 3.0, 1.0, 1.0,
                                                    1.0, 1.0, 3.0, 6.0, 3.0, 1.0, 1.0,
                                                         1.0, 1.0, 3.0, 6.0, 3.0, 1.0,
                                                              1.0, 1.0, 3.0, 6.0, 2.0,
                                                                   1.0, 1.0, 2.0, 5.0 };


  fpm[0] = 1;
  Emin = 3;
  Emax = 7;
  m0 = 11;
  M = 11;
*/

//  dfeast_scsrev(
//      &UPLO,   /* IN: UPLO = 'F', stores the full matrix */
//      &N,      /* IN: Size of the problem */
//      val,     /* IN: CSR matrix A, values of non-zero elements */
//      rows,    /* IN: CSR matrix A, index of the first non-zero element in row */
//      cols,    /* IN: CSR matrix A, columns indices for each non-zero element */
//      fpm,     /* IN/OUT: Array is used to pass parameters to Intel MKL Extended Eigensolvers */
//      &epsout, /* OUT: Relative error of on the trace */
//      &loop,   /* OUT: Contains the number of refinement loop executed */
//      &Emin,   /* IN: Lower bound of search interval */
//      &Emax,   /* IN: Upper bound of search interval */
//      &m0,     /* IN: The initial guess for subspace dimension to be used. */
//      eigval,       /* OUT: The first M entries of Eigenvalues */
//      eigVec,       /* IN/OUT: The first M entries of Eigenvectors */
//      &M,      /* OUT: The total number of eigenvalues found in the interval */
//      res.data(),     /* OUT: The first M components contain the relative residual vector */
//      &info    /* OUT: Error code */
//      );

//cout<<"Routine dfeast_scsrev returns code of ERROR: " <<info <<std::endl;
//return false;

  std::cout<<"\nEntering dfeast_scsrev routine." <<std::endl;
  std::cout<<"Default subspace size: " <<fpm[4] <<std::endl; 
  std::cout<<"Problem size: " <<N <<std::endl;
  std::cout<<"Available memory: ";

  struct sysinfo si;
  sysinfo(&si);
  si.freeram+=si.bufferram;
  std::cout<<int(si.freeram>>20) <<std::endl;

  fpm[0] = 1;
  fpm[4] = 1;

  /* Step 2. Solve the standard Ax = ex eigenvalue problem. */
  dfeast_scsrev(
      &UPLO,   /* IN: UPLO = 'F', stores the full matrix */
      &N,      /* IN: Size of the problem */
      A.values(),     /* IN: CSR matrix A, values of non-zero elements */
      A.row_data(),    /* IN: CSR matrix A, index of the first non-zero element in row */
      A.column_data(),    /* IN: CSR matrix A, columns indices for each non-zero element */
      fpm,     /* IN/OUT: Array is used to pass parameters to Intel MKL Extended Eigensolvers */
      &epsout, /* OUT: Relative error of on the trace */
      &loop,   /* OUT: Contains the number of refinement loop executed */
      &Emin,   /* IN: Lower bound of search interval */
      &Emax,   /* IN: Upper bound of search interval */
      &m0,     /* IN: The initial guess for subspace dimension to be used. */
      eigval,       /* OUT: The first M entries of Eigenvalues */
      eigVec,       /* IN/OUT: The first M entries of Eigenvectors */
      &M,      /* OUT: The total number of eigenvalues found in the interval */
      res.data(),     /* OUT: The first M components contain the relative residual vector */
      &info    /* OUT: Error code */
      );
  if ( info != 0 )
  {
      std::cerr<<"Routine dfeast_scsrev returns code of ERROR: " <<info <<std::endl;
      std::cout<<"Routine dfeast_scsrev returns code of ERROR: " <<info <<std::endl;
      return false;
  }
  m0 = M;
  return true;
#else
APP_ABORT("Error: sparseEigenSystem only implemented with MKL library. n");
  return false;
#endif

}

bool sparseEigenSystem(ComplexSpMat &A, int& m0, RealType *eigval, ComplexType* eigVec, double Emin )
{

  if(A.cols() != A.rows()) {
    std::cerr<<"Problems in sparseEigenSystem: A matrix not squared. \n" <<std::endl;
    return false;
  }

  if(!A.isCompressed()) {
    std::cerr<<"Problems in sparseEigenSystem: A matrix not compressed. \n" <<std::endl;
    return false;
  }

#if defined(HAVE_MKL)

  char UPLO('F');
  int N = A.rows();

  /* Declaration of FEAST variables */
  MKL_INT       fpm[128];      /* Array to pass parameters to Intel MKL Extended Eigensolvers */

  //double        Emin, Emax;    /* Lower/upper bound of search interval [Emin,Emax] */
  double        Emax = 1e6;    /* Lower/upper bound of search interval [Emin,Emax] */

  double        epsout;        /* Relative error of the trace */
  MKL_INT       loop;          /* Number of refinement loop */
  //MKL_INT       M0 = m0_;            /* Initial guess for subspace dimension to be used */
  MKL_INT       M;             /* Total number of eigenvalues found in the interval */

  /* Declaration of local variables */
  MKL_INT       info;          /* Errors */

  std::vector<double> res(m0);  

  /* Step 1. Call  FEASTINIT to define the default values for the input FEAST parameters */
  feastinit(
      fpm /* OUT: Array is used to pass parameters to Intel MKL Extended Eigensolvers */
      );

  std::cout<<"Entering zfeast_hcsrev routine. \n" <<std::endl;

  /* Step 2. Solve the standard Ax = ex eigenvalue problem. */
  zfeast_hcsrev(
      &UPLO,   /* IN: UPLO = 'F', stores the full matrix */
      &N,      /* IN: Size of the problem */
      A.values(),     /* IN: CSR matrix A, values of non-zero elements */
      A.row_data(),    /* IN: CSR matrix A, index of the first non-zero element in row */
      A.column_data(),    /* IN: CSR matrix A, columns indices for each non-zero element */
      fpm,     /* IN/OUT: Array is used to pass parameters to Intel MKL Extended Eigensolvers */
      &epsout, /* OUT: Relative error of on the trace */
      &loop,   /* OUT: Contains the number of refinement loop executed */
      &Emin,   /* IN: Lower bound of search interval */
      &Emax,   /* IN: Upper bound of search interval */
      &m0,     /* IN: The initial guess for subspace dimension to be used. */
      eigval,       /* OUT: The first M entries of Eigenvalues */
      eigVec,       /* IN/OUT: The first M entries of Eigenvectors */
      &M,      /* OUT: The total number of eigenvalues found in the interval */
      res.data(),     /* OUT: The first M components contain the relative residual vector */
      &info    /* OUT: Error code */
      );
  if ( info != 0 )
  {
      std::cerr<<"Routine zfeast_hcsrev returns code of ERROR: " <<info <<std::endl;
      return false;
  }
  m0 = M;
  return true;

#else
APP_ABORT("Error: sparseEigenSystem only implemented with MKL library. n");
  return false;
#endif

}


template
void product_SpMatV<ComplexSpMat>(const int nrows, const ComplexSpMat& A, const ComplexType* B, ComplexType* C);
template
void product_SpMatV<ComplexSMSpMat>(const int nrows, const ComplexSMSpMat& A, const ComplexType* B, ComplexType* C);

template
void product_SpMatV<ComplexSpMat>(const int M, const int K,
             const ComplexType& alpha,
             const ComplexSpMat& A,
             const ComplexType* B,
             const ComplexType& beta,
             ComplexType* C );

template
void product_SpMatTV<ComplexSpMat>(const int M, const int K,
             const ComplexType& alpha,
             const ComplexSpMat& A,
             const ComplexType* B,
             const ComplexType& beta,
             ComplexType* C );

template
void product_SpMatM<ComplexSpMat>(const int M, const int N, const int K,
             const ComplexType& alpha,
             const ComplexSpMat& A,
             const ComplexType* B, int ldb,
             const ComplexType& beta,
             ComplexType* C, int ldc );

template
void product_SpMatM<RealSpMat>(const int M, const int N, const int K,
             const RealType& alpha,
             const RealSpMat& A,
             const RealType* B, int ldb,
             const RealType& beta,
             RealType* C, int ldc );

template
void product_SpMatV<RealSMSpMat>(const int M, const int K,
             const RealType& alpha,
             const RealSMSpMat& A,
             const RealType* B,
             const RealType& beta,
             RealType* C );

template
void product_SpMatV<ComplexSMSpMat>(const int M, const int K,
             const ComplexType& alpha,
             const ComplexSMSpMat& A,
             const ComplexType* B,
             const ComplexType& beta,
             ComplexType* C );

template
void product_SpMatTV<ComplexSMSpMat>(const int M, const int K,
             const ComplexType& alpha,
             const ComplexSMSpMat& A,
             const ComplexType* B,
             const ComplexType& beta,
             ComplexType* C );

template
void product_SpMatM<ComplexSMSpMat>(const int M, const int N, const int K,
             const ComplexType& alpha,
             const ComplexSMSpMat& A,
             const ComplexType* B, int ldb,
             const ComplexType& beta,
             ComplexType* C, int ldc );

template
void product_SpMatM<RealSMSpMat>(const int M, const int N, const int K,
             const RealType& alpha,
             const RealSMSpMat& A,
             const RealType* B, int ldb,
             const RealType& beta,
             RealType* C, int ldc );

template
void product_SpMatM<SMSparseMatrix<float>>(const int M, const int N, const int K,
             const float& alpha,
             const SMSparseMatrix<float>& A,
             const float* B, int ldb,
             const float& beta,
             float* C, int ldc );

template
ComplexType product_SpVSpV(const int n1, const int* indx1, const ComplexType* A1, const int n2, const int* indx2, const ComplexType* A2);
template
RealType product_SpVSpV(const int n1, const int* indx1, const RealType* A1, const int n2, const int* indx2, const RealType* A2);

} // namespace SparseMatrixOperators

} // namespace qmcplusplus


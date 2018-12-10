
#include<iterator>
#include<tuple>
#include<cassert>
#include"AFQMC/Numerics/SparseMatrixOperations.h"

#include "AFQMC/config.h"
#ifdef PRINT_FREE_MEMORY
#include "sys/sysinfo.h"
#endif

/*
#if defined(HAVE_MKL)
 #include "mkl.h"
 #include "mkl_service.h"
 #include "mkl_solvers_ee.h"
#elif defined(HAVE_ESSL)

#endif
*/

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
//  feastinit(
//      fpm /* OUT: Array is used to pass parameters to Intel MKL Extended Eigensolvers */
//      );

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

#ifdef PRINT_FREE_MEMORY
  struct sysinfo si;
  sysinfo(&si);
  si.freeram+=si.bufferram;
  std::cout<<int(si.freeram>>20) <<std::endl;
#endif

  fpm[0] = 1;
  fpm[4] = 1;

  /* Step 2. Solve the standard Ax = ex eigenvalue problem. */
//  dfeast_scsrev(
//      &UPLO,   /* IN: UPLO = 'F', stores the full matrix */
//      &N,      /* IN: Size of the problem */
//      A.values(),     /* IN: CSR matrix A, values of non-zero elements */
//      A.row_data(),    /* IN: CSR matrix A, index of the first non-zero element in row */
//      A.column_data(),    /* IN: CSR matrix A, columns indices for each non-zero element */
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
//  feastinit(
//      fpm /* OUT: Array is used to pass parameters to Intel MKL Extended Eigensolvers */
//      );

  std::cout<<"Entering zfeast_hcsrev routine. \n" <<std::endl;

  /* Step 2. Solve the standard Ax = ex eigenvalue problem. */
//   zfeast_hcsrev(
//      &UPLO,   /* IN: UPLO = 'F', stores the full matrix */
//      &N,      /* IN: Size of the problem */
//      A.values(),     /* IN: CSR matrix A, values of non-zero elements */
//      A.row_data(),    /* IN: CSR matrix A, index of the first non-zero element in row */
//      A.column_data(),    /* IN: CSR matrix A, columns indices for each non-zero element */
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
ComplexType product_SpVSpV(const int n1, const int* indx1, const ComplexType* A1, const int n2, const int* indx2, const ComplexType* A2);
template
RealType product_SpVSpV(const int n1, const int* indx1, const RealType* A1, const int n2, const int* indx2, const RealType* A2);
template
RealType product_SpVSpVc(const int n1, const int* indx1, const RealType* A1, const int n2, const int* indx2, const RealType* A2);

} // namespace SparseMatrixOperators

} // namespace qmcplusplus


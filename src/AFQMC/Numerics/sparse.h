
#ifndef AFQMC_SPARSE_H
#define AFQMC_SPARSE_H

#include "AFQMC/Numerics/spblas.h"
#include<cassert>
#include<complex>

struct mySPBLAS
{
 
  template<typename T>
  inline static
  void csrmv(const char transa, const int M, const int K, const T alpha, const char *matdescra, const T* A, const int* indx, const int *pntrb, const int *pntre, const T* x, const T beta, T *y  )
  {
    assert(matdescra[0]=='G' && (matdescra[3]=='C' || matdescra[3]=='F'));
    int disp = (matdescra[3]=='C')?0:-1;
    if(transa=='n' || transa=='N') {  
      for(int nr=0; nr<M; nr++,y++,pntrb++,pntre++) {
        (*y) *= beta;
        for(int i=*pntrb; i<*pntre; i++) {
          if(*(indx+i)+disp >= K) continue;
          *y += alpha * (*(A+i)) * ( *( x + (*(indx+i)) + disp) );
        }
      }
    } else if(transa=='t' || transa=='T') {
      for(int k=0; k<K; k++) 
        (*(y+k)) *= beta;
      for(int nr=0; nr<M; nr++,pntrb++,pntre++,x++) {
        for(int i=*pntrb; i<*pntre; i++) {
          if(*(indx+i)+disp >= K) continue;
          *(y+(*(indx+i))+disp) += alpha * (*(A+i)) * (*x);
        }
      }
    } else if(transa=='h' || transa=='H') {
      for(int k=0; k<K; k++)
        (*(y+k)) *= beta;
      for(int nr=0; nr<M; nr++,pntrb++,pntre++,x++) {
        for(int i=*pntrb; i<*pntre; i++) {
          if(*indx+disp >= K) continue;
          *(y+(*(indx+i))+disp) += alpha * (*(A+i)) * (*x);
        }
      }
    }
  }

  template<typename T>
  inline static
  void csrmv(const char transa, const int M, const int K, const std::complex<T> alpha, const char *matdescra, const std::complex<T>* A, const int* indx, const int *pntrb, const int *pntre, const std::complex<T>* x, const std::complex<T> beta, std::complex<T> *y  )
  {
    assert(matdescra[0]=='G' && (matdescra[3]=='C')); // || matdescra[3]=='F'));
    int disp = (matdescra[3]=='C')?0:-1;
    if(transa=='n' || transa=='N') {
      for(int nr=0; nr<M; nr++,y++,pntrb++,pntre++) {
        (*y) *= beta;
        for(int i=*pntrb; i<*pntre; i++) {
          if(*(indx+i)+disp >= K) continue;
          *y += alpha * (*(A+i)) * ( *( x + (*(indx+i)) + disp) );
        }
      }
    } else if(transa=='t' || transa=='T') {
      for(int k=0; k<K; k++)
        (*(y+k)) *= beta;
      for(int nr=0; nr<M; nr++,pntrb++,pntre++,x++) {
        for(int i=*pntrb; i<*pntre; i++) {
          if(*(indx+i)+disp >= K) continue;
          *(y+(*(indx+i))+disp) += alpha * (*(A+i)) * (*x);
        }
      }
    } else if(transa=='h' || transa=='H') {
      for(int k=0; k<K; k++)
        (*(y+k)) *= beta;
      for(int nr=0; nr<M; nr++,pntrb++,pntre++,x++) {
        for(int i=*pntrb; i<*pntre; i++) {
          if(*indx+disp >= K) continue;
          *(y+(*(indx+i))+disp) += alpha * std::conj(*(A+i)) * (*x);
        }
      }
    }
  }

  template<typename T>
  inline static
  void csrmm(const char transa, const int M, const int N, const int K, const T alpha, const char *matdescra, const T *A, const int *indx, const int *pntrb, const int *pntre, const T *B, const int ldb, const T beta, T *C, const int ldc)
  {
    assert(matdescra[0]=='G' && (matdescra[3]=='C')); // || matdescra[3]=='F'));
    int disp = (matdescra[3]=='C')?0:-1;
    if(transa=='n' || transa=='N') {
      for(int nr=0; nr<M; nr++,pntrb++,pntre++,C+=ldc) {
        for(int i=0; i<N; i++)
          (*(C+i)) *= beta;
        for(int i=*pntrb; i<*pntre; i++) {
          if(*(indx+i)+disp >= K) continue;
          // at this point *(A+i) is A_rc, c=*(indx+i)+disp, *C is C(r,0)
          // C(r,:) = A_rc * B(c,:)
          const T* Bc = B+ldb*(*(indx+i)+disp);  
          T* Cr = C; 
          T Arc = alpha*(*(A+i)); 
          for(int k=0; k<N; k++, Cr++, Bc++) 
            *Cr += Arc * (*Bc);
        }
      }
    } else if(transa=='t' || transa=='T') {
      // not optimal, but simple
      for(int i=0; i<K; i++)
       for(int j=0; j<N; j++)
        (*(C+i*ldc+j)) *= beta;
      for(int nr=0; nr<M; nr++,pntrb++,pntre++) {
        for(int i=*pntrb; i<*pntre; i++, B+=ldb) {
          if(*(indx+i)+disp >= K) continue;
          // at this point *(A+i) is A_rc, c=*(indx+i)+disp
          // C(c,:) = A_rc * B(r,:)
          const T* Br = B;
          T* Cc = C+ldc*(*(indx+i)+disp);
          T Arc = alpha*(*(A+i));
          for(int k=0; k<N; k++, Cc++, Br++)
            *Cc += Arc * (*Br);
        }
      }
    } else if(transa=='h' || transa=='H') {
      // not optimal, but simple
      for(int i=0; i<K; i++)
       for(int j=0; j<N; j++)
        (*(C+i*ldc+j)) *= beta;
      for(int nr=0; nr<M; nr++,pntrb++,pntre++) {
        for(int i=*pntrb; i<*pntre; i++, B+=ldb) {
          if(*(indx+i)+disp >= K) continue;
          // at this point *(A+i) is A_rc, c=*(indx+i)+disp
          // C(c,:) = A_rc * B(r,:)
          const T* Br = B;
          T* Cc = C+ldc*(*(indx+i)+disp);
          T Arc = alpha*(*(A+i));
          for(int k=0; k<N; k++, Cc++, Br++)
            *Cc += Arc * (*Br);
        }
      }
    }
  }

  template<typename T>
  inline static
  void csrmm(const char transa, const int M, const int N, const int K, const std::complex<T> alpha, const char *matdescra, const std::complex<T> *A, const int *indx, const int *pntrb, const int *pntre, const std::complex<T> *B, const int ldb, const std::complex<T> beta, std::complex<T> *C, const int ldc)
  {
    assert(matdescra[0]=='G' && (matdescra[3]=='C')); // || matdescra[3]=='F'));
    int disp = (matdescra[3]=='C')?0:-1;
    if(transa=='n' || transa=='N') {
      for(int nr=0; nr<M; nr++,pntrb++,pntre++,C+=ldc) {
        for(int i=0; i<N; i++)
          (*(C+i)) *= beta;
        for(int i=*pntrb; i<*pntre; i++) {
          if(*(indx+i)+disp >= K) continue;
          // at this point *(A+i) is A_rc, c=*(indx+i)+disp, *C is C(r,0)
          // C(r,:) = A_rc * B(c,:)
          const std::complex<T>* Bc = B+ldb*(*(indx+i)+disp);  
          std::complex<T>* Cr = C; 
          std::complex<T> Arc = alpha*(*(A+i)); 
          for(int k=0; k<N; k++, Cr++, Bc++) 
            *Cr += Arc * (*Bc);
        }
      }
    } else if(transa=='t' || transa=='T') {
      // not optimal, but simple
      for(int i=0; i<K; i++)
       for(int j=0; j<N; j++)
        (*(C+i*ldc+j)) *= beta;
      for(int nr=0; nr<M; nr++,pntrb++,pntre++) {
        for(int i=*pntrb; i<*pntre; i++, B+=ldb) {
          if(*(indx+i)+disp >= K) continue;
          // at this point *(A+i) is A_rc, c=*(indx+i)+disp
          // C(c,:) = A_rc * B(r,:)
          const std::complex<T>* Br = B;
          std::complex<T>* Cc = C+ldc*(*(indx+i)+disp);
          std::complex<T> Arc = alpha*(*(A+i));
          for(int k=0; k<N; k++, Cc++, Br++)
            *Cc += Arc * (*Br);
        }
      }
    } else if(transa=='h' || transa=='H') {
      // not optimal, but simple
      for(int i=0; i<K; i++)
       for(int j=0; j<N; j++)
        (*(C+i*ldc+j)) *= beta;
      for(int nr=0; nr<M; nr++,pntrb++,pntre++) {
        for(int i=*pntrb; i<*pntre; i++, B+=ldb) {
          if(*(indx+i)+disp >= K) continue;
          // at this point *(A+i) is A_rc, c=*(indx+i)+disp
          // C(c,:) = A_rc * B(r,:)
          const std::complex<T>* Br = B;
          std::complex<T>* Cc = C+ldc*(*(indx+i)+disp);
          std::complex<T> Arc = alpha*std::conj(*(A+i));
          for(int k=0; k<N; k++, Cc++, Br++)
            *Cc += Arc * (*Br);
        }
      }
    }
  }

};

struct SPBLAS
{

  static const char TRANS = 'N';
  static const char NOTRANS = 'N';

  inline static
  void csrmv(const char transa, const int M, const int K, const float alpha, const char *matdescra, const float *A, const int* indx, const int *pntrb, const int *pntre, const float *x, const float beta, float *y  )
  {
#if defined(HAVE_MKL)
    mkl_scsrmv(transa,M,K,alpha,matdescra,A,indx,pntrb,pntre,x,beta,y);
#else
    mySPBLAS::csrmv(transa,M,K,alpha,matdescra,A,indx,pntrb,pntre,x,beta,y);
#endif
  }

  inline static
  void csrmv(const char transa, const int M, const int K, const double alpha, const char *matdescra, const double *A, const int* indx, const int *pntrb, const int *pntre, const double *x, const double beta, double *y  )
  {
#if defined(HAVE_MKL)
    mkl_dcsrmv(transa,M,K,alpha,matdescra,A,indx,pntrb,pntre,x,beta,y);
#else
    mySPBLAS::csrmv(transa,M,K,alpha,matdescra,A,indx,pntrb,pntre,x,beta,y);
#endif
  }

  inline static
  void csrmv(const char transa, const int M, const int K, const std::complex<float> alpha, const char *matdescra, const std::complex<float> *A, const int* indx, const int *pntrb, const int *pntre, const std::complex<float> *x, const std::complex<float> beta, std::complex<float> *y  )
  {
#if defined(HAVE_MKL)
    mkl_ccsrmv(transa,M,K,alpha,matdescra,A,indx,pntrb,pntre,x,beta,y);
#else
    mySPBLAS::csrmv(transa,M,K,alpha,matdescra,A,indx,pntrb,pntre,x,beta,y);
#endif
  }

  inline static
  void csrmv(const char transa, const int M, const int K, const std::complex<double> alpha, const char *matdescra, const std::complex<double> *A, const int* indx, const int *pntrb, const int *pntre, const std::complex<double> *x, const std::complex<double> beta, std::complex<double> *y  )
  {
#if defined(HAVE_MKL)
    mkl_zcsrmv(transa,M,K,alpha,matdescra,A,indx,pntrb,pntre,x,beta,y);
#else
    mySPBLAS::csrmv(transa,M,K,alpha,matdescra,A,indx,pntrb,pntre,x,beta,y);
#endif
  }

  inline static
  void csrmm(const char transa, const int M, const int N, const int K, const float alpha, const char *matdescra, const float *A, const int *indx, const int *pntrb, const int *pntre, const float *B, const int ldb, const float beta, float *C, const int ldc)
  {
#if defined(HAVE_MKL)
    mkl_scsrmm(transa,M,N,K,alpha,matdescra,A,indx,pntrb,pntre,B,ldb,beta,C,ldc);
#else
    mySPBLAS::csrmm(transa,M,N,K,alpha,matdescra,A,indx,pntrb,pntre,B,ldb,beta,C,ldc);
#endif
  }

  inline static
  void csrmm(const char transa, const int M, const int N, const int K, const std::complex<float> alpha, const char *matdescra, const std::complex<float> *A, const int *indx, const int *pntrb, const int *pntre, const std::complex<float> *B, const int ldb, const std::complex<float> beta, std::complex<float> *C, const int ldc)
  {
#if defined(HAVE_MKL)
    mkl_ccsrmm(transa,M,N,K,alpha,matdescra,A,indx,pntrb,pntre,B,ldb,beta,C,ldc);
#else
    mySPBLAS::csrmm(transa,M,N,K,alpha,matdescra,A,indx,pntrb,pntre,B,ldb,beta,C,ldc);
#endif
  }

  inline static
  void csrmm(const char transa, const int M, const int N, const int K, const double alpha, const char *matdescra, const double *A, const int *indx, const int *pntrb, const int *pntre, const double *B, const int ldb, const double beta, double *C, const int ldc)
  {
#if defined(HAVE_MKL)
    mkl_dcsrmm(transa,M,N,K,alpha,matdescra,A,indx,pntrb,pntre,B,ldb,beta,C,ldc);
#else
    mySPBLAS::csrmm(transa,M,N,K,alpha,matdescra,A,indx,pntrb,pntre,B,ldb,beta,C,ldc);
#endif
  }

  inline static
  void csrmm(const char transa, const int M, const int N, const int K, const std::complex<double> alpha, const char *matdescra, const std::complex<double> *A, const int *indx, const int *pntrb, const int *pntre, const std::complex<double> *B, const int ldb, const std::complex<double> beta, std::complex<double> *C, const int ldc)
  {
#if defined(HAVE_MKL)
    mkl_zcsrmm(transa,M,N,K,alpha,matdescra,A,indx,pntrb,pntre,B,ldb,beta,C,ldc);
#else
    mySPBLAS::csrmm(transa,M,N,K,alpha,matdescra,A,indx,pntrb,pntre,B,ldb,beta,C,ldc);
#endif
  }

};



#endif

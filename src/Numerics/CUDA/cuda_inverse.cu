#include <cstdio>
#include <unistd.h>
#include <sstream>
#include <vector>
#include <iostream>

template<typename T, int BS>
__global__ void
block_inverse (float* A, int N, int stride)
{
  __shared__ unsigned int ipiv[BS];
  __shared__ unsigned int kb;
  __shared__ T maxval[BS], mask[BS], pivotInv;
  __shared__ T Arowk[BS], Acolk[BS];
  ipiv[threadIdx.x] = threadIdx.x;
  mask[threadIdx.x] = (T)1.0;
  __syncthreads();

  unsigned int tid = threadIdx.x;

  for (int k=0; k<N; k++) {
    // First, find locate of maximum of kth column, excluding the
    // first k rows through the mask.
    maxval[tid] = mask[tid] * fabsf(A[tid*stride + k]);
    __syncthreads();
    for (int bs = BS>>1; bs>0; bs=bs>>1) {
      if (tid < bs) 
	maxval[tid] =  max(maxval[tid], maxval[tid+bs]);
      __syncthreads();
    }
    if ((mask[tid] * fabsf(A[tid*stride + k])) > 0.999* maxval[0]) {
      kb = tid;
      pivotInv = (T)1.0/A[tid*stride + k];
    }
    __syncthreads();
    // HACK HACK HACK
    //kb = k;
    //pivotInv = 1.0f/A[k*stride + k];
    //__syncthreads();

    // Now kb holds pivot row and pivot the value
    
    // Swap rows
    T tmp = A[k*stride+tid];
    A[k*stride +tid] = A[kb*stride+tid];
    A[kb*stride+tid] = tmp;
    
    // Swap pivot
    if (tid == 0) {
      int itmp = ipiv[kb];
      ipiv[kb] = ipiv[k];
      ipiv[k]  = itmp;
    }
    __syncthreads();

    // Col k update
    if (tid != k)
      A[stride*tid+k] = -pivotInv*A[stride*tid+k];
    else
      A[stride*k+k] = (T)0.0;
    __syncthreads();

    // Rank-1 update
    Arowk[tid] = A[stride*k   + tid];
    Acolk[tid] = A[stride*tid +   k];
    __syncthreads();
    for (int i=0; i<N; i++) 
      A[i*stride+tid] += Arowk[tid]*Acolk[i];
    __syncthreads();

    // Row k update
    if (tid != k) 
      A[k*stride+tid] *= pivotInv;
    else {
      A[k*stride+k] = pivotInv;
      mask[k] = 0.0;
    }
    __syncthreads();
  }
  // Finally, do backward pivoting
  for (int i=0; i<N; i++) {
    Arowk[tid] = A[i*stride+tid];
    __syncthreads();
    A[i*stride+ipiv[tid]] = Arowk[tid];
  }
}


template<typename T, int BS>
__device__ T
block_inverse2 (T A[BS][BS+1])
{
  __shared__ unsigned int ipiv[BS];
  __shared__ unsigned int kb;
  __shared__ T maxval[BS], mask[BS], pivotInv;
  __shared__ T Arowk[BS], Acolk[BS];
  bool write = threadIdx.y == 0;


  ipiv[threadIdx.x] = threadIdx.x;
  mask[threadIdx.x] = (T)1.0;
  __syncthreads();

  unsigned int tid = threadIdx.x;

  __shared__ T det;
  if (tid == 0)
    det = (T)1.0;


  for (int k=0; k<BS; k++) {
    // First, find locate of maximum of kth column, excluding the
    // first k rows through the mask.
    if (write) 
      maxval[tid] = mask[tid] * fabsf(A[tid][k]);
    __syncthreads();

    for (int bs = BS>>1; bs>0; bs=bs>>1) {
      if (tid < bs && write) 
	maxval[tid] =  max(maxval[tid], maxval[tid+bs]);
      __syncthreads();
    }

    if ((mask[tid] * fabsf(A[tid][k])) == maxval[0] && write) {
      kb = tid;
      pivotInv = (T)1.0/A[tid][k];
      if (kb == k)	det *= A[tid][k];
      else              det *= -A[tid][k];
    }


    __syncthreads();

    // Now kb holds pivot row and pivot the value
    
    // Swap rows
    if (write) {
      T tmp = A[k][tid];
      A[k][tid] = A[kb][tid];
      A[kb][tid] = tmp;
    }
    
    // Swap pivot
    if (tid == 0 && write) {
      int itmp = ipiv[kb];
      ipiv[kb] = ipiv[k];
      ipiv[k]  = itmp;
    }
    __syncthreads();

    // Col k update
    if (write) {
      if (tid != k)
	A[tid][k] = -pivotInv*A[tid][k];
      else
	A[k][k] = (T)0.0;
    }
    __syncthreads();

    // Rank-1 update
    Arowk[tid] = A[k][tid];
    Acolk[tid] = A[tid][k];
    __syncthreads();
    for (int i=0; i<BS; i+=blockDim.y) 
      A[i+threadIdx.y][tid] += Arowk[tid]*Acolk[i+threadIdx.y];
    __syncthreads();

    // Row k update
    if (write) {
      if (tid != k) 
	A[k][tid] *= pivotInv;
      else {
	A[k][k] = pivotInv;
	mask[k] = 0.0;
      }
    }
    __syncthreads();
  }
  // Finally, do backward pivoting
  for (int i=0; i<BS; i++) {
    if (write)
      Arowk[tid] = A[i][tid];
    __syncthreads();
    if (write)
      A[i][ipiv[tid]] = Arowk[tid];
  }
  return det;
}


template<typename T, int BS>
__device__ void block_mul2 (T A[BS][BS+1],
			    T B[BS][BS+1],
			    T C[BS][BS+1])
{
  int tid = threadIdx.x;
  for (int row=0; row<BS; row++)
    C[row][tid] = (T)0.0;
  __syncthreads();

  for (int k=0; k<BS; k++)
    for (int i=0; i<BS; i++)
      C[i][tid] += A[i][k]*B[k][tid];
}  


template<typename T, int BS>
__device__ void block_mul_add2 (T A[BS][BS+1],
				T B[BS][BS+1],
				T *C, int Cstride)
{
  int tid = threadIdx.x;
  __shared__ T Crow[2][BS];

  for (int i=0; i<BS; i+=blockDim.y) {
    Crow[threadIdx.y][tid] = C[(i+threadIdx.y)*Cstride + tid];
    for (int k=0; k<BS; k++) 
      Crow[threadIdx.y][tid] += A[i+threadIdx.y][k]*B[k][tid];
    C[(i+threadIdx.y)*Cstride + tid] = Crow[threadIdx.y][tid];
  }
}  

template<typename T, int BS>
__device__ void block_mul_set2 (T A[BS][BS+1],
				T B[BS][BS+1],
				T *C, int Cstride)
{
  int tid = threadIdx.x;
  __shared__ T Crow[2][BS];


  for (int i=0; i<BS; i+=blockDim.y) {
    Crow[threadIdx.y][tid] = (T)0.0;
    for (int k=0; k<BS; k++) 
      Crow[threadIdx.y][tid] += A[i+threadIdx.y][k]*B[k][tid];
    C[(i+threadIdx.y)*Cstride + tid] = Crow[threadIdx.y][tid];
  }
}  



template<typename T, int BS>
__device__ T
block_inverse1 (T A[BS][BS+1])
{
  __shared__ unsigned int ipiv[BS];
  __shared__ unsigned int kb;
  __shared__ T maxval[BS], mask[BS], pivotInv;
  __shared__ T Arowk[BS], Acolk[BS];
  ipiv[threadIdx.x] = threadIdx.x;
  mask[threadIdx.x] = (T)1.0;
  __syncthreads();

  unsigned int tid = threadIdx.x;

  __shared__ T det;
  if (tid == 0)
    det = (T)1.0;


  for (int k=0; k<BS; k++) {
    // First, find locate of maximum of kth column, excluding the
    // first k rows through the mask.
    maxval[tid] = mask[tid] * fabsf(A[tid][k]);
    __syncthreads();

    for (int bs = BS>>1; bs>0; bs=bs>>1) {
      if (tid < bs) 
	maxval[tid] =  max(maxval[tid], maxval[tid+bs]);
      __syncthreads();
    }

    if ((mask[tid] * fabsf(A[tid][k])) == maxval[0]) {
      kb = tid;
      pivotInv = (T)1.0/A[tid][k];
      if (kb == k)	det *= A[tid][k];
      else              det *= -A[tid][k];
    }


    __syncthreads();

    // Now kb holds pivot row and pivot the value
    
    // Swap rows
    T tmp = A[k][tid];
    A[k][tid] = A[kb][tid];
    A[kb][tid] = tmp;
    
    // Swap pivot
    if (tid == 0) {
      int itmp = ipiv[kb];
      ipiv[kb] = ipiv[k];
      ipiv[k]  = itmp;
    }
    __syncthreads();

    // Col k update
    if (tid != k)
      A[tid][k] = -pivotInv*A[tid][k];
    else
      A[k][k] = (T)0.0;
    __syncthreads();

    // Rank-1 update
    Arowk[tid] = A[k][tid];
    Acolk[tid] = A[tid][k];
    __syncthreads();
    for (int i=0; i<BS; i++) 
      A[i][tid] += Arowk[tid]*Acolk[i];
    __syncthreads();

    // Row k update
    if (tid != k) 
      A[k][tid] *= pivotInv;
    else {
      A[k][k] = pivotInv;
      mask[k] = 0.0;
    }
    __syncthreads();
  }
  // Finally, do backward pivoting
  for (int i=0; i<BS; i++) {
    Arowk[tid] = A[i][tid];
    __syncthreads();
    A[i][ipiv[tid]] = Arowk[tid];
  }
  return det;
}


template<typename T, int BS>
__device__ void block_mul (T A[BS][BS+1],
			   T B[BS][BS+1],
			   T C[BS][BS+1])
{
  int tid = threadIdx.x;
  for (int row=0; row<BS; row++)
    C[row][tid] = (T)0.0;
  __syncthreads();

  for (int k=0; k<BS; k++)
    for (int i=0; i<BS; i++)
      C[i][tid] += A[i][k]*B[k][tid];
}  


template<typename T, int BS>
__device__ void block_mul_add (T A[BS][BS+1],
			       T B[BS][BS+1],
			       T *C, int Cstride)
{
  int tid = threadIdx.x;
  __shared__ T Crow[BS];

  for (int i=0; i<BS; i++) {
    Crow[tid] = C[i*Cstride + tid];
    for (int k=0; k<BS; k++) 
      Crow[tid] += A[i][k]*B[k][tid];
    C[i*Cstride + tid] = Crow[tid];
  }
}  

template<typename T, int BS>
__device__ void block_mul_set (T A[BS][BS+1],
			       T B[BS][BS+1],
			       T *C, int Cstride)
{
  int tid = threadIdx.x;
  __shared__ T Crow[BS];


  for (int i=0; i<BS; i++) {
    Crow[tid] = (T)0.0;
    for (int k=0; k<BS; k++) 
      Crow[tid] += A[i][k]*B[k][tid];
    C[i*Cstride + tid] = Crow[tid];
  }
}  




template<typename T, int BS>
__global__ void
inverse (T* A, T* work, int N, int stride)
{
  T *Atmp = work;
  T *pivot_tmp = work+N*stride;

  __shared__ T pivot[BS][BS+1], in[BS][BS+1];
  int NB = N/BS;
  if (N%BS) NB++;
  int tid = threadIdx.x;


  for (int kb=0; kb<NB; kb++) {
    // load pivot block
    int row = kb*BS;
    for (int j=0; j<BS; j++)
      if (row+tid < N)
	pivot[j][tid] = A[(row+j)*stride + row+tid];
    
    // invert pivot
    block_inverse1<T,BS> (pivot);

    // Column scaling
    int col = kb*BS;
    for (int jb=0; jb < NB; jb++) {
      int row = jb*BS;
      if (kb != jb) {
    	for (int j=0; j<BS; j++)
    	  in[j][tid] = -A[(row+j)*stride + col + tid];
    	block_mul_set<T,BS>(in, pivot, A+row*stride+col, stride);
      }
      else {
    	for (int j=0; j<BS; j++)
    	  A[(row+j)*stride + col+tid] = (T)0.0;
      }
    }	

    // Save pivot to global memory here!
    // We use it for temporary space in the rank-1 update
    for (int j=0; j<BS; j++)
      pivot_tmp[j*BS+tid] = pivot[j][tid];


    // Copy Ato Atmp
    for (int ib=0; ib<NB; ib++)
      for (int row=0; row<N; row++)
    	Atmp[row*stride+ib*BS+tid] =  A[row*stride+ib*BS+tid];
    
    // Rank-1 update
    for (int ib=0; ib < NB; ib++) {
      for (int i=0; i<BS; i++)
    	in[i][tid] = A[(ib*BS+i)*stride + kb*BS + tid];
      for (int jb=0; jb<NB; jb++) {
    	for (int i=0; i<BS; i++) {
    	  pivot[i][tid] = A[(kb*BS+i)*stride + jb*BS + tid];
    	  // Atmp[(ib*BS+i)*stride + (jb*BS+tid)] = 
    	  //   A[(ib*BS+i)*stride + (jb*BS+tid)];
    	}
    	block_mul_add<T,BS>(in, pivot,  Atmp+(ib*BS)*stride + jb*BS,
    			    stride);
      }
    }
    // Copy Atmp back to A
    for (int ib=0; ib<NB; ib++)
      for (int row=0; row<N; row++)
    	A[row*stride+ib*BS+tid] =  Atmp[row*stride+ib*BS+tid];

    // Restore pivot from global memory here!
    for (int j=0; j<BS; j++)
      pivot[j][tid] = pivot_tmp[j*BS+tid];

    // Row-scaling
    for (int ib=0; ib<NB; ib++) {
      int row = kb*BS;
      int col = ib*BS;
      if (kb != ib) {
    	for (int j=0; j<BS; j++)
    	  in[j][tid] = A[(row+j)*stride + col+tid];
    	block_mul_set<T,BS>(pivot, in, A+row*stride+col, stride);
      }
      else {
    	for (int j=0; j<BS; j++) 
    	  A[(row+j)*stride + col+tid] = pivot[j][tid];
      }
    }	
  }
}



template<typename T, int BS>
__global__ void
inverse_many (T **A_list, T **work_list, int N, int stride)
{
  int tid = threadIdx.x;
  __shared__ T *A, *work;
  if (tid == 0 && threadIdx.y == 0) {
    A    = A_list[blockIdx.x];
    work = work_list[blockIdx.x];
  }
  __syncthreads();

  T *Atmp = work;
  T *pivot_tmp = work+N*stride;

  __shared__ T pivot[BS][BS+1], in[BS][BS+1];
  int NB = N/BS;
  if (N%BS) NB++;


  for (int kb=0; kb<NB; kb++) {
    // load pivot block
    int row = kb*BS;
    for (int j=0; j<BS; j++)
      if (row+tid < N && threadIdx.y == 0)
	pivot[j][tid] = A[(row+j)*stride + row+tid];
    
    // invert pivot
    block_inverse2<T,BS> (pivot);

    // Column scaling
    int col = kb*BS;
    for (int jb=0; jb < NB; jb++) {
      int row = jb*BS;
      if (kb != jb) {
	if (threadIdx.y == 0)
	  for (int j=0; j<BS; j++)
	    in[j][tid] = -A[(row+j)*stride + col + tid];
	__syncthreads();
    	block_mul_set2<T,BS>(in, pivot, A+row*stride+col, stride);
      }
      else if (threadIdx.y == 0)
    	for (int j=0; j<BS; j++)
    	  A[(row+j)*stride + col+tid] = (T)0.0;
    }	

    // Save pivot to global memory here!
    // We use it for temporary space in the rank-1 update
    if (threadIdx.y == 0) {
      for (int j=0; j<BS; j++)
	pivot_tmp[j*BS+tid] = pivot[j][tid];

      // Copy Ato Atmp
      for (int ib=0; ib<NB; ib++)
	for (int row=0; row<N; row++)
	  Atmp[row*stride+ib*BS+tid] =  A[row*stride+ib*BS+tid];
    }
    
    __syncthreads();
    // Rank-1 update
    for (int ib=0; ib < NB; ib++) {
      if (threadIdx.y == 0)
	for (int i=0; i<BS; i++)
	  in[i][tid] = A[(ib*BS+i)*stride + kb*BS + tid];
      for (int jb=0; jb<NB; jb++) {
	if (threadIdx.y == 0)
	  for (int i=0; i<BS; i++) 
	    pivot[i][tid] = A[(kb*BS+i)*stride + jb*BS + tid];
	__syncthreads();
    	block_mul_add2<T,BS>(in, pivot,  Atmp+(ib*BS)*stride + jb*BS,
    			    stride);
	__syncthreads();
      }
    }
    // Copy Atmp back to A
    if (threadIdx.y == 0)
      for (int ib=0; ib<NB; ib++)
	for (int row=0; row<N; row++)
	  A[row*stride+ib*BS+tid] =  Atmp[row*stride+ib*BS+tid];

    // Restore pivot from global memory here!
    if (threadIdx.y == 0)
      for (int j=0; j<BS; j++)
	pivot[j][tid] = pivot_tmp[j*BS+tid];

    // Row-scaling
    for (int ib=0; ib<NB; ib++) {
      int row = kb*BS;
      int col = ib*BS;
      if (kb != ib) {
	if (threadIdx.y == 0)
	  for (int j=0; j<BS; j++)
	    in[j][tid] = A[(row+j)*stride + col+tid];
	__syncthreads();
    	block_mul_set2<T,BS>(pivot, in, A+row*stride+col, stride);
      }
      else {
	if (threadIdx.y == 0)
	  for (int j=0; j<BS; j++) 
	    A[(row+j)*stride + col+tid] = pivot[j][tid];
      }
    }	
  }
}

#define MAX_BLOCKS 32

template<typename T, int BS>
__global__ void
inverse_many_pivot (T **A_list, T **work_list, int N, int stride)
{
  int tid = threadIdx.x;
  __shared__ T *A, *work;
  T maxdet, blockdet, det;
  __shared__ int ipiv[MAX_BLOCKS];

  if (tid == 0) {
    A    = A_list[blockIdx.x];
    work = work_list[blockIdx.x];
    det = (T)1.0;
  }
  ipiv[tid] = tid;
  ipiv[tid+BS] = tid+BS;
  __syncthreads();

  T *Atmp = work;
  T *pivot_tmp = work+N*stride;

  __shared__ T pivot[BS][BS+1], in[MAX_BLOCKS][BS+1];


  int NB = N/BS;
  if (N%BS) NB++;

  for (int kb=0; kb<NB; kb++) {
    int imax = kb;
    maxdet = (T)0.0;
    // Find pivot block
    for (int block=kb; block<NB; block++) {
      // load pivot block
      int row = block*BS;
      for (int j=0; j<BS; j++)
    	if (row+tid < N)
    	  in[j][tid] = A[(row+j)*stride + kb*BS + tid];
      __syncthreads();
      // invert pivot
      blockdet = block_inverse1<T,BS> (in);
      __syncthreads();
      if (fabs(blockdet) > fabs(maxdet)) {
      //if (block == kb) {
    	imax = block;
    	maxdet = blockdet;
    	for (int j=0; j<BS; j++)
    	  pivot[j][tid] = in[j][tid];
      }
    }
	
    // Now, swap row blocks
    for (int j=0; j<BS; j++) {
      int rowa = kb   * BS + j;
      int rowb = imax * BS + j;
      for (int n=0; n<NB; n++) {
    	int col = n*BS + tid;
    	T tmp = A[rowa*stride + col];
    	__syncthreads();
    	A[rowa*stride + col] = A[rowb*stride + col];
    	__syncthreads();
    	A[rowb*stride + col] = tmp;
      }
    }
    if (tid == 0) {
      int tmp = ipiv[kb];
      ipiv[kb] = ipiv[imax];
      ipiv[imax] = tmp;
      if (imax == kb)
    	det *= blockdet;
      else
    	det *= -blockdet;
    }

    // Column scaling
    int col = kb*BS;
    for (int jb=0; jb < NB; jb++) {
      int row = jb*BS;
      if (kb != jb) {
    	for (int j=0; j<BS; j++)
    	  in[j][tid] = -A[(row+j)*stride + col + tid];
    	block_mul_set<T,BS>(in, pivot, A+row*stride+col, stride);
      }
      else {
    	for (int j=0; j<BS; j++)
    	  A[(row+j)*stride + col+tid] = (T)0.0;
      }
    }	

    // Save pivot to global memory here!
    // We use it for temporary space in the rank-1 update
    for (int j=0; j<BS; j++)
      pivot_tmp[j*BS+tid] = pivot[j][tid];


    // Copy Ato Atmp
    for (int ib=0; ib<NB; ib++)
      for (int row=0; row<N; row++)
    	Atmp[row*stride+ib*BS+tid] =  A[row*stride+ib*BS+tid];
    
    // Rank-1 update
    for (int ib=0; ib < NB; ib++) {
      for (int i=0; i<BS; i++)
    	in[i][tid] = A[(ib*BS+i)*stride + kb*BS + tid];
      for (int jb=0; jb<NB; jb++) {
    	for (int i=0; i<BS; i++) 
    	  pivot[i][tid] = A[(kb*BS+i)*stride + jb*BS + tid];
    	block_mul_add<T,BS>(in, pivot,  Atmp+(ib*BS)*stride + jb*BS,
    			    stride);
      }
    }
    // Copy Atmp back to A
    for (int ib=0; ib<NB; ib++)
      for (int row=0; row<N; row++)
    	A[row*stride+ib*BS+tid] =  Atmp[row*stride+ib*BS+tid];

    // Restore pivot from global memory here!
    for (int j=0; j<BS; j++)
      pivot[j][tid] = pivot_tmp[j*BS+tid];

    // Row-scaling
    for (int ib=0; ib<NB; ib++) {
      int row = kb*BS;
      int col = ib*BS;
      if (kb != ib) {
    	for (int j=0; j<BS; j++)
    	  in[j][tid] = A[(row+j)*stride + col+tid];
    	block_mul_set<T,BS>(pivot, in, A+row*stride+col, stride);
      }
      else {
    	for (int j=0; j<BS; j++) 
    	  A[(row+j)*stride + col+tid] = pivot[j][tid];
      }
    }	
  }

  // Note:  the following assumes that N <= BS*BS
  // Finally, reverse pivoting
  for (int row=0; row<N; row++) {
    for (int block=0; block < NB; block++) 
      in[block][tid] = A[row*stride + BS*block + tid];
    __syncthreads();
    for (int block=0; block<NB; block++)
      A[row*stride + BS*ipiv[block] + tid] = in[block][tid];
  }
}


template<typename T, int BS>
__global__ void
inverse_many_naive (T **A_list, T **work_list, int N, int stride)
{
  int tid = threadIdx.x;
  __shared__ T *A;
  if (tid == 0) 
    A    = A_list[blockIdx.x];
  __syncthreads();

  __shared__ T colk[BS], rowk[BS];

  for (int k=0; k<N; k++) {
    rowk[tid] = A[k*stride+tid];
    __syncthreads();
    T pivInv = 1.0/rowk[k];
    __syncthreads();
    if (tid == k) rowk[k] = T();
    // Column scaling
    colk[tid] = (tid==k) ? T() : -pivInv*A[tid*stride+k];
    __syncthreads();
    A[tid*stride+k] = colk[tid];
    __syncthreads();
    // Rank-1 update
    for (int j=0; j<N; j++)
      A[j*stride+tid] += colk[j]*rowk[tid];
    __syncthreads();
    // Row scaling
    A[k*stride + tid] = pivInv * ((tid==k) ? 1.0 : rowk[tid]);
    __syncthreads();
  }
}




template<typename T, int BS>
__global__ void
inverse_many_naive_pivot (T **A_list, T **work_list, int N, int stride)
{
  int tid = threadIdx.x;
  __shared__ T *A;
  if (tid == 0) 
    A    = A_list[blockIdx.x];
  __syncthreads();
  __shared__ int kbar, ipiv[BS];  
  __shared__ T colk[BS], rowk[BS];
  __shared__ short imax[BS];

  ipiv[tid] = tid;
  __syncthreads();
  for (int k=0; k<N; k++) {
    // Find location of largest element in the column at or below k
    colk[tid] = (tid < k) ? 0.0 : fabs(A[tid*stride+k]);
    rowk[tid] = colk[tid];
    __syncthreads();
    int skip = 1<<((int)ceil(log2((double)BS)-1.0e-6)-1);
    imax[tid] = tid;
    __syncthreads();
    for (; skip>0; skip>>=1) {
      if (tid < skip && (tid+skip)<N) 
    	// colk[tid] = max(colk[tid],colk[tid+skip]);
	if (colk[tid+skip] > colk[tid]) {
	  imax[tid] = imax[tid+skip];
	  colk[tid] = colk[tid+skip];
	}
      __syncthreads();
    }
    if (tid == 0) {
      kbar = imax[0];
      int i = ipiv[kbar];
      ipiv[kbar] = ipiv[k];
      ipiv[k] = i;
    }
    
    // if (rowk[tid] == colk[0]) {
    //   kbar = tid;
    //   int i = ipiv[tid];
    //   ipiv[tid] = ipiv[k];
    //   ipiv[k]   = i;
    // }
    __syncthreads();
    // Swap rows
    rowk[tid] = A[kbar*stride + tid];
    colk[tid] = A[k   *stride + tid];
    __syncthreads();
    A[k   *stride + tid] = rowk[tid];
    A[kbar*stride + tid] = colk[tid];
    __syncthreads();
    T pivInv = 1.0/rowk[k];
    __syncthreads();
    if (tid == k) rowk[k] = T();
    // Column scaling
    colk[tid] = (tid==k) ? T() : -pivInv*A[tid*stride+k];
    __syncthreads();
    A[tid*stride+k] = colk[tid];
    __syncthreads();
    // Rank-1 update
    for (int j=0; j<N; j++)
      A[j*stride+tid] += colk[j]*rowk[tid];
    __syncthreads();
    // Row scaling
    A[k*stride + tid] = pivInv * ((tid==k) ? 1.0 : rowk[tid]);
    __syncthreads();
  }
  // Now, permute columns one row at a time in shared memory
  for (int k=0; k<N; k++) {
    rowk[tid] = A[k*stride+tid];
    __syncthreads();
    colk[ipiv[tid]] = rowk[tid];
    __syncthreads();
    A[k*stride+tid] = colk[tid];
    __syncthreads();
  }
}



template<typename T, int BS>
__global__ void
complex_inverse_many_naive_pivot (T **A_list, T **work_list, int N, int stride)
{
  int tid = threadIdx.x;
  __shared__ T *A;
  unsigned str = 2*stride;
  if (tid == 0) 
    A    = A_list[blockIdx.x];
  __syncthreads();
  __shared__ int kbar, ipiv[BS];  
  __shared__ T colk[2*BS], rowk[2*BS];
  __shared__ short imax[BS];
  
  ipiv[tid] = tid;
  __syncthreads();
  for (int k=0; k<N; k++) {
    // Find location of largest element in the column at or below k
    T re = A[tid*str+2*k];
    T im = A[tid*str+2*k+1];
    colk[tid] = (tid < k) ? 0.0 : re*re + im*im;
    __syncthreads();
    int skip = 1<<((int)ceil(log2((double)BS)-1.0e-6)-1);
    imax[tid] = tid;
    __syncthreads();
    for (; skip>0; skip>>=1) {
      if (tid < skip && (tid+skip)<N) 
    	// colk[tid] = max(colk[tid],colk[tid+skip]);
	if (colk[tid+skip] > colk[tid]) {
	  imax[tid] = imax[tid+skip];
	  colk[tid] = colk[tid+skip];
	}
      __syncthreads();
    }
    if (tid == 0) {
      kbar = imax[0];
      int i = ipiv[kbar];
      ipiv[kbar] = ipiv[k];
      ipiv[k] = i;
    }
    __syncthreads();
    // Swap rows
    rowk[tid   ] = A[kbar*str + tid];
    rowk[tid+BS] = A[kbar*str + tid + BS];
    colk[tid   ] = A[k   *str + tid];
    colk[tid+BS] = A[k   *str + tid + BS];
    __syncthreads();
    A[k   *str + tid   ] = rowk[tid];
    A[k   *str + tid+BS] = rowk[tid+BS];
    A[kbar*str + tid   ] = colk[tid];
    A[kbar*str + tid+BS] = colk[tid+BS];
    __syncthreads();
    re = rowk[2*k+0];
    im = rowk[2*k+1];
    T nrmInv = 1.0/(re*re + im*im);
    T pivInv_re =  nrmInv * re;
    T pivInv_im = -nrmInv * im;
    __syncthreads();
    if (tid == k) {
      rowk[2*k+0] = T();
      rowk[2*k+1] = T();
    }
    // Column scaling
    re = A[tid*str+2*k+0];
    im = A[tid*str+2*k+1];
    colk[2*tid+0] = -(pivInv_re*re - pivInv_im*im);
    colk[2*tid+1] = -(pivInv_im*re + pivInv_im*re);
    
    if (tid == k) 
      colk[2*tid+0] = colk[2*tid+1] = T();
    __syncthreads();
    A[tid*str+k   ] = colk[tid];
    A[tid*str+k+BS] = colk[tid+BS];
    __syncthreads();
    // Rank-1 update
    for (int j=0; j<N; j++) {
      re = colk[2*j+0]*rowk[2*tid+0] - colk[2*j+1]*rowk[2*tid+1];
      im = colk[2*j+0]*rowk[2*tid+1] + colk[2*j+1]*rowk[2*tid+0];
      A[j*str+2*tid  ] += re;
      A[j*str+2*tid+1] += im;
    }
    __syncthreads();
    // Row scaling
    re = pivInv_re * rowk[2*tid+0] - pivInv_im * rowk[2*tid+1];
    im = pivInv_re * rowk[2*tid+1] + pivInv_im * rowk[2*tid+0];
    if (tid == k) 
      re = im = T();
    A[k*str + 2*tid+0] = re;
    A[k*str + 2*tid+1] = im;
    __syncthreads();
  }
  // Now, permute columns one row at a time in shared memory
  for (int k=0; k<N; k++) {
    rowk[tid   ] = A[k*str+tid   ];
    rowk[tid+BS] = A[k*str+tid+BS];
    __syncthreads();
    colk[2*ipiv[tid]+0] = rowk[2*tid+0];
    colk[2*ipiv[tid]+1] = rowk[2*tid+1];
    __syncthreads();
    A[k*str+tid   ] = colk[tid   ];
    A[k*str+tid+BS] = colk[tid+BS];
    __syncthreads();
  }
}







#define CONVERT_BS 256


template<typename Tdest, typename Tsrc>
__global__ void
convert (Tdest **dest_list, Tsrc **src_list, int len)
{
  __shared__ Tsrc *mysrc;
  __shared__ Tdest *mydest;
  if (threadIdx.x ==0) {
    mysrc = src_list[blockIdx.y];
    mydest = dest_list[blockIdx.y];
  }
  __syncthreads();
  int i = blockIdx.x * CONVERT_BS + threadIdx.x;
  if (i < len)
    mydest[i] = (Tdest)mysrc[i];

}


template<typename Tdest, typename Tsrc>
__global__ void
convert (Tdest **dest_list, Tsrc **src_list, 
	 int dest_rows, int dest_cols, int dest_rowstride,
	 int src_rows,  int src_cols,  int src_rowstride)
{
  __shared__ Tsrc *mysrc;
  __shared__ Tdest *mydest;
  if (threadIdx.x ==0) {
    mysrc = src_list[blockIdx.y];
    mydest = dest_list[blockIdx.y];
  }
  __syncthreads();
  int i = blockIdx.x * CONVERT_BS + threadIdx.x;
  int row = i / dest_rowstride;
  int col = i - row*dest_rowstride;
  if (row < dest_rows && col < dest_rowstride) {
    if (col < src_cols && row < src_rows)
      mydest[i] = (Tdest)mysrc[row*src_rowstride + col];
    else
      mydest[i] = (row == col) ? (Tdest)1.0 : (Tdest)0.0;
  }
}





#define INVERSE_BS 16

void
cuda_inverse_many (float **Alist_d, float **worklist_d,
		   int N, int num_mats)
{
  dim3 dimBlock(INVERSE_BS,2);
  dim3 dimGrid(num_mats);
  
  inverse_many_pivot<float,INVERSE_BS><<<dimGrid,dimBlock>>> 
    (Alist_d, worklist_d, N, N);
}


size_t
cuda_inverse_many_worksize(int N)
{
  return (N * N + INVERSE_BS*INVERSE_BS);
}


size_t
cuda_inverse_many_double_worksize(int N)
{
  int N_double = ((N+INVERSE_BS-1)/INVERSE_BS) * INVERSE_BS;
  
  return 2*(2*N_double*N_double + INVERSE_BS*INVERSE_BS);
}

void
cuda_inverse_many_double (float *Alist_d[], float *worklist_d[],
			  int N, int num_mats)
{
  dim3 dimBlockConvert (CONVERT_BS);
  dim3 dimGridConvert(N*N/CONVERT_BS, num_mats);
  if (N*N % CONVERT_BS)
    dimGridConvert.x++;
  convert<<<dimGridConvert,dimBlockConvert>>> 
    ((double**)worklist_d, Alist_d, N*N);

  float *Alist_new[num_mats], *Alist_h[num_mats];
  float *worklist_h[num_mats];
  double *worklist_double_h[num_mats];

  cudaMemcpy (worklist_h, worklist_d, num_mats*sizeof(float*),
	      cudaMemcpyDeviceToHost);
  cudaMemcpy (Alist_h, Alist_d, num_mats*sizeof(float*),
	      cudaMemcpyDeviceToHost);

  for (int i=0; i<num_mats; i++) {
    Alist_new[i] = worklist_h[i];
    worklist_double_h[i] = (double*)(worklist_h[i]) +N*N;
  }
  cudaMemcpy (worklist_d, worklist_double_h, num_mats*sizeof(double*),
	      cudaMemcpyHostToDevice);
  cudaMemcpy (Alist_d, Alist_new, num_mats*sizeof(double*),
	      cudaMemcpyHostToDevice);


  dim3 dimBlock(INVERSE_BS,2);
  dim3 dimGrid(num_mats);
  
  inverse_many_pivot<double,INVERSE_BS><<<dimGrid,dimBlock>>> 
    ((double**)Alist_d, (double**)worklist_d, N, N);

  cudaThreadSynchronize();
  cudaError_t err = cudaGetLastError();
  if (err != cudaSuccess) {
    fprintf (stderr, "CUDA error in inverse_many_pivot<double,%d>:\n  %s\n",
	     INVERSE_BS, cudaGetErrorString(err));
    abort();
  }

  cudaMemcpy (Alist_d, Alist_h, num_mats*sizeof(float*),
	      cudaMemcpyHostToDevice);

  cudaMemcpy (worklist_d, worklist_h, num_mats*sizeof(float*),
	      cudaMemcpyHostToDevice);

  convert<<<dimGridConvert,dimBlockConvert>>> 
    (Alist_d, (double**) worklist_d, N*N);

}



// VERY slow implementation used only for debugging.
template<typename T1, typename T2, int BS>
__global__ void 
check_inv (T1 **A, T2 **B, int N, int Astride, int Bstride)
{
  int tid = threadIdx.x;
  __shared__ T1 *myA;
  __shared__ T2 *myB;
  if (tid == 0) {
    myA = A[blockIdx.x];
    myB = B[blockIdx.x];
  }
  __syncthreads();
  __shared__ double AB[BS];
  int NB = (N+BS-1)/BS;
  bool error = false;
  for (int row=0; row<N; row++)
    for (int col=0; col<N; col++) {
      AB[tid] = 0.0;
       int off = tid;
       for (int block=0; block<NB; block++) {
   	 if (off < N)
   	   AB[tid] += myA[row*Astride+off]*myB[off*Bstride+col];
   	off += BS;
       }
      // Now do reduction
      int skip = 1<<((int)ceil(log2((double)BS)-1.0e-6)-1);
      for (; skip>0; skip>>=1) {
  	if (tid < skip)
  	  AB[tid] += AB[tid+skip];
  	__syncthreads();
      }
    
      double expected = (row==col) ? 1.0 : 0.0;
      error = error || (fabs(AB[tid]-expected) > 1.0e-6);
    }
  if (tid == 0) 
    A[blockIdx.x] = (T1*)(error ? 1 : 0);
  __syncthreads();
}




void
cuda_inverse_many_double (float *Alist_d[], float *worklist_d[],
			  int N, int row_stride, int num_mats)
{
  int N_double = ((N + INVERSE_BS-1)/INVERSE_BS)*INVERSE_BS;
  
  dim3 dimBlockConvert (CONVERT_BS);
  dim3 dimGridConvert((N_double*N_double+(CONVERT_BS-1))/CONVERT_BS, 
		      num_mats);
  convert<<<dimGridConvert,dimBlockConvert>>> 
    ((double**)worklist_d, Alist_d, 
     N_double, N_double, N_double,
     N, N, row_stride);

  // We need to replace the pointers to the single-precision A
  // with the double-precision version we just converted.  We
  // Also need generate a new set of workspace pointers to point
  // to the region after the double-precision A.

  float *Alist_new[num_mats], *Alist_h[num_mats];
  float *worklist_h[num_mats];
  double *worklist_double_h[num_mats];
  float *bad_inverse[num_mats];

  // Save the original pointer lists on the host
  cudaMemcpy (worklist_h, worklist_d, num_mats*sizeof(float*),
	      cudaMemcpyDeviceToHost);
  cudaMemcpy (Alist_h, Alist_d, num_mats*sizeof(float*),
	      cudaMemcpyDeviceToHost);

  // Create new pointers as discussed above
  for (int i=0; i<num_mats; i++) {
    Alist_new[i] = worklist_h[i];
    worklist_double_h[i] = (double*)(worklist_h[i]) +N_double*N_double;
  }
  cudaMemcpy (worklist_d, worklist_double_h, num_mats*sizeof(double*),
	      cudaMemcpyHostToDevice);
  cudaMemcpy (Alist_d, Alist_new, num_mats*sizeof(double*),
	      cudaMemcpyHostToDevice);


  // Do the inversion in double-precision
  dim3 dimGrid(num_mats);
  
  // This appears to cause NANs for certain matrix sizes on occasion.
  // Check the pivoting algorithm.
//   dim3 dimBlock(INVERSE_BS,2);
//   inverse_many<double,INVERSE_BS><<<dimGrid,dimBlock>>> 
//     ((double**)Alist_d, (double**)worklist_d, N_double, N_double);
//   dim3 dimBlock(INVERSE_BS);
//   inverse_many_pivot<double,INVERSE_BS><<<dimGrid,dimBlock>>> 
//     ((double**)Alist_d, (double**)worklist_d, N_double, N_double);


  int NB = (N+15)/16;
  int BS=0;
  dim3 dimBlock(NB*16);
  switch (NB) {
    case 1:
      inverse_many_naive_pivot<double,16><<<dimGrid,dimBlock>>> 
	((double**)Alist_d, (double**)worklist_d, N, N_double);
      BS=16;
      break;
    case 2:
      inverse_many_naive_pivot<double,32><<<dimGrid,dimBlock>>> 
	((double**)Alist_d, (double**)worklist_d, N, N_double);
      BS=32;
      break;
    case 3:
      inverse_many_naive_pivot<double,48><<<dimGrid,dimBlock>>> 
	((double**)Alist_d, (double**)worklist_d, N, N_double);
      BS=48;
      break;
    case 4:
      inverse_many_naive_pivot<double,64><<<dimGrid,dimBlock>>> 
	((double**)Alist_d, (double**)worklist_d, N, N_double);
      BS=64;
      break;
    case 5:
      inverse_many_naive_pivot<double,80><<<dimGrid,dimBlock>>> 
	((double**)Alist_d, (double**)worklist_d, N, N_double);
      BS=80;
      break;
    case 6:
      inverse_many_naive_pivot<double,96><<<dimGrid,dimBlock>>> 
	((double**)Alist_d, (double**)worklist_d, N, N_double);
      BS=96;
      break;
    case 7:
      inverse_many_naive_pivot<double,112><<<dimGrid,dimBlock>>> 
	((double**)Alist_d, (double**)worklist_d, N, N_double);
      BS=112;
      break;
    case 8:
      inverse_many_naive_pivot<double,128><<<dimGrid,dimBlock>>> 
	((double**)Alist_d, (double**)worklist_d, N, N_double);
      BS=128;
      break;
    case 9:
      inverse_many_naive_pivot<double,144><<<dimGrid,dimBlock>>> 
	((double**)Alist_d, (double**)worklist_d, N, N_double);
      BS=144;
      break;
    case 10:
      inverse_many_naive_pivot<double,160><<<dimGrid,dimBlock>>> 
	((double**)Alist_d, (double**)worklist_d, N, N_double);
      BS=160;
      break;
    case 11:
      inverse_many_naive_pivot<double,176><<<dimGrid,dimBlock>>> 
	((double**)Alist_d, (double**)worklist_d, N, N_double);
      BS=176;
      break;
    case 12:
      inverse_many_naive_pivot<double,192><<<dimGrid,dimBlock>>> 
	((double**)Alist_d, (double**)worklist_d, N, N_double);
      BS=192;
      break;
    case 13:
      inverse_many_naive_pivot<double,208><<<dimGrid,dimBlock>>> 
	((double**)Alist_d, (double**)worklist_d, N, N_double);
      BS = 208;
      break;
    case 14:
      inverse_many_naive_pivot<double,224><<<dimGrid,dimBlock>>> 
	((double**)Alist_d, (double**)worklist_d, N, N_double);
      BS=224;
      break;
    case 15:
      inverse_many_naive_pivot<double,240><<<dimGrid,dimBlock>>> 
	((double**)Alist_d, (double**)worklist_d, N, N_double);
      BS=240;
      break;
    case 16:
      inverse_many_naive_pivot<double,256><<<dimGrid,dimBlock>>> 
	((double**)Alist_d, (double**)worklist_d, N, N_double);
      BS=256;
      break;
    case 17:
      inverse_many_naive_pivot<double,272><<<dimGrid,dimBlock>>> 
	((double**)Alist_d, (double**)worklist_d, N, N_double);
      BS=272;
      break;
    case 18:
      inverse_many_naive_pivot<double,288><<<dimGrid,dimBlock>>> 
	((double**)Alist_d, (double**)worklist_d, N, N_double);
      BS=288;
      break;
    case 19:
      inverse_many_naive_pivot<double,304><<<dimGrid,dimBlock>>> 
	((double**)Alist_d, (double**)worklist_d, N, N_double);
      BS=304;
      break;
    case 20:
      inverse_many_naive_pivot<double,320><<<dimGrid,dimBlock>>> 
	((double**)Alist_d, (double**)worklist_d, N, N_double);
      BS=320;
      break;
    case 21:
      inverse_many_naive_pivot<double,336><<<dimGrid,dimBlock>>> 
	((double**)Alist_d, (double**)worklist_d, N, N_double);
      BS=336;
      break;
    case 22:
      inverse_many_naive_pivot<double,352><<<dimGrid,dimBlock>>>
	((double**)Alist_d, (double**)worklist_d, N, N_double);
      BS = 352;
      break;
    case 23:
      inverse_many_naive_pivot<double,368><<<dimGrid,dimBlock>>>
	((double**)Alist_d, (double**)worklist_d, N, N_double);
      BS = 368;
      break;
    case 24:
      inverse_many_naive_pivot<double,384><<<dimGrid,dimBlock>>>
	((double**)Alist_d, (double**)worklist_d, N, N_double);
      BS = 384;
      break;
    case 25:
      inverse_many_naive_pivot<double,400><<<dimGrid,dimBlock>>>
	((double**)Alist_d, (double**)worklist_d, N, N_double);
      BS = 400;
      break;
    case 26:
      inverse_many_naive_pivot<double,416><<<dimGrid,dimBlock>>>
	((double**)Alist_d, (double**)worklist_d, N, N_double);
      BS = 416;
      break;
    case 27:
      inverse_many_naive_pivot<double,432><<<dimGrid,dimBlock>>>
	((double**)Alist_d, (double**)worklist_d, N, N_double);
      BS = 432;
      break;
    case 28:
      inverse_many_naive_pivot<double,448><<<dimGrid,dimBlock>>>
	((double**)Alist_d, (double**)worklist_d, N, N_double);
      BS = 448;
      break;
    case 29:
      inverse_many_naive_pivot<double,464><<<dimGrid,dimBlock>>>
	((double**)Alist_d, (double**)worklist_d, N, N_double);
      BS = 464;
      break;
    case 30:
      inverse_many_naive_pivot<double,480><<<dimGrid,dimBlock>>>
	((double**)Alist_d, (double**)worklist_d, N, N_double);
      BS = 480;
      break;
    case 31:
      inverse_many_naive_pivot<double,496><<<dimGrid,dimBlock>>>
	((double**)Alist_d, (double**)worklist_d, N, N_double);
      BS = 496;
      break;
    case 32:
      inverse_many_naive_pivot<double,512><<<dimGrid,dimBlock>>>
	((double**)Alist_d, (double**)worklist_d, N, N_double);
      BS = 512;
      break;
    case 33:
      inverse_many_naive_pivot<double,528><<<dimGrid,dimBlock>>>
	((double**)Alist_d, (double**)worklist_d, N, N_double);
      BS = 528;
      break;
    case 34:
      inverse_many_naive_pivot<double,544><<<dimGrid,dimBlock>>>
	((double**)Alist_d, (double**)worklist_d, N, N_double);
      BS = 544;
      break;
    case 35:
      inverse_many_naive_pivot<double,560><<<dimGrid,dimBlock>>>
	((double**)Alist_d, (double**)worklist_d, N, N_double);
      BS = 560;
      break;
    case 36:
      inverse_many_naive_pivot<double,576><<<dimGrid,dimBlock>>>
	((double**)Alist_d, (double**)worklist_d, N, N_double);
      BS = 576;
      break;
    case 37:
      inverse_many_naive_pivot<double,592><<<dimGrid,dimBlock>>>
	((double**)Alist_d, (double**)worklist_d, N, N_double);
      BS = 592;
      break;
    case 38:
      inverse_many_naive_pivot<double,608><<<dimGrid,dimBlock>>>
	((double**)Alist_d, (double**)worklist_d, N, N_double);
      BS = 608;
    break;
    case 39:
      inverse_many_naive_pivot<double,624><<<dimGrid,dimBlock>>>
	((double**)Alist_d, (double**)worklist_d, N, N_double);
      BS = 624;
      break;
    case 40:
      inverse_many_naive_pivot<double,640><<<dimGrid,dimBlock>>>
	((double**)Alist_d, (double**)worklist_d, N, N_double);
      BS = 640;
      break;
    case 41:
      inverse_many_naive_pivot<double,656><<<dimGrid,dimBlock>>>
	((double**)Alist_d, (double**)worklist_d, N, N_double);
      BS = 656;
      break;
    case 42:
      inverse_many_naive_pivot<double,672><<<dimGrid,dimBlock>>>
	((double**)Alist_d, (double**)worklist_d, N, N_double);
      BS = 672;
      break;
    case 43:
      inverse_many_naive_pivot<double,688><<<dimGrid,dimBlock>>>
	((double**)Alist_d, (double**)worklist_d, N, N_double);
      BS = 688;
      break;
    case 44:
      inverse_many_naive_pivot<double,704><<<dimGrid,dimBlock>>>
	((double**)Alist_d, (double**)worklist_d, N, N_double);
      BS = 704;
      break;
    case 45:
      inverse_many_naive_pivot<double,720><<<dimGrid,dimBlock>>>
	((double**)Alist_d, (double**)worklist_d, N, N_double);
      BS = 720;
      break;
    case 46:
      inverse_many_naive_pivot<double,736><<<dimGrid,dimBlock>>>
	((double**)Alist_d, (double**)worklist_d, N, N_double);
      BS = 736;
      break;
    // case 47:
    //   inverse_many_naive_pivot<double,752><<<dimGrid,dimBlock>>>
    // 	((double**)Alist_d, (double**)worklist_d, N, N_double);
    //   BS = 752;
    //   break;
    // case 48:
    //   inverse_many_naive_pivot<double,768><<<dimGrid,dimBlock>>>
    // 	((double**)Alist_d, (double**)worklist_d, N, N_double);
    //   BS = 768;
    //   break;
    // case 49:
    //   inverse_many_naive_pivot<double,784><<<dimGrid,dimBlock>>>
    // 	((double**)Alist_d, (double**)worklist_d, N, N_double);
    //   BS = 784;
    //   break;
    // case 50:
    //   inverse_many_naive_pivot<double,800><<<dimGrid,dimBlock>>>
    // 	((double**)Alist_d, (double**)worklist_d, N, N_double);
    //   BS = 800;
    //   break;

    default:
      fprintf (stderr, "N=%d is larger than maximum 736 in cuda_inverse_many_double.\n");
  };

  cudaThreadSynchronize();
  cudaError_t err = cudaGetLastError();
  if (err != cudaSuccess) {
    char buff[500];
    gethostname(buff, 500);
    fprintf (stderr, "CUDA error in inverse_many<double,%d>:\n  %s\n",
	     BS, cudaGetErrorString(err));
    fprintf (stderr, "Hostname = %s\n", buff);

    abort();
  }


#ifdef CHECK_INVERSES
  // Check inverses for correctness
  // Copy A matrix pointers to worklist_d
  cudaMemcpy (worklist_d, Alist_h, num_mats*sizeof(double*),
  	      cudaMemcpyHostToDevice);
  // Call kernel to check (A^(-1) * A - I).  It returns 0 or 1 in
  // the Alist_d pointer.
  dim3 checkBlock(32);
  dim3 checkGrid(num_mats);
  check_inv<double,float,32><<<checkGrid,checkBlock>>>((double**)Alist_d, worklist_d, N, N_double, row_stride);
  cudaThreadSynchronize();
  cudaError_t checkerr = cudaGetLastError();
  if (checkerr != cudaSuccess) {
    char buff[500];
    gethostname(buff, 500);
    fprintf (stderr, "CUDA error in check_inv<double,float,32>:\n  %s\n",
	     cudaGetErrorString(checkerr));
    fprintf (stderr, "Hostname = %s\n", buff);

    abort();
  }

  cudaMemcpy (bad_inverse, Alist_d, num_mats*sizeof(double*),
	      cudaMemcpyDeviceToHost);
  for (int mat=0; mat<num_mats; mat++)
    if (bad_inverse[mat]) {
      char name[1000];
      gethostname(name, 1000);
      fprintf (stderr, "Offending hostname = %s\n", name);
      std::cerr << "bad inverse for matrix " << mat << std::endl;
      std::vector<float>  Amat(N*row_stride);
      cudaMemcpy (&(Amat[0]), Alist_h[mat], N*row_stride*sizeof(float), cudaMemcpyDeviceToHost);
      std::ostringstream matName;
      matName << "BadMat_" << mat << ".dat";
      FILE *fout = fopen (matName.str().c_str(), "w");
      for (int row=0; row<N; row++) {
	for (int col=0; col<N; col++) {
	  //fprintf (stderr, "row=%d col=%d\n", row, col);
	  fprintf (fout, "%24.16e ", Amat[row*row_stride+col]);
	}
	fprintf (fout, "\n");
      }
      fclose (fout);

      std::vector<double> Ainv (N*N_double);
      cudaMemcpy (&(Ainv[0]), worklist_h[mat], N*N_double*sizeof(double), cudaMemcpyDeviceToHost);
      std::ostringstream invName;
      invName << "BadInv_" << mat << ".dat";
      fout = fopen (invName.str().c_str(), "w");
      for (int row=0; row<N; row++) {
	for (int col=0; col<N; col++) {
	  //fprintf (stderr, "row=%d col=%d\n", row, col);
	  fprintf (fout, "%24.16e ", Ainv[row*N_double+col]);
	}
	fprintf (fout, "\n");
      }
      fclose (fout);
    }

#endif

  // Copy original pointer lists back to device
  cudaMemcpy (Alist_d, Alist_h, num_mats*sizeof(float*),
	      cudaMemcpyHostToDevice);

  cudaMemcpy (worklist_d, worklist_h, num_mats*sizeof(float*),
	      cudaMemcpyHostToDevice);

  dim3 dimGridConvert2((N*row_stride+(CONVERT_BS-1))/CONVERT_BS, num_mats);

  // Convert back to single precision.
  convert<<<dimGridConvert2,dimBlockConvert>>> 
    (Alist_d, (double**) worklist_d, 
     N, N, row_stride,
     N_double, N_double, N_double);
}




void
cuda_inverse_many_complex_double (float *Alist_d[], float *worklist_d[],
				  int N, int row_stride, int num_mats)
{
  int N_double = ((N + INVERSE_BS-1)/INVERSE_BS)*INVERSE_BS;
  
  dim3 dimBlockConvert (CONVERT_BS);
  dim3 dimGridConvert((N_double*N_double+(CONVERT_BS-1))/CONVERT_BS, 
		      num_mats);
  convert<<<dimGridConvert,dimBlockConvert>>> 
    ((double**)worklist_d, Alist_d, 
     N_double, N_double, N_double,
     N, N, row_stride);

  // We need to replace the pointers to the single-precision A
  // with the double-precision version we just converted.  We
  // Also need generate a new set of workspace pointers to point
  // to the region after the double-precision A.

  float *Alist_new[num_mats], *Alist_h[num_mats];
  float *worklist_h[num_mats];
  double *worklist_double_h[num_mats];
  float *bad_inverse[num_mats];

  // Save the original pointer lists on the host
  cudaMemcpy (worklist_h, worklist_d, num_mats*sizeof(float*),
	      cudaMemcpyDeviceToHost);
  cudaMemcpy (Alist_h, Alist_d, num_mats*sizeof(float*),
	      cudaMemcpyDeviceToHost);

  // Create new pointers as discussed above
  for (int i=0; i<num_mats; i++) {
    Alist_new[i] = worklist_h[i];
    worklist_double_h[i] = (double*)(worklist_h[i]) +N_double*N_double;
  }
  cudaMemcpy (worklist_d, worklist_double_h, num_mats*sizeof(double*),
	      cudaMemcpyHostToDevice);
  cudaMemcpy (Alist_d, Alist_new, num_mats*sizeof(double*),
	      cudaMemcpyHostToDevice);


  // Do the inversion in double-precision
  dim3 dimGrid(num_mats);
  
  int NB = (N+15)/16;
  int BS=0;
  dim3 dimBlock(NB*16);
  switch (NB) {
    case 1:
      complex_inverse_many_naive_pivot<double,16><<<dimGrid,dimBlock>>> 
	((double**)Alist_d, (double**)worklist_d, N, N_double);
      BS=16;
      break;
    case 2:
      complex_inverse_many_naive_pivot<double,32><<<dimGrid,dimBlock>>> 
	((double**)Alist_d, (double**)worklist_d, N, N_double);
      BS=32;
      break;
    case 3:
      complex_inverse_many_naive_pivot<double,48><<<dimGrid,dimBlock>>> 
	((double**)Alist_d, (double**)worklist_d, N, N_double);
      BS=48;
      break;
    case 4:
      complex_inverse_many_naive_pivot<double,64><<<dimGrid,dimBlock>>> 
	((double**)Alist_d, (double**)worklist_d, N, N_double);
      BS=64;
      break;
    case 5:
      complex_inverse_many_naive_pivot<double,80><<<dimGrid,dimBlock>>> 
	((double**)Alist_d, (double**)worklist_d, N, N_double);
      BS=80;
      break;
    case 6:
      complex_inverse_many_naive_pivot<double,96><<<dimGrid,dimBlock>>> 
	((double**)Alist_d, (double**)worklist_d, N, N_double);
      BS=96;
      break;
    case 7:
      complex_inverse_many_naive_pivot<double,112><<<dimGrid,dimBlock>>> 
	((double**)Alist_d, (double**)worklist_d, N, N_double);
      BS=112;
      break;
    case 8:
      complex_inverse_many_naive_pivot<double,128><<<dimGrid,dimBlock>>> 
	((double**)Alist_d, (double**)worklist_d, N, N_double);
      BS=128;
      break;
    case 9:
      complex_inverse_many_naive_pivot<double,144><<<dimGrid,dimBlock>>> 
	((double**)Alist_d, (double**)worklist_d, N, N_double);
      BS=144;
      break;
    case 10:
      complex_inverse_many_naive_pivot<double,160><<<dimGrid,dimBlock>>> 
	((double**)Alist_d, (double**)worklist_d, N, N_double);
      BS=160;
      break;
    case 11:
      complex_inverse_many_naive_pivot<double,176><<<dimGrid,dimBlock>>> 
	((double**)Alist_d, (double**)worklist_d, N, N_double);
      BS=176;
      break;
    case 12:
      complex_inverse_many_naive_pivot<double,192><<<dimGrid,dimBlock>>> 
	((double**)Alist_d, (double**)worklist_d, N, N_double);
      BS=192;
      break;
    case 13:
      complex_inverse_many_naive_pivot<double,208><<<dimGrid,dimBlock>>> 
	((double**)Alist_d, (double**)worklist_d, N, N_double);
      BS = 208;
      break;
    case 14:
      complex_inverse_many_naive_pivot<double,224><<<dimGrid,dimBlock>>> 
	((double**)Alist_d, (double**)worklist_d, N, N_double);
      BS=224;
      break;
    case 15:
      complex_inverse_many_naive_pivot<double,240><<<dimGrid,dimBlock>>> 
	((double**)Alist_d, (double**)worklist_d, N, N_double);
      BS=240;
      break;
    case 16:
      complex_inverse_many_naive_pivot<double,256><<<dimGrid,dimBlock>>> 
	((double**)Alist_d, (double**)worklist_d, N, N_double);
      BS=256;
      break;
    case 17:
      complex_inverse_many_naive_pivot<double,272><<<dimGrid,dimBlock>>> 
	((double**)Alist_d, (double**)worklist_d, N, N_double);
      BS=272;
      break;
    case 18:
      complex_inverse_many_naive_pivot<double,288><<<dimGrid,dimBlock>>> 
	((double**)Alist_d, (double**)worklist_d, N, N_double);
      BS=288;
      break;
    case 19:
      complex_inverse_many_naive_pivot<double,304><<<dimGrid,dimBlock>>> 
	((double**)Alist_d, (double**)worklist_d, N, N_double);
      BS=304;
      break;
    case 20:
      complex_inverse_many_naive_pivot<double,320><<<dimGrid,dimBlock>>> 
	((double**)Alist_d, (double**)worklist_d, N, N_double);
      BS=320;
      break;
    case 21:
      complex_inverse_many_naive_pivot<double,336><<<dimGrid,dimBlock>>> 
	((double**)Alist_d, (double**)worklist_d, N, N_double);
      BS=336;
      break;
    case 22:
      complex_inverse_many_naive_pivot<double,352><<<dimGrid,dimBlock>>>
	((double**)Alist_d, (double**)worklist_d, N, N_double);
      BS = 352;
      break;
    case 23:
      complex_inverse_many_naive_pivot<double,368><<<dimGrid,dimBlock>>>
	((double**)Alist_d, (double**)worklist_d, N, N_double);
      BS = 368;
      break;
    case 24:
      complex_inverse_many_naive_pivot<double,384><<<dimGrid,dimBlock>>>
	((double**)Alist_d, (double**)worklist_d, N, N_double);
      BS = 384;
      break;
    case 25:
      complex_inverse_many_naive_pivot<double,400><<<dimGrid,dimBlock>>>
	((double**)Alist_d, (double**)worklist_d, N, N_double);
      BS = 400;
      break;
    case 26:
      complex_inverse_many_naive_pivot<double,416><<<dimGrid,dimBlock>>>
	((double**)Alist_d, (double**)worklist_d, N, N_double);
      BS = 416;
      break;
    // case 27:
    //   complex_inverse_many_naive_pivot<double,432><<<dimGrid,dimBlock>>>
    // 	((double**)Alist_d, (double**)worklist_d, N, N_double);
    //   BS = 432;
    //   break;
    // case 28:
    //   complex_inverse_many_naive_pivot<double,448><<<dimGrid,dimBlock>>>
    // 	((double**)Alist_d, (double**)worklist_d, N, N_double);
    //   BS = 448;
    //   break;
    // case 29:
    //   complex_inverse_many_naive_pivot<double,464><<<dimGrid,dimBlock>>>
    // 	((double**)Alist_d, (double**)worklist_d, N, N_double);
    //   BS = 464;
    //   break;
    // case 30:
    //   complex_inverse_many_naive_pivot<double,480><<<dimGrid,dimBlock>>>
    // 	((double**)Alist_d, (double**)worklist_d, N, N_double);
    //   BS = 480;
    //   break;
    // case 31:
    //   complex_inverse_many_naive_pivot<double,496><<<dimGrid,dimBlock>>>
    // 	((double**)Alist_d, (double**)worklist_d, N, N_double);
    //   BS = 496;
    //   break;
    // case 32:
    //   complex_inverse_many_naive_pivot<double,512><<<dimGrid,dimBlock>>>
    // 	((double**)Alist_d, (double**)worklist_d, N, N_double);
    //   BS = 512;
    //   break;
    // case 33:
    //   complex_inverse_many_naive_pivot<double,528><<<dimGrid,dimBlock>>>
    // 	((double**)Alist_d, (double**)worklist_d, N, N_double);
    //   BS = 528;
    //   break;
    // case 34:
    //   complex_inverse_many_naive_pivot<double,544><<<dimGrid,dimBlock>>>
    // 	((double**)Alist_d, (double**)worklist_d, N, N_double);
    //   BS = 544;
    //   break;
    // case 35:
    //   complex_inverse_many_naive_pivot<double,560><<<dimGrid,dimBlock>>>
    // 	((double**)Alist_d, (double**)worklist_d, N, N_double);
    //   BS = 560;
    //   break;
    // case 36:
    //   complex_inverse_many_naive_pivot<double,576><<<dimGrid,dimBlock>>>
    // 	((double**)Alist_d, (double**)worklist_d, N, N_double);
    //   BS = 576;
    //   break;
    // case 37:
    //   complex_inverse_many_naive_pivot<double,592><<<dimGrid,dimBlock>>>
    // 	((double**)Alist_d, (double**)worklist_d, N, N_double);
    //   BS = 592;
    //   break;
    // case 38:
    //   complex_inverse_many_naive_pivot<double,608><<<dimGrid,dimBlock>>>
    // 	((double**)Alist_d, (double**)worklist_d, N, N_double);
    //   BS = 608;
    // break;
    // case 39:
    //   complex_inverse_many_naive_pivot<double,624><<<dimGrid,dimBlock>>>
    // 	((double**)Alist_d, (double**)worklist_d, N, N_double);
    //   BS = 624;
    //   break;
    // case 40:
    //   complex_inverse_many_naive_pivot<double,640><<<dimGrid,dimBlock>>>
    // 	((double**)Alist_d, (double**)worklist_d, N, N_double);
    //   BS = 640;
    //   break;
    // case 41:
    //   complex_inverse_many_naive_pivot<double,656><<<dimGrid,dimBlock>>>
    // 	((double**)Alist_d, (double**)worklist_d, N, N_double);
    //   BS = 656;
    //   break;
    // case 42:
    //   complex_inverse_many_naive_pivot<double,672><<<dimGrid,dimBlock>>>
    // 	((double**)Alist_d, (double**)worklist_d, N, N_double);
    //   BS = 672;
    //   break;
    // case 43:
    //   complex_inverse_many_naive_pivot<double,688><<<dimGrid,dimBlock>>>
    // 	((double**)Alist_d, (double**)worklist_d, N, N_double);
    //   BS = 688;
    //   break;
    // case 44:
    //   complex_inverse_many_naive_pivot<double,704><<<dimGrid,dimBlock>>>
    // 	((double**)Alist_d, (double**)worklist_d, N, N_double);
    //   BS = 704;
    //   break;
    // case 45:
    //   complex_inverse_many_naive_pivot<double,720><<<dimGrid,dimBlock>>>
    // 	((double**)Alist_d, (double**)worklist_d, N, N_double);
    //   BS = 720;
    //   break;
    // case 46:
    //   complex_inverse_many_naive_pivot<double,736><<<dimGrid,dimBlock>>>
    // 	((double**)Alist_d, (double**)worklist_d, N, N_double);
    //   BS = 736;
    //   break;
    // case 47:
    //   complex_inverse_many_naive_pivot<double,752><<<dimGrid,dimBlock>>>
    // 	((double**)Alist_d, (double**)worklist_d, N, N_double);
    //   BS = 752;
    //   break;
    // case 48:
    //   complex_inverse_many_naive_pivot<double,768><<<dimGrid,dimBlock>>>
    // 	((double**)Alist_d, (double**)worklist_d, N, N_double);
    //   BS = 768;
    //   break;
    // case 49:
    //   complex_inverse_many_naive_pivot<double,784><<<dimGrid,dimBlock>>>
    // 	((double**)Alist_d, (double**)worklist_d, N, N_double);
    //   BS = 784;
    //   break;
    // case 50:
    //   complex_inverse_many_naive_pivot<double,800><<<dimGrid,dimBlock>>>
    // 	((double**)Alist_d, (double**)worklist_d, N, N_double);
    //   BS = 800;
    //   break;

    default:
      fprintf (stderr, "N=%d is larger than maximum 416 in cuda_complex_inverse_many_double.\n");
  };

  cudaThreadSynchronize();
  cudaError_t err = cudaGetLastError();
  if (err != cudaSuccess) {
    char buff[500];
    gethostname(buff, 500);
    fprintf (stderr, "CUDA error in complex_inverse_many<double,%d>:\n  %s\n",
	     BS, cudaGetErrorString(err));
    fprintf (stderr, "Hostname = %s\n", buff);

    abort();
  }

#ifdef CHECK_INVERSES
  // Check inverses for correctness
  // Copy A matrix pointers to worklist_d
  cudaMemcpy (worklist_d, Alist_h, num_mats*sizeof(double*),
  	      cudaMemcpyHostToDevice);
  // Call kernel to check (A^(-1) * A - I).  It returns 0 or 1 in
  // the Alist_d pointer.
  dim3 checkBlock(32);
  dim3 checkGrid(num_mats);
  check_inv<double,float,32><<<checkGrid,checkBlock>>>((double**)Alist_d, worklist_d, N, N_double, row_stride);
  cudaThreadSynchronize();
  cudaError_t checkerr = cudaGetLastError();
  if (checkerr != cudaSuccess) {
    char buff[500];
    gethostname(buff, 500);
    fprintf (stderr, "CUDA error in check_inv<double,float,32>:\n  %s\n",
	     cudaGetErrorString(checkerr));
    fprintf (stderr, "Hostname = %s\n", buff);

    abort();
  }

  cudaMemcpy (bad_inverse, Alist_d, num_mats*sizeof(double*),
	      cudaMemcpyDeviceToHost);
  for (int mat=0; mat<num_mats; mat++)
    if (bad_inverse[mat]) {
      char name[1000];
      gethostname(name, 1000);
      fprintf (stderr, "Offending hostname = %s\n", name);
      std::cerr << "bad inverse for matrix " << mat << std::endl;
      std::vector<float>  Amat(N*row_stride);
      cudaMemcpy (&(Amat[0]), Alist_h[mat], N*row_stride*sizeof(float), cudaMemcpyDeviceToHost);
      std::ostringstream matName;
      matName << "BadMat_" << mat << ".dat";
      FILE *fout = fopen (matName.str().c_str(), "w");
      for (int row=0; row<N; row++) {
	for (int col=0; col<N; col++) {
	  //fprintf (stderr, "row=%d col=%d\n", row, col);
	  fprintf (fout, "%24.16e ", Amat[row*row_stride+col]);
	}
	fprintf (fout, "\n");
      }
      fclose (fout);

      std::vector<double> Ainv (N*N_double);
      cudaMemcpy (&(Ainv[0]), worklist_h[mat], N*N_double*sizeof(double), cudaMemcpyDeviceToHost);
      std::ostringstream invName;
      invName << "BadInv_" << mat << ".dat";
      fout = fopen (invName.str().c_str(), "w");
      for (int row=0; row<N; row++) {
	for (int col=0; col<N; col++) {
	  //fprintf (stderr, "row=%d col=%d\n", row, col);
	  fprintf (fout, "%24.16e ", Ainv[row*N_double+col]);
	}
	fprintf (fout, "\n");
      }
      fclose (fout);
    }

#endif
  // Copy original pointer lists back to device
  cudaMemcpy (Alist_d, Alist_h, num_mats*sizeof(float*),
	      cudaMemcpyHostToDevice);

  cudaMemcpy (worklist_d, worklist_h, num_mats*sizeof(float*),
	      cudaMemcpyHostToDevice);

  dim3 dimGridConvert2((N*row_stride+(CONVERT_BS-1))/CONVERT_BS, num_mats);

  // Convert back to single precision.
  convert<<<dimGridConvert2,dimBlockConvert>>> 
    (Alist_d, (double**) worklist_d, 
     N, N, row_stride,
     N_double, N_double, N_double);
}



void
cuda_inverse_many_double (double *Alist_d[], double *worklist_d[],
			  int N, int num_mats)
{
  dim3 dimBlock(INVERSE_BS,2);
  dim3 dimGrid(num_mats);
  
  inverse_many_pivot<double,INVERSE_BS><<<dimGrid,dimBlock>>> 
    (Alist_d, worklist_d, N, N);
}




void
cuda_inverse_many (double *Alist_d[], double *worklist_d[],
		   int N, int num_mats)
{
  dim3 dimBlock(INVERSE_BS);
  dim3 dimGrid(num_mats);
  
  inverse_many_pivot<double,INVERSE_BS><<<dimGrid,dimBlock>>> 
    (Alist_d, worklist_d, N, N);
}



//////////////////////////////////////////////////////
//                  Test routines                   //
//////////////////////////////////////////////////////



#ifdef CUDA_TEST_MAIN

void 
test_inverse()
{
  int N = 32;
  dim3 dimBlock(INVERSE_BS,2);
  dim3 dimGrid(1);

  float *A_d, *work_d;
  int lwork = N*N + INVERSE_BS * INVERSE_BS;


  cudaMalloc((void**)&A_d, N*N*sizeof(float));
  cudaMalloc((void**)&work_d, lwork*sizeof(float));
  
  float A[N*N], Ainv[N*N];
  for (int i=0; i<N*N; i++)
    A[i] = drand48();
  cudaMemcpy(A_d, A, N*N*sizeof(float), cudaMemcpyHostToDevice);
  
  inverse<float,INVERSE_BS><<<dimGrid,dimBlock>>> (A_d, work_d, N, N);
  cudaThreadSynchronize();
  cudaError_t err = cudaGetLastError();
  if (err != cudaSuccess) {
    fprintf (stderr, "CUDA error in block_inverse:\n  %s\n",
	     cudaGetErrorString(err));
    abort();
  }

  // Copy Ainv back to host memory
  
  cudaMemcpy(Ainv, A_d, N*N*sizeof(float), cudaMemcpyDeviceToHost);

  float error = 0.0;
  for (int i=0; i<N; i++)
    for (int j=0; j<N; j++) {
      float val = 0.0;
      for (int k=0; k<N; k++)
	val += Ainv[i*N+k]*A[k*N+j];
      float diff = (i==j) ? (1.0f-val) : val;
      error += diff*diff;
    }
  fprintf (stderr, "error = %1.8e\n", sqrt(error/(double)(N*N)));

}



void 
test_inverse_many()
{
  int numMats = 10000;

  int N = 128;

  int lwork = N*N + INVERSE_BS * INVERSE_BS;
  fprintf (stderr, "lwork = %d\n", lwork);

  float **Alist, **worklist;
  float **Alist_d, **worklist_d;

  Alist    = (float**)malloc(numMats*sizeof(float*));
  worklist = (float**)malloc(numMats*sizeof(float*));
  cudaMalloc((void**)&Alist_d,    numMats*sizeof(float*));
  cudaMalloc((void**)&worklist_d, numMats*sizeof(float*));

  float A[N*N];
  for (int i=0; i<N*N; i++)
    A[i] = drand48();

  for (int mat=0; mat<numMats; mat++) {
    cudaMalloc ((void**)&(Alist[mat]),    N*N*sizeof(float));
    cudaMalloc ((void**)&(worklist[mat]), lwork*sizeof(float));
    cudaMemcpy(Alist[mat], A, N*N*sizeof(float), cudaMemcpyHostToDevice);
  }

  cudaMemcpy(Alist_d   ,    Alist, numMats*sizeof(float*), 
	     cudaMemcpyHostToDevice);
  cudaMemcpy(worklist_d, worklist, numMats*sizeof(float*), 
	     cudaMemcpyHostToDevice);
  
  dim3 dimBlock(INVERSE_BS,2);
  dim3 dimGrid(numMats);

  clock_t start = clock();
  for (int i=0; i<1; i++) {
    inverse_many_pivot<float,INVERSE_BS><<<dimGrid,dimBlock>>> 
      (Alist_d, worklist_d, N, N);
//     inverse_many<float,INVERSE_BS><<<dimGrid,dimBlock>>> 
//       (Alist_d, worklist_d, N, N);
    cudaThreadSynchronize();
  }
  clock_t end = clock();
  
  double time = (double)(end-start)/(double)CLOCKS_PER_SEC
    / (double)numMats;
  double rate = 1.0/time;
  fprintf (stderr, "Rate is %1.3f matrix inversions per second.\n",
	   rate);


  cudaError_t err = cudaGetLastError();
  if (err != cudaSuccess) {
    fprintf (stderr, "CUDA error in block_inverse:\n  %s\n",
	     cudaGetErrorString(err));
    abort();
  }

  // Copy Ainv back to host memory
  float Ainv[N*N];
  cudaMemcpy(Ainv, Alist[10], N*N*sizeof(float), cudaMemcpyDeviceToHost);

  float error = 0.0;
  for (int i=0; i<N; i++)
    for (int j=0; j<N; j++) {
      float val = 0.0;
      for (int k=0; k<N; k++)
	val += Ainv[i*N+k]*A[k*N+j];
      float diff = (i==j) ? (1.0f-val) : val;
      error += diff*diff;
    }
  fprintf (stderr, "error = %1.8e\n", sqrt(error/(double)(N*N)));

}


void 
test_inverse_many_double()
{
  int numMats = 1000;

  int N = 128;

  int lwork = N*N + INVERSE_BS * INVERSE_BS;
  fprintf (stderr, "lwork = %d\n", lwork);

  double **Alist, **worklist;
  double **Alist_d, **worklist_d;

  Alist    = (double**)malloc(numMats*sizeof(double*));
  worklist = (double**)malloc(numMats*sizeof(double*));
  cudaMalloc((void**)&Alist_d,    numMats*sizeof(double*));
  cudaMalloc((void**)&worklist_d, numMats*sizeof(double*));

  double A[N*N];
  for (int i=0; i<N*N; i++)
    A[i] = drand48();

  for (int mat=0; mat<numMats; mat++) {
    cudaMalloc ((void**)&(Alist[mat]),    N*N*sizeof(double));
    cudaMalloc ((void**)&(worklist[mat]), lwork*sizeof(double));
    cudaMemcpy(Alist[mat], A, N*N*sizeof(double), cudaMemcpyHostToDevice);
  }

  cudaMemcpy(Alist_d   ,    Alist, numMats*sizeof(double*), 
	     cudaMemcpyHostToDevice);
  cudaMemcpy(worklist_d, worklist, numMats*sizeof(double*), 
	     cudaMemcpyHostToDevice);
  
  dim3 dimBlock(INVERSE_BS,2);
  dim3 dimGrid(numMats);

  clock_t start = clock();
  for (int i=0; i<1; i++) {
    inverse_many_pivot<double,INVERSE_BS><<<dimGrid,dimBlock>>> 
      (Alist_d, worklist_d, N, N);
//     inverse_many<double,INVERSE_BS><<<dimGrid,dimBlock>>> 
//       (Alist_d, worklist_d, N, N);
    cudaThreadSynchronize();
  }
  clock_t end = clock();
  
  double time = (double)(end-start)/(double)CLOCKS_PER_SEC
    / (double)numMats;
  double rate = 1.0/time;
  fprintf (stderr, "Rate is %1.3f matrix inversions per second.\n",
	   rate);


  cudaError_t err = cudaGetLastError();
  if (err != cudaSuccess) {
    fprintf (stderr, "CUDA error in block_inverse:\n  %s\n",
	     cudaGetErrorString(err));
    abort();
  }

  // Copy Ainv back to host memory
  double Ainv[N*N];
  cudaMemcpy(Ainv, Alist[10], N*N*sizeof(double), cudaMemcpyDeviceToHost);

  double error = 0.0;
  for (int i=0; i<N; i++)
    for (int j=0; j<N; j++) {
      double val = 0.0;
      for (int k=0; k<N; k++)
	val += Ainv[i*N+k]*A[k*N+j];
      double diff = (i==j) ? (1.0f-val) : val;
      error += diff*diff;
    }
  fprintf (stderr, "error = %1.8e\n", sqrt(error/(double)(N*N)));

}




void 
test_inverse_many_double_conv()
{
  
  srand48((long) 12394);
  int numMats = 100;

  int N = 244;
  int row_stride = 256;

  int lwork = cuda_inverse_many_double_worksize(N);
  fprintf (stderr, "lwork = %d\n", lwork);

  float **Alist, **worklist;
  float **Alist_d, **worklist_d;

  Alist    = (float**)malloc(numMats*sizeof(float*));
  worklist = (float**)malloc(numMats*sizeof(float*));
  cudaMalloc((void**)&Alist_d,    numMats*sizeof(float*));
  cudaMalloc((void**)&worklist_d, numMats*sizeof(float*));

  float *A = (float*)malloc(sizeof(float)*numMats*N*row_stride);
  for (int j=0; j<numMats; j++)
    for (int i=0; i<N*row_stride; i++)
      A[j*N*row_stride+i] = 1.0*(drand48()-0.5);

  for (int mat=0; mat<numMats; mat++) {
    cudaMalloc ((void**)&(Alist[mat]),    N*row_stride*sizeof(float));
    cudaMalloc ((void**)&(worklist[mat]), lwork*sizeof(float));
    cudaMemcpy(Alist[mat], &A[mat*N*row_stride], N*row_stride*sizeof(float), cudaMemcpyHostToDevice);
  }
  
  cudaMemcpy(Alist_d   ,    Alist, numMats*sizeof(float*), 
	     cudaMemcpyHostToDevice);
  cudaMemcpy(worklist_d, worklist, numMats*sizeof(float*), 
	     cudaMemcpyHostToDevice);
  
  dim3 dimBlock(INVERSE_BS,2);
  dim3 dimGrid(numMats);
  
  clock_t start = clock();
  for (int i=0; i<1; i++) {
    cuda_inverse_many_double (Alist_d, worklist_d, N, row_stride, numMats);
    //    cuda_inverse_many_double (Alist_d, worklist_d, N, numMats);
    cudaThreadSynchronize();
  }

  clock_t end = clock();
  
  double t = (double)(end-start)/(double)CLOCKS_PER_SEC / (double)numMats;
  double rate = 1.0/t;
  fprintf (stderr, "Rate is %1.3f matrix inversions per second.\n",
	   rate);


  cudaError_t err = cudaGetLastError();
  if (err != cudaSuccess) {
    fprintf (stderr, "CUDA error in block_inverse:\n  %s\n",
	     cudaGetErrorString(err));
    abort();
  }

  // Copy Ainv back to host memory
  for (int mat=0; mat<numMats; mat++) {
    float Ainv[N*row_stride];
    cudaMemcpy(Ainv, Alist[mat], N*row_stride*sizeof(float), cudaMemcpyDeviceToHost);
    
    double error = 0.0;
    for (int i=0; i<N; i++)
      for (int j=0; j<N; j++) {
	double val = 0.0;
	for (int k=0; k<N; k++)
	  val += Ainv[i*row_stride+k]*A[mat*N*row_stride+k*row_stride+j];
	double diff = (i==j) ? (1.0f-val) : val;
	error += diff*diff;
      }
    fprintf (stderr, "error = %1.8e\n", sqrt(error/(double)(N*N)));
  }

}






#include <stdio.h>

main()
{
  //test_inverse_many();
  test_inverse_many_double_conv();
  //test_inverse_many_double();

  // int N=32;
  // float A[N*N], Acopy[N*N];
  // float *A_d;
  
  // for (int i=0; i<N; i++)
  //   for (int j=0; j<N; j++)
  //     A[N*i+j] = Acopy[N*i+j] = (float) drand48();

  // cudaMalloc ((void**)&A_d, N*N*sizeof(float));
  // cudaMemcpy (A_d, A, N*N*sizeof(float),
  // 	      cudaMemcpyHostToDevice);

  // dim3 dimBlock(N);
  // dim3 dimGrid(1);
  // block_inverse<float,32><<<dimGrid,dimBlock>>> (A_d, N, N);

  // cudaThreadSynchronize();
  // cudaError_t err = cudaGetLastError();
  // if (err != cudaSuccess) {
  //   fprintf (stderr, "CUDA error in block_inverse:\n  %s\n",
  // 	     cudaGetErrorString(err));
  //   abort();
  // }

  // cudaMemcpy (A, A_d, N*N*sizeof(float),
  // 	      cudaMemcpyDeviceToHost);

  // float nrm = 0.0;
  // for (int i=0; i<N; i++)
  //   for (int j=0; j<N; j++) {
  //     float val = 0.0;
  //     for (int k=0; k<N; k++)
  // 	val += A[i*N+k] * Acopy[k*N+j];
  //     float diff = (i==j) ? 1.0-val : val;
  //     nrm += diff*diff;
  //   }
  // fprintf (stderr, "Error = %1.6e\n", sqrt(nrm/(double)(N*N)));
}
#endif

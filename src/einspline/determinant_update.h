#define DET_BLOCK_SIZE 64

#include <unistd.h>
#include <stdlib.h>


// The first kernel just computes AinvT * u and also stores the kth
// col of Ainv in global memory
__global__ static void
update_inverse_cuda1 (float *A_g[], float *Ainv_g[], float *u_g[],
                      float *Ainv_delta_g[], float *Ainv_colk_g[],
                      int N, int rowstride, int k)
{
  __shared__ float *A, *Ainv, *u, *Ainv_delta, *Ainv_colk;
  if (threadIdx.x==0)
  {
    A        = A_g[blockIdx.y];
    Ainv     = Ainv_g[blockIdx.y];
    u         = u_g[blockIdx.y];
    Ainv_delta    = Ainv_delta_g[blockIdx.y];
    Ainv_colk = Ainv_colk_g[blockIdx.y];
  }
  __syncthreads();
  // Store the product Ainv * u in shared memory
  __shared__ float Ainv_delta_shared[DET_BLOCK_SIZE],
             Ainv_colk_shared[DET_BLOCK_SIZE], u_shared[DET_BLOCK_SIZE],
             uold_shared[DET_BLOCK_SIZE];
  Ainv_delta_shared[threadIdx.x] = 0.0;
  int col = blockIdx.x*DET_BLOCK_SIZE + threadIdx.x;
  int numblocks = N / DET_BLOCK_SIZE;
  // If the column I need to pull from Ainv is in this thread block
  // domain, do the following
  if (blockIdx.x*DET_BLOCK_SIZE <= k && k < (blockIdx.x+1)*DET_BLOCK_SIZE)
  {
    for (int block=0; block<numblocks; block++)
    {
      u_shared[threadIdx.x] = u[block*DET_BLOCK_SIZE+threadIdx.x];
      uold_shared[threadIdx.x] =
        A[k*rowstride + block*DET_BLOCK_SIZE+threadIdx.x];
      // Write new row into A matrix
      A[k*rowstride + block*DET_BLOCK_SIZE+threadIdx.x] = u_shared[threadIdx.x];
      __syncthreads();
      for (int i=0; i<DET_BLOCK_SIZE; i++)
      {
        int row = block*DET_BLOCK_SIZE + i;
        float a = Ainv[row*rowstride+col];
        if (col == k)
          Ainv_colk_shared[i] = a;
        Ainv_delta_shared[threadIdx.x] += a*(u_shared[i]-uold_shared[i]);
      }
      __syncthreads();
      Ainv_colk[block*DET_BLOCK_SIZE+threadIdx.x] = Ainv_colk_shared[threadIdx.x];
    }
  }
  else
  {
    for (int block=0; block<numblocks; block++)
    {
      u_shared[threadIdx.x] = u[block*DET_BLOCK_SIZE+threadIdx.x];
      uold_shared[threadIdx.x] =
        A[k*rowstride + block*DET_BLOCK_SIZE+threadIdx.x];
      // Write new row into A matrix
      A[k*rowstride + block*DET_BLOCK_SIZE+threadIdx.x] = u_shared[threadIdx.x];
      __syncthreads();
      for (int i=0; i<DET_BLOCK_SIZE; i++)
      {
        int row = block*DET_BLOCK_SIZE + i;
        Ainv_delta_shared[threadIdx.x] +=
          Ainv[row*rowstride+col]*(u_shared[i]- uold_shared[i]);
      }
    }
  }
  __syncthreads();
  // Write the data back to global memory
  Ainv_delta[col]    = Ainv_delta_shared[threadIdx.x];
}

__global__ static void
update_inverse_cuda2 (float *Ainv_g[], float *u_g[], float *Ainv_delta_g[],
                      float *Ainv_colk_g[], int N, int rowstride, int k)
{
  __shared__ float *Ainv, *Ainv_delta, *Ainv_colk;
  if (threadIdx.x==0)
  {
    Ainv     = Ainv_g[blockIdx.y];
    Ainv_delta    = Ainv_delta_g[blockIdx.y];
    Ainv_colk = Ainv_colk_g[blockIdx.y];
  }
  __syncthreads();
  __shared__ float Ainv_delta_shared[DET_BLOCK_SIZE];
  __shared__ float  Ainv_colk_shared[DET_BLOCK_SIZE];
  int col = blockIdx.x*DET_BLOCK_SIZE + threadIdx.x;
  // Read the data back from global memory
  Ainv_delta_shared[threadIdx.x] = Ainv_delta[col];
  Ainv_colk_shared[threadIdx.x] = Ainv_colk[col];
  __shared__ float prefact;
  if (threadIdx.x == 0)
    prefact = -1.0f/(1.0f+Ainv_delta[k]);
  __syncthreads();
  int numblocks = N / DET_BLOCK_SIZE;
  for (int block=0; block<numblocks; block++)
  {
    Ainv_colk_shared[threadIdx.x] =
      prefact*Ainv_colk[block*DET_BLOCK_SIZE+threadIdx.x];
    __syncthreads();
    for (int i=0; i<DET_BLOCK_SIZE; i++)
    {
      int row = block*DET_BLOCK_SIZE + i;
      Ainv[row*rowstride+col] +=
        Ainv_delta_shared[threadIdx.x]*Ainv_colk_shared[i];
    }
  }
}

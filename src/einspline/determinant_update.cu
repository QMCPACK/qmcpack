#define BLOCK_SIZE 64

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>


// The first kernel just computes Ainv * u and also stores the kth
// row of Ainv in global memory
__global__ static void
update_inverse_cuda1 (float *Ainv_g[], float *u_g[], float *AinvT_u_g[],
		      float *Ainv_colk_g[], int N, int rowstride, int k)
{
  __shared__ float *Ainv, *u, *AinvT_u, *Ainv_colk;
  if (threadIdx.x==0) {
    Ainv     = Ainv_g[blockIdx.y];
    u         = u_g[blockIdx.y];
    AinvT_u    = AinvT_u_g[blockIdx.y];
    Ainv_colk = Ainv_colk_g[blockIdx.y];
  }

  __syncthreads();

  // Store the product Ainv * u in shared memory
  __shared__ float AinvT_u_shared[BLOCK_SIZE], Ainv_colk_shared[BLOCK_SIZE];
  __shared__ float u_shared[BLOCK_SIZE];
  AinvT_u_shared[threadIdx.x] = 0.0;
  int col = blockIdx.x*BLOCK_SIZE + threadIdx.x;
  int numblocks = N / BLOCK_SIZE;

  if (blockIdx.x*BLOCK_SIZE <= k && k < (blockIdx.x+1)*BLOCK_SIZE) {
    for (int block=0; block<numblocks; block++) {
      u_shared[threadIdx.x] = u[block*BLOCK_SIZE+threadIdx.x];
      __syncthreads();
      for (int i=0; i<BLOCK_SIZE; i++) {
	int row = block*BLOCK_SIZE + i;
	
	float ainv = Ainv[row*rowstride+col];
	if (col == k)
	  Ainv_colk_shared[i] = ainv;
	AinvT_u_shared[threadIdx.x] += ainv*u_shared[i];
      }
      __syncthreads();
      Ainv_colk[block*BLOCK_SIZE+threadIdx.x] = Ainv_colk_shared[threadIdx.x];
    }
  }
  else {
    for (int block=0; block<numblocks; block++) {
      u_shared[threadIdx.x] = u[block*BLOCK_SIZE+threadIdx.x];
      __syncthreads();
      for (int i=0; i<BLOCK_SIZE; i++) {
	int row = block*BLOCK_SIZE + i;
	AinvT_u_shared[threadIdx.x] += Ainv[row*rowstride+col]*u_shared[i];
      }
    }
  }

  __syncthreads();
  
  // Write the data back to global memory
  AinvT_u[col]    = AinvT_u_shared[threadIdx.x];
}

__global__ static void
update_inverse_cuda2 (float *Ainv_g[], float *u_g[], float *AinvT_u_g[],
		      float *Ainv_colk_g[], int N, int rowstride, int k)
{
  __shared__ float *Ainv, *AinvT_u, *Ainv_colk;
  if (threadIdx.x==0) {
    Ainv     = Ainv_g[blockIdx.y];
    AinvT_u    = AinvT_u_g[blockIdx.y];
    Ainv_colk = Ainv_colk_g[blockIdx.y];
  }
  __syncthreads();

  __shared__ float AinvT_u_shared[BLOCK_SIZE];
  __shared__ float  Ainv_colk_shared[BLOCK_SIZE];
  int col = blockIdx.x*BLOCK_SIZE + threadIdx.x;
  // Read the data back from global memory
  AinvT_u_shared[threadIdx.x] = AinvT_u[col];
  Ainv_colk_shared[threadIdx.x] = Ainv_colk[col];
  __shared__ float prefact;
  if (threadIdx.x == 0)
    prefact = -1.0f/(1.0f+AinvT_u[k]);
  __syncthreads();
		   
  int numblocks = N / BLOCK_SIZE;
  for (int block=0; block<numblocks; block++) {
    Ainv_colk_shared[threadIdx.x] = prefact*Ainv_colk[block*BLOCK_SIZE+threadIdx.x];
    __syncthreads();
    for (int i=0; i<BLOCK_SIZE; i++) {
      int row = block*BLOCK_SIZE + i;
      Ainv[row*rowstride+col] += AinvT_u_shared[threadIdx.x]*Ainv_colk_shared[i];
    }
  }
}

#define NMAX 128

__global__ static void
update_inverse_cuda (float *Ainv, float *u, int N, int rowstride, int k)
{
  __shared__ float A_k[NMAX], u_shared[NMAX], Ainv_u[NMAX], Ainv_shared[NMAX];
  A_k[threadIdx.x] = Ainv[k*rowstride+threadIdx.x];
  u_shared[threadIdx.x] = u[threadIdx.x];

  // First, compute k'th element of Ainv_u
  Ainv_u[threadIdx.x] = u_shared[threadIdx.x] * A_k[threadIdx.x];
  __syncthreads();
  for (int n=N>>1; n>0; n = n>>1) {
    float a;
    if (threadIdx.x < n) 
      a = Ainv_u[2*threadIdx.x] + Ainv_u[2*threadIdx.x+1];
    __syncthreads();
    Ainv_u[threadIdx.x] = a;
    __syncthreads();
  }
  float prefact = -1.0f/(1.0f + Ainv_u[0]);

  for (int row=0; row<N; row++) {
    Ainv_shared[threadIdx.x] = Ainv[row*rowstride+threadIdx.x];
    __syncthreads();
    Ainv_u[threadIdx.x] = u_shared[threadIdx.x] * Ainv_shared[threadIdx.x];
    for (int n=N>>1; n>0; n = n>>1) {
      float a;
      if (threadIdx.x < n) 
	a = Ainv_u[2*threadIdx.x] + Ainv_u[2*threadIdx.x+1];
      __syncthreads();
      Ainv_u[threadIdx.x] = a;
      __syncthreads();
    }
    __syncthreads();
    // Now Ainv_u[0] has the row'th element of Ainv_u.
    Ainv[row*rowstride + threadIdx.x] = 
      Ainv_shared[threadIdx.x] + prefact*Ainv_u[0]*A_k[threadIdx.x];
  }

}



// __global__ static void
// update_inverse_cuda (float *AinvT, float *u, int N, int rowstride, int k)
// {
//   // Store the product Ainv * u in shared memory
//   __shared__ float Ainv_u[BLOCK_SIZE], Ainv_u_k[BLOCK_SIZE];
//   Ainv_u[threadIdx.x] = 0.0;
//   __syncthreads();

//   for (int row=0; row < N; row++)
//     Ainv_u[threadIdx.x] += AinvT[row*rowstride+threadIdx.x]*u[row];
  
//   // Compute lambda = [A^(-1)]_k dot u
//   float lambda = 0.0;
//   for (int i=0; i<N; i += BLOCK_SIZE) {
//     Ainv_u_k[threadIdx.x] = AinvT[i+threadIdx.x] * u[i+threadIdx.x];
//     __syncthreads();
//     for (int j=BLOCK_SIZE>>1; j!=0; j >>=1) {
//       if (threadIdx.x < j)
// 	Ainv_u_k[threadIdx.x] = Ainv_u_k[2*threadIdx.x] + Ainv_u_k[2*threadIdx.x+1];
//       lambda += Ainv_u_k[0];
//     }
//     float prefact = 1.0/(1.0+lambda);
//   }

//   // Now, subtract off outer product
// }





void
update_inverse (float *AinvT, float *u, int N, int k)
{
  float Ainv_u[128], Ainv_rowk[128];
  
  for (int i=0; i<N; i++) {
    Ainv_u[i] = 0.0f;
    Ainv_rowk[i] = AinvT[N*i+k];
    for (int j=0; j<N; j++)
      Ainv_u[i] += AinvT[j*N+i] * u[j];
  }

  float prefact = 1.0/(1.0+Ainv_u[k]);

  for (int i=0; i<N; i++)
    for (int j=0; j<N; j++)
      AinvT[j*N+i] -= prefact * Ainv_u[i]*Ainv_rowk[j];
}



// Replaces A with its inverse by gauss-jordan elimination with full pivoting
// Adapted from Numerical Recipes in C
void GJInverse (double *A, int n)
{
  const int maxSize = 2000;

  if (n == 2) { // Special case for 2x2
    double a=A[0]; double b=A[1];
    double c=A[2]; double d=A[3];
    double detInv = 1.0/(a*d-b*c);
    A[0] = d*detInv;
    A[1] = -b*detInv;
    A[2] = -c*detInv;
    A[3] =  a*detInv;
    return;
  }

  int colIndex[maxSize], rowIndex[maxSize], ipiv[maxSize];
  double big, pivInv;
  int icol, irow;
  
  for (int j=0; j<n; j++)
    ipiv[j] = -1;

  for (int i=0; i<n; i++) {
    big = 0.0;
    for (int j=0; j<n; j++) 
      if (ipiv[j] != 0)
	for (int k=0; k<n; k++) {
	  if (ipiv[k] == -1) {
	    if (fabs(A[n*j+k]) >= big) {
	      big = fabs(A[n*j+k]);
	      irow = j; 
	      icol = k;
	    }
	  }
	  else if (ipiv[k] > 0) {
	    fprintf (stderr, "GJInverse: Singular matrix!\n");
	    exit(1);
	  }
	}
    ++(ipiv[icol]); 
    
    if (irow != icol) 
      for (int l=0; l<n; l++) {
	double tmp = A[n*irow+l];
	A[n*irow+l] = A[n*icol+l];
	A[n*icol+l] = tmp;
	// swap (A[n*irow+l], A[n*icol+l]);
      }
			     
    
    rowIndex[i] = irow;
    colIndex[i] = icol;
    if (A[n*icol+icol] == 0.0) { 
      fprintf (stderr, "GJInverse: Singular matrix!\n");
      exit(1);
    }
    pivInv = 1.0/A[n*icol+icol];
    A[n*icol+icol] = 1.0;
    for (int l=0; l<n; l++)
      A[n*icol+l] *= pivInv;
    for (int ll=0; ll<n; ll++)
      if (ll != icol) {
	double dum = A[n*ll+icol];
	A[n*ll+icol] = 0.0;
	for (int l=0; l<n; l++)
	  A[n*ll+l] -= A[n*icol+l]*dum;
      }
  }
  // Now unscramble the permutations
  for (int l=n-1; l>=0; l--) {
    if (rowIndex[l] != colIndex[l])
      for (int k=0; k<n ; k++) {
	double tmp = A[n*k+rowIndex[l]];
	A[n*k+rowIndex[l]] = A[n*k+colIndex[l]];
	A[n*k+colIndex[l]] = tmp;
	// swap (A(k,rowIndex[l]),A(k, colIndex[l]));
      }
  }
}


#define MAT_SIZE 128
#define NUM_MATS 1000

main()
{
  int N = MAT_SIZE;
  double *A, *Ainv;
  int numMats = NUM_MATS;
  float *Ainv_h, *u_h;
  float *Ainv_d, *Ainv_u_d, *Ainv_colk_d, *u_d;


  A       = (double*)malloc (N*N*sizeof(double));
  Ainv    = (double*)malloc (N*N*sizeof(double));
  Ainv_h  = (float*) malloc (N*N*sizeof(float));
  u_h     = (float*) malloc (N*sizeof(float));
  cudaMalloc((void**)&Ainv_d,  N*N*sizeof(float));
  cudaMalloc((void**)&Ainv_d, N*N*sizeof(float));
  cudaMalloc((void**)&u_d, N*sizeof(float));
  cudaMalloc((void**)&Ainv_u_d, N*sizeof(float));
  cudaMalloc((void**)&Ainv_colk_d, N*sizeof(float));
  
  float **AinvList, **Ainv_uList,   
    **Ainv_colkList, **uList;

  AinvList     = (float**)malloc(NUM_MATS*sizeof(float*));
  Ainv_uList    = (float**)malloc(NUM_MATS*sizeof(float*));
  Ainv_colkList = (float**)malloc(NUM_MATS*sizeof(float*));
  uList         = (float**)malloc(NUM_MATS*sizeof(float*));

  float **AinvList_d, **Ainv_uList_d, **Ainv_colkList_d, **uList_d;
  cudaMalloc((void**)&AinvList_d,      numMats*sizeof(float*));
  cudaMalloc((void**)&Ainv_uList_d,    numMats*sizeof(float*));
  cudaMalloc((void**)&Ainv_colkList_d, numMats*sizeof(float*));
  cudaMalloc((void**)&uList_d,         numMats*sizeof(float*));

  fprintf (stderr, "N = %d\n", N);

  
  for (int mat=0; mat<numMats; mat++) {
    cudaMalloc((void**)&(AinvList[mat])  ,   N*N*sizeof(float));
    cudaMalloc((void**)&(Ainv_uList[mat]) ,   N*sizeof(float));
    cudaMalloc((void**)&(Ainv_colkList[mat]), N*sizeof(float));
    cudaMalloc((void**)&(uList[mat])        , N*sizeof(float));
  }

  fprintf (stderr, "N = %d\n", N);


  cudaMemcpy (AinvList_d, AinvList, numMats*sizeof(float*), 
	      cudaMemcpyHostToDevice);
  cudaMemcpy (Ainv_uList_d, Ainv_uList, numMats*sizeof(float*), 
	      cudaMemcpyHostToDevice);
  cudaMemcpy (Ainv_colkList_d, Ainv_colkList, numMats*sizeof(float*), 
	      cudaMemcpyHostToDevice);
  cudaMemcpy (uList_d, uList, numMats*sizeof(float*), 
	      cudaMemcpyHostToDevice);
  
  srand48((long int) 12341313);

  fprintf (stderr, "N = %d\n", N);

  for (int mat=0; mat<numMats; mat++) {
    if (mat == 0 ) {
      for (int i=0; i<N; i++) {
	u_h[i] = drand48();
	for (int j=0; j<N; j++) 
	  A[i*N+j] = Ainv[i*N+j] = drand48();
      }
      GJInverse(Ainv, N);
      for (int i=0; i<N; i++)
	for (int j=0; j<N; j++) 
	  Ainv_h[i*N+j] = (float)Ainv[i*N+j];
    }

    cudaMemcpy (AinvList[mat], Ainv_h, N*N*sizeof(float), 
		cudaMemcpyHostToDevice);
    cudaMemcpy (uList[mat], u_h, N*sizeof(float), cudaMemcpyHostToDevice);
  }

  dim3 dimBlock2(BLOCK_SIZE);
  dim3 dimGrid2(N/BLOCK_SIZE, NUM_MATS);

  int row = 1;


  fprintf (stderr, "Before updates.\n");
  clock_t upStart = clock();
  for (int i=0; i<1; i++) {
    update_inverse_cuda1<<<dimGrid2,dimBlock2>>>
      (AinvList_d, uList_d, Ainv_uList_d, Ainv_colkList_d, N, N, row);
    update_inverse_cuda2<<<dimGrid2,dimBlock2>>>
      (AinvList_d, uList_d, Ainv_uList_d, Ainv_colkList_d, N, N, row);
  }
  clock_t upEnd = clock();
  double uptime = (double)(upEnd - upStart)/(double)CLOCKS_PER_SEC;
  double uprate = (double)N*10*NUM_MATS/uptime;
  fprintf (stderr, "%1.2f updates per second.\n", uprate);
  fprintf (stderr, "%1.3f generations per second.\n", 10.0/uptime);

  cudaMemcpy (Ainv_h, AinvList[1], N*N*sizeof(float),cudaMemcpyDeviceToHost);

  for (int i=0; i<N; i++)
    A[row*N+i] += u_h[i];
  for (int i=0; i<N; i++)
    for (int j=0; j<N; j++) {
      double ident = 0.0;
      for (int k=0; k<N; k++)
  	ident += Ainv_h[i*N+k]*A[k*N+j];
      if ((i==j && fabs(ident - 1.0) > 1.0e-4) ||
  	  (i!=j && fabs(ident) > 1.0e-4))
	fprintf (stderr, "Error in matrix inverse, (%d, %d) = %1.8f\n", i, j, ident);
    }
  fprintf (stderr, "Finished.\n");


//   cudaMemcpy (AinvT_h, AinvT_d, N*N*sizeof(float), cudaMemcpyDeviceToHost);


//   for (int i=0; i<N; i++) {
//     u_h[i] = drand48();
//     for (int j=0; j<N; j++) 
//       A[i*N+j] = Ainv[i*N+j] = drand48();
//   }
  
//   GJInverse(Ainv, N);

//   for (int i=0; i<N; i++)
//     for (int j=0; j<N; j++) {
//       double ident = 0.0;
//       for (int k=0; k<N; k++)
// 	ident += Ainv[i*N+k]*A[k*N+j];
//       if ((i==j && fabs(ident - 1.0) > 1.0e-8) ||
// 	  (i!=j && fabs(ident) > 1.0e-8))
// 	fprintf (stderr, "Error in matrix inverse.\n");
//     }

//   for (int i=0; i<N; i++)
//     for (int j=0; j<N; j++) {
//       AinvT_h[j*N+i] = (float)Ainv[i*N+j];
//       Ainv_h[i*N+j]  = (float)Ainv[i*N+j];
//     }

//   cudaMemcpy (Ainv_d, Ainv_h, N*N*sizeof(float), cudaMemcpyHostToDevice);
//   cudaMemcpy (AinvT_d, AinvT_h, N*N*sizeof(float), cudaMemcpyHostToDevice);
//   cudaMemcpy (u_d, u_h, N*sizeof(float), cudaMemcpyHostToDevice);

//   int col = 1;

//   update_inverse (AinvT_h, u_h, N, col);

//   for (int i=0; i<N; i++)
//     A[i*N+col] += u_h[i];

//   for (int i=0; i<N; i++)
//     for (int j=0; j<N; j++) {
//       double ident = 0.0;
//       for (int k=0; k<N; k++)
// 	ident += AinvT_h[k*N+i]*A[k*N+j];
//       if ((i==j && fabs(ident - 1.0) > 1.0e-4) ||
// 	  (i!=j && fabs(ident) > 1.0e-4))
// 	fprintf (stderr, "Error in matrix inverse, (%d, %d) = %1.8f\n", i, j, ident);
//     }

// //   clock_t host_start = clock();
// //   for (int i=0; i<100000; i++) 
// //     update_inverse (AinvT_h, u_h, N, col);
// //   clock_t host_end = clock();
// //   double host_time = (double)(host_end - host_start)/(double)(CLOCKS_PER_SEC);
// //   double host_rate = 1.0e5/host_time;
// //   fprintf (stderr, "Host rate = %1.8f updates per seconds.\n", host_rate);


//   dim3 dimBlock2(BLOCK_SIZE);
//   dim3 dimGrid2(N/BLOCK_SIZE);

//   update_inverse_cuda1<<<dimGrid2,dimBlock2>>>
//     (AinvT_d, u_d, Ainv_u_d, Ainv_rowk_d, N, N, col);
//   update_inverse_cuda2<<<dimGrid2,dimBlock2>>>
//     (AinvT_d, u_d, Ainv_u_d, Ainv_rowk_d, N, N, col);
//     cudaMemcpy (AinvT_h, AinvT_d, N*N*sizeof(float), cudaMemcpyDeviceToHost);

//  fprintf (stderr, "2 kernel Device test:  ");
//   bool passed = true;
//   for (int i=0; i<N; i++)
//     for (int j=0; j<N; j++) {
//       double ident = 0.0;
//       for (int k=0; k<N; k++)
// 	ident += AinvT_h[k*N+i]*A[k*N+j];
//       if ((i==j && fabs(ident - 1.0) > 1.0e-4) ||
// 	  (i!=j && fabs(ident) > 1.0e-4)) {
// 	fprintf (stderr, "Error in matrix inverse, (%d, %d) = %1.8f\n", i, j, ident);
// 	passed = false;
//       }
//     }
//   if (passed)
//     fprintf (stderr, "Passed.\n");
//   else
//     fprintf (stderr, "Failed.\n");


//   dim3 dimBlock1(MAT_SIZE);
//   dim3 dimGrid1(1);
//   update_inverse_cuda<<<dimGrid1, dimBlock1>>> (Ainv_d, u_d, N, N, col);

//   cudaMemcpy (Ainv_h, Ainv_d, N*N*sizeof(float), cudaMemcpyDeviceToHost);


//   fprintf (stderr, "1-kernel Device test:  ");
//   passed = true;
//   for (int i=0; i<N; i++)
//     for (int j=0; j<N; j++) {
//       double ident = 0.0;
//       for (int k=0; k<N; k++)
// 	//ident += AinvT_h[k*N+i]*A[k*N+j];
// 	ident += Ainv_h[i*N+k]*A[k*N+j];
//       if ((i==j && fabs(ident - 1.0) > 1.0e-4) ||
// 	  (i!=j && fabs(ident) > 1.0e-4)) {
// 	fprintf (stderr, "Error in matrix inverse, (%d, %d) = %1.8f\n", i, j, ident);
// 	passed = false;
//       }
//     }
//   if (passed)
//     fprintf (stderr, "Passed.\n");
//   else
//     fprintf (stderr, "Failed.\n");
    
//   dim3 dimGrid1000(1000);

//   clock_t start = clock();
//   for (int i=0; i<1000; i++)
//     update_inverse_cuda<<<dimGrid1000,dimBlock1>>>
//       (AinvT_d, u_d, N, N, col);
//   clock_t end = clock();

//   double time = (double)(end-start)/(double)CLOCKS_PER_SEC;
//   double rate = 1.0e6/time;

//   fprintf (stderr, "Device rate = %1.8f updates per seconds.\n", rate);



//   // dim3 dimGrid3(N/BLOCK_SIZE, 1000);
//   // dim3 dimGrid4(N/BLOCK_SIZE, 1000);

//   // clock_t start = clock();
//   // for (int i=0; i<1000; i++) {
//   //   update_inverse_cuda1<<<dimGrid3,dimBlock>>>
//   //     (AinvT_d, u_d, Ainv_u_d, Ainv_rowk_d, N, N, col);
//   //   update_inverse_cuda2<<<dimGrid4,dimBlock>>>
//   //     (AinvT_d, u_d, Ainv_u_d, Ainv_rowk_d, N, N, col);
//   // }
//   // clock_t end = clock();

//   // double time = (double)(end-start)/(double)CLOCKS_PER_SEC;
//   // double rate = 1.0e6/time;

//   // fprintf (stderr, "Device rate = %1.8f updates per seconds.\n", rate);


}

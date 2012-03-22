#include "cuda_walker.h"
#include "determinant_update.h"
#include <unistd.h>

cuda_determinant::cuda_determinant() : 
  N(0), A(NULL), Ainv(NULL), Ainv_delta(NULL), Ainv_colk(0),
  new_row(NULL), delta(0)
{

};

cuda_determinant::cuda_determinant(int n)
{
  resize(N);
}

void
cuda_determinant::resize(int n)
{
  N = n;
  cudaMalloc((void**)&A         , N*N*sizeof(float));
  cudaMalloc((void**)&Ainv      , N*N*sizeof(float));
  cudaMalloc((void**)&Ainv_delta, 1*N*sizeof(float));
  cudaMalloc((void**)&Ainv_colk , 1*N*sizeof(float));
  cudaMalloc((void**)&new_row   , 1*N*sizeof(float));
  cudaMalloc((void**)&delta     , 1*N*sizeof(float));
}

void
cuda_walker::resize(int nup, int ndown) 
{
  N[0] = nup; N[1] = ndown;
  dets[0].resize(N[0]);
  dets[1].resize(N[1]);
}



cuda_population::cuda_population() : MaxPop(1000)
{
  A_vec.resize(MaxPop);
  Ainv_vec.resize(MaxPop);
  delta_vec.resize(MaxPop);
  Ainv_delta_vec.resize(MaxPop);
  Ainv_colk_vec.resize(MaxPop);
  ratio_vec.resize(MaxPop);
  pos_vec.resize(3*MaxPop);


  cudaMalloc((void**) &A_list_d,          MaxPop*sizeof(float*));
  cudaMalloc((void**) &Ainv_list_d,       MaxPop*sizeof(float*));
  cudaMalloc((void**) &Ainv_delta_list_d, MaxPop*sizeof(float*));
  cudaMalloc((void**) &Ainv_colk_list_d,  MaxPop*sizeof(float*));
  cudaMalloc((void**) &delta_list_d,      MaxPop*sizeof(float*));
  cudaMalloc((void**) &ratios_d,          MaxPop*sizeof(float));
  cudaMalloc((void**) &pos_d,           4*MaxPop*sizeof(float));
}


__global__ static void
update_inverse_cuda1 (float *A_g[], float *Ainv_g[], float *u_g[], float *Ainv_delta_g[],
		      float *Ainv_colk_g[], int N, int rowstride, int k);
__global__ static void
update_inverse_cuda2 (float *Ainv_g[], float *u_g[], float *Ainv_delta_g[],
		      float *Ainv_colk_g[], int N, int rowstride, int k);


void
cuda_population::calc_new_row(int elec)
{
  int detnum = (elec < num_elecs[0]) ? 0 : 1;
  int N = num_elecs[detnum];
  for (int wi=0; wi<walkers.size(); wi++) {
    cuda_walker &w = walkers[wi];
    cuda_determinant &det = w.dets[detnum];
    pos_vec[4*wi+0] = w.R[3*elec+0];
    pos_vec[4*wi+1] = w.R[3*elec+1];
    pos_vec[4*wi+2] = w.R[3*elec+2];
    delta_vec[wi] = det.delta;
  }
  cudaMemcpy(pos_d, &(pos_vec[0]), 4*walkers.size()*sizeof(float), 
	     cudaMemcpyHostToDevice);
  cudaMemcpy(delta_list_d, &(delta_vec[0]), walkers.size()*sizeof(float*), 
	     cudaMemcpyHostToDevice);

  dim3 dimBlock(SPLINE_BLOCK_SIZE);
  dim3 dimGrid (N/SPLINE_BLOCK_SIZE, walkers.size());
  
  eval_multi_multi_UBspline_3d_s_cuda<<<dimGrid,dimBlock>>>
    (pos_d, multi_spline->gridInv, multi_spline->coefs,
     delta_list_d, multi_spline->stride);

}


void 
cuda_population::update_determinants(int elec)
{
  int index=0;
  int detnum = (elec < num_elecs[0]) ? 0 : 1;
  int N = num_elecs[detnum];
  int row = (elec < num_elecs[0]) ? elec : elec - num_elecs[0];
  for (int wi=0; wi<walkers.size(); wi++) {
    cuda_walker &w = walkers[wi];
    cuda_determinant &det = w.dets[detnum];
    if (w.accept) {
      A_vec[index]          = det.A;
      Ainv_vec[index]       = det.Ainv;
      Ainv_delta_vec[index] = det.Ainv_delta;
      Ainv_colk_vec[index]  = det.Ainv_colk;
      delta_vec[index]      = det.delta;
      index++;
    }
  }
  int num_accept = index;

  cudaMemcpy (A_list_d, &(A_vec[0]), 
	      num_accept*sizeof(float*), cudaMemcpyHostToDevice);
  cudaMemcpy (Ainv_list_d, &(Ainv_vec[0]), 
	      num_accept*sizeof(float*), cudaMemcpyHostToDevice);
  cudaMemcpy (Ainv_delta_list_d, &(Ainv_delta_vec[0]),
	      num_accept*sizeof(float*), cudaMemcpyHostToDevice);
  cudaMemcpy (Ainv_colk_list_d, &(Ainv_colk_vec[0]), 
	      num_accept*sizeof(float*), cudaMemcpyHostToDevice);
  cudaMemcpy (delta_list_d, &(delta_vec[0]), 
	      num_accept*sizeof(float*), cudaMemcpyHostToDevice);

  dim3 dimBlock(DET_BLOCK_SIZE);
  dim3 dimGrid (N/DET_BLOCK_SIZE, num_accept);
  
  update_inverse_cuda1<<<dimGrid,dimBlock>>>
    (A_list_d, Ainv_list_d, delta_list_d, Ainv_delta_list_d, 
     Ainv_colk_list_d, N, N, row);
  update_inverse_cuda2<<<dimGrid,dimBlock>>>
    (Ainv_list_d, delta_list_d, Ainv_delta_list_d, 
     Ainv_colk_list_d, N, N, row);
};

#define RATIO_BLOCK_SIZE 128

__global__ void
calc_ratios1 (float *Ainv_list[], float *new_row_list[],
	      float *Ainv_tran, float *new_row_tran,
	      int N, int k, int row_stride, int num_mats)
{
  int col = threadIdx.x + blockIdx.x*RATIO_BLOCK_SIZE;
  if (col < num_mats) {
    float* Ainv = Ainv_list[col];
    float *new_row = new_row_list[col];
    for (int row=0; row<N; row++) {
      // __shared__ new_row_tran_shared[RATIO_BLOCK_SIZE];
      // __shared__ Ainv_tran_shared[RATIO_BLOCK_SIZE];
      new_row_tran[row_stride*row + col] = new_row[row];
      Ainv_tran[row_stride*row+col] = Ainv[row*N + k];
    }
  }
}


__global__ void
calc_ratios (float *Ainv_list[], float *new_row_list[], 
	     float *ratio, int N, int row_stride, int elec)
{
  int tid = threadIdx.x;

  int col = /*blockIdx.x*RATIO_BLOCK_SIZE * */tid;
  __shared__ float *Ainv, *new_row;

  if (col < N) {
    if (tid == 0) {
      Ainv = Ainv_list[blockIdx.x];
      new_row = new_row_list[blockIdx.x];
    }
    __syncthreads();
    __shared__ float new_row_shared[RATIO_BLOCK_SIZE];
    
    new_row_shared[tid] = new_row[tid];
    
    __shared__ float Ainv_colk_shared[RATIO_BLOCK_SIZE];
    // This is *highly* uncoallesced, but we just have to eat it to allow
    // other kernels to operate quickly.
    Ainv_colk_shared[tid] = Ainv[col*row_stride + elec];
    __syncthreads();

    __shared__ float Ainv_new_row[RATIO_BLOCK_SIZE];
    Ainv_new_row[tid] = Ainv_colk_shared[tid] * new_row_shared[tid];
    
    __syncthreads();
    // Now, we have to dot
    for (unsigned int s=RATIO_BLOCK_SIZE/2; s>0; s>>=1) {
      if (tid < s)
	Ainv_new_row[tid] += Ainv_new_row[tid + s];
      __syncthreads();
    }
    if (tid == 0)      ratio[blockIdx.x] = Ainv_new_row[0];
  }
}


__global__ void
calc_ratios2 (float *Ainv_list[], float *new_row_list[], 
	      float *ratio, int N, int row_stride, int elec)
{
  int tid = threadIdx.x;
  __shared__ float *Ainv, *new_row;
  if (tid == 0) {
    Ainv = Ainv_list[blockIdx.x];
    new_row = new_row_list[blockIdx.x];
  }
  __syncthreads();

  int numBlocks = N/RATIO_BLOCK_SIZE;
  float sum=0.0;
  for (int block=0; block<numBlocks; block++) {
    int row = block*RATIO_BLOCK_SIZE + tid;
    __shared__ float new_row_shared[RATIO_BLOCK_SIZE];
    new_row_shared[tid] = new_row[block*RATIO_BLOCK_SIZE+tid];
    __syncthreads();
    for (int i=0; i<RATIO_BLOCK_SIZE; i++) 
      if (tid==0)
	sum += Ainv[row*row_stride + elec] * new_row_shared[i];
    
  }
  if (tid==0)
    ratio[blockIdx.x] = sum;
}

extern "C" void 
dgetrf_(int *m, int *n, double A[], int *lda, int ipiv[], int *info);

double 
Determinant (double *A, int N)
{
  double LU[N*N];
  int ipiv[N];
  int info;
  for (int i=0; i<N*N; i++)
    LU[i] = A[i];
  // Do LU factorization
  dgetrf_ (&N, &N, LU, &N, ipiv, &info);
  double det = 1.0;
  int numPerm = 0;
  for (int i=0; i<N; i++) {
    det *= LU[i*N+i];
    numPerm += (ipiv[i] != (i+1));
  }
  if (numPerm & 1)
    det *= -1.0;
  
  return det;
}


template<typename T> void 
GJInverse (T *A, int n)
{
  const int maxSize = 2000;

  if (n == 2) { // Special case for 2x2
    T a=A[0]; T b=A[1];
    T c=A[2]; T d=A[3];
    T detInv = 1.0/(a*d-b*c);
    A[0] = d*detInv;
    A[1] = -b*detInv;
    A[2] = -c*detInv;
    A[3] =  a*detInv;
    return;
  }

  int colIndex[maxSize], rowIndex[maxSize], ipiv[maxSize];
  T big, pivInv;
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
	T tmp = A[n*irow+l];
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
	T dum = A[n*ll+icol];
	A[n*ll+icol] = 0.0;
	for (int l=0; l<n; l++)
	  A[n*ll+l] -= A[n*icol+l]*dum;
      }
  }
  // Now unscramble the permutations
  for (int l=n-1; l>=0; l--) {
    if (rowIndex[l] != colIndex[l])
      for (int k=0; k<n ; k++) {
	T tmp = A[n*k+rowIndex[l]];
	A[n*k+rowIndex[l]] = A[n*k+colIndex[l]];
	A[n*k+colIndex[l]] = tmp;
	// swap (A(k,rowIndex[l]),A(k, colIndex[l]));
      }
  }
}




void
test_ratio()
{
  //const int N = RATIO_BLOCK_SIZE;
  const int N = 128;
  const int numWalkers = 1024;
  const int elec = 15;
  float **AinvList, **uList;
  float **AinvList_d, **uList_d, *ratio_d;

  AinvList = (float**) malloc(numWalkers*sizeof(float*));
  uList    = (float**) malloc(numWalkers*sizeof(float*));

  for (int i=0; i<numWalkers; i++) {
    cudaMalloc((void**)&(AinvList[i]), N*N*sizeof(float));
    cudaMalloc((void**)&(uList[i]),      N*sizeof(float));
  }

  fprintf (stderr, "N = %d\n", N);
    
  cudaMalloc((void**)&(AinvList_d), numWalkers*sizeof(float*));
  cudaMalloc((void**)&(uList_d),    numWalkers*sizeof(float*));
  cudaMalloc((void**)&ratio_d,      numWalkers*sizeof(float));

  cudaMemcpy (AinvList_d, AinvList, numWalkers*sizeof(float*), cudaMemcpyHostToDevice);
  cudaMemcpy (   uList_d,    uList, numWalkers*sizeof(float*), cudaMemcpyHostToDevice);

  dim3 dimBlock(RATIO_BLOCK_SIZE);
  dim3 dimGrid(numWalkers);

  double *A   = (double*)malloc(N*N*sizeof(double));
  float *Ainv = (float*) malloc(N*N*sizeof(float));
  float *u    = (float*) malloc(N*sizeof(float));
  for (int i=0; i<N; i++) {
    u[i] = drand48();
    for (int j=0; j<N; j++) 
      A[i*N+j] = Ainv[i*N+j] = (float)drand48();
  }

  GJInverse(Ainv, N);
  double det1 = Determinant (A, N);
  for (int i=0; i<N; i++)
    A[elec*N+i] = u[i];
  double det2 = Determinant (A, N);
  fprintf (stderr, "Host ratio = %1.8f\n", det2/det1);

  for (int wi=0; wi<numWalkers; wi++) {
    cudaMemcpy (AinvList[wi], Ainv, N*N*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy (uList[wi],       u, 1*N*sizeof(float), cudaMemcpyHostToDevice);
  }

  clock_t start = clock();
  for (int i=0; i<10*N; i++) 
    calc_ratios<<<dimGrid,dimBlock>>> (AinvList_d, uList_d, ratio_d, N, N, elec);
  clock_t end = clock();

  float ratio;
  cudaMemcpy (&ratio, &(ratio_d[331]), sizeof(float), cudaMemcpyDeviceToHost);
  fprintf (stderr, "Device ratio = %1.8f\n", ratio);

  
  double time = (double)(end-start)/(double)CLOCKS_PER_SEC;
  double rate = 10.0/time;
  fprintf (stderr, "Rate = %1.3f generations per second.\n", rate);

}


void
test_ratio1()
{
  const int N = 128;
  const int numWalkers = 1024;
  float **AinvList, **uList;
  float **AinvList_d, **uList_d, *ratio_d;
  float *Ainv_tran, *ratio_tran;

  AinvList = (float**) malloc(numWalkers*sizeof(float*));
  uList    = (float**) malloc(numWalkers*sizeof(float*));
  cudaMalloc ((void**) &Ainv_tran, N*numWalkers);
  cudaMalloc ((void**) &ratio_tran, N*numWalkers);

  for (int i=0; i<numWalkers; i++) {
    cudaMalloc((void**)&(AinvList[i]), N*N*sizeof(float));
    cudaMalloc((void**)&(uList[i]),      N*sizeof(float));
  }

  fprintf (stderr, "N = %d\n", N);
    
  cudaMalloc((void**)&(AinvList_d), numWalkers*sizeof(float*));
  cudaMalloc((void**)&(uList_d),    numWalkers*sizeof(float*));
  cudaMalloc((void**)&ratio_d,      numWalkers*sizeof(float));

  cudaMemcpy (AinvList_d, AinvList, numWalkers*sizeof(float*), cudaMemcpyHostToDevice);
  cudaMemcpy (   uList_d,    uList, numWalkers*sizeof(float*), cudaMemcpyHostToDevice);

  dim3 dimBlock(RATIO_BLOCK_SIZE);
  dim3 dimGrid(numWalkers/RATIO_BLOCK_SIZE);

  clock_t start = clock();
  for (int i=0; i<10*N; i++) 
    calc_ratios1<<<dimGrid,dimBlock>>> (AinvList_d, uList_d, Ainv_tran, ratio_tran,
					N, 1, numWalkers, numWalkers);
  clock_t end = clock();
  
  double time = (double)(end-start)/(double)CLOCKS_PER_SEC;
  double rate = 10.0/time;
  fprintf (stderr, "Rate = %1.3f generations per second.\n", rate);

}



main()
{
  test_ratio();
}

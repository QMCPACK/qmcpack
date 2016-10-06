//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
template<typename T, int BS>
__device__
T min_dist_all (T x, T y, T z, 
		T L[3][3], T Linv[3][3], T images[27][3])
{
  int tid = threadIdx.x;

  T u0 = Linv[0][0]*x + Linv[0][1]*y + Linv[0][2]*z;  
  T u1 = Linv[1][0]*x + Linv[1][1]*y + Linv[1][2]*z;
  T u2 = Linv[2][0]*x + Linv[2][1]*y + Linv[2][2]*z;

  u0 -= rintf(u0);
  u1 -= rintf(u1);
  u2 -= rintf(u2);

  x = L[0][0]*u0 + L[0][1]*u1 + L[0][2]*u2;
  y = L[1][0]*u0 + L[1][1]*u1 + L[1][2]*u2;
  z = L[2][0]*u0 + L[2][1]*u1 + L[2][2]*u2;

  __shared__ T dist2[27];
  if (tid < 27) {
    x += images[tid][0];
    y += images[tid][1];
    z += images[tid][2];
    dist2[tid] = x*x + y*y + z*z;
  }
  __syncthreads();
  for (int s=BS>>1; s>0; s>>=1) {
    if (tid < s && (tid+s) < 27)
      dist2[tid] = min(dist2[tid+s],dist2[tid]);
    __syncthreads();
  }

  return sqrtf(dist2[0]);
}

template<typename T>
__device__
T min_dist_only (T x, T y, T z, T L[3][3], T Linv[3][3])
{
  T u0 = Linv[0][0]*x + Linv[0][1]*y + Linv[0][2]*z;  
  T u1 = Linv[1][0]*x + Linv[1][1]*y + Linv[1][2]*z;
  T u2 = Linv[2][0]*x + Linv[2][1]*y + Linv[2][2]*z;

  u0 -= rintf(u0);
  u1 -= rintf(u1);
  u2 -= rintf(u2);

  x = L[0][0]*u0 + L[0][1]*u1 + L[0][2]*u2;
  y = L[1][0]*u0 + L[1][1]*u1 + L[1][2]*u2;
  z = L[2][0]*u0 + L[2][1]*u1 + L[2][2]*u2;

  return sqrtf(x*x + y*y + z*z);
}



// This should be okay for all but the smallest of primitive cells.
// That is, the r_c for the core should be smaller than the simulation
// cell radius
template<typename T>
__device__
T min_dist (T &x, T &y, T &z, T L[3][3], T Linv[3][3])
{
  T u0 = Linv[0][0]*x + Linv[0][1]*y + Linv[0][2]*z;  
  T u1 = Linv[1][0]*x + Linv[1][1]*y + Linv[1][2]*z;
  T u2 = Linv[2][0]*x + Linv[2][1]*y + Linv[2][2]*z;

  u0 -= rintf(u0);
  u1 -= rintf(u1);
  u2 -= rintf(u2);

  x = L[0][0]*u0 + L[0][1]*u1 + L[0][2]*u2;
  y = L[1][0]*u0 + L[1][1]*u1 + L[1][2]*u2;
  z = L[2][0]*u0 + L[2][1]*u1 + L[2][2]*u2;

  return sqrtf(x*x + y*y + z*z);
}




template<typename T, int BS>
__global__ void
find_core_electrons_PBC_kernel(T **R, int numElec,
			       T *I, int firstIon, int lastIon,
			       T rcut, T *L_global, T *Linv_global,
			       int2 **pairs, T **dist, int *numPairs)
{
  int tid = threadIdx.x;
  __shared__ T *myR, *mydist;
  __shared__ int2 *mypairs;
  if (tid == 0) {
    myR     =     R[blockIdx.x];
    mydist  =  dist[blockIdx.x];
    mypairs = pairs[blockIdx.x];
  }
  __shared__ T L[3][3], Linv[3][3];
  if (tid < 9) {
    L[0][tid] = L_global[tid];
    Linv[0][tid] = Linv_global[tid];
  }

  __syncthreads();

//   int i0 = tid / 9;
//   int i1 = (tid - 9*i0)/3;
//   int i2 = (tid - 9*i0 - 3*i1);
//   __syncthreads();


  int numIon = lastIon - firstIon + 1;
  int numElecBlocks = numElec/BS + ((numElec % BS) ? 1 : 0);
  int numIonBlocks  = numIon /BS + ((numIon  % BS) ? 1 : 0);

  __shared__ T r[BS][3];
  __shared__ T i[BS][3];
  __shared__ T d[BS];
  __shared__ int2 blockpairs[BS];
  __shared__ T blockdist[BS];
  int npairs=0, index=0, blockNum=0;


  for (int iBlock=0; iBlock<numIonBlocks; iBlock++) {
    for (int dim=0; dim<3; dim++) 
      if (dim*BS+tid < 3*numIon)
	i[0][dim*BS+tid] = I[3*BS*iBlock + 3*firstIon + dim*BS+tid];
    int ionEnd = ((iBlock+1)*BS < numIon) ? BS : (numIon - iBlock*BS);

    for (int eBlock=0; eBlock<numElecBlocks; eBlock++) {
      int elecEnd = ((eBlock+1)*BS < numElec) ? BS : (numElec - eBlock*BS);
      for (int dim=0; dim<3; dim++) 
	if (dim*BS+tid < 3*numElec)
	  r[0][dim*BS+tid] = myR[3*BS*eBlock + dim*BS+tid];
      for (int ion=0; ion<ionEnd; ion++) {
	d[tid] = min_dist_only(r[tid][0]-i[ion][0], r[tid][1]-i[ion][1],
			       r[tid][2]-i[ion][2], L, Linv);
	for (int elec=0; elec<elecEnd; elec++) {
	  if (d[elec] < rcut) {
	    if (index == BS) {
	      mypairs[blockNum*BS+tid] = blockpairs[tid];
	      mydist[blockNum*BS+tid]  = blockdist[tid];
	      blockNum++;
	      index = 0;
	    }

	    if (tid == 0) {
	      blockpairs[index].x = iBlock*BS+ion;
	      blockpairs[index].y = eBlock*BS+elec;
	      blockdist[index]    = d[tid];
	    }
	    index++;
	    npairs++;
	  }
	}
      }
      __syncthreads();
    }
    __syncthreads();
  }
  // Write pairs and distances remaining the final block
  if (tid < index) {
    mypairs[blockNum*BS+tid] = blockpairs[tid];
    mydist[blockNum*BS+tid]  = blockdist[tid];
  }
  if (tid == 0)
    numPairs[blockIdx.x] = npairs;
}



void
find_core_electrons_PBC (float *R[], int numElec, 
			 float I[], int firstIon, int lastIon,
			 float rcut, float L[], float Linv[], 
			 int2 *pairs[], float *dist[], 
			 int numPairs[], int numWalkers)
{
  const int BS = 32;
  
  dim3 dimBlock(BS);
  dim3 dimGrid(numWalkers);
  
  find_core_electrons_PBC_kernel<float,BS><<<dimGrid,dimBlock>>> 
    (R, numElec, I, firstIon, lastIon, rcut, L, Linv, pairs, dist, numPairs);
}

void
find_core_electrons_PBC (double *R[], int numElec, 
			 double I[], int firstIon, int lastIon,
			 double rcut, double L[], double Linv[], 
			 int2 *pairs[], double *dist[], 
			 int numPairs[], int numWalkers)
{
  const int BS = 32;
  
  dim3 dimBlock(BS);
  dim3 dimGrid(numWalkers);
  
  find_core_electrons_PBC_kernel<double,BS><<<dimGrid,dimBlock>>> 
    (R, numElec, I, firstIon, lastIon, rcut, L, Linv, pairs, dist, numPairs);
}



template<typename T, int BS>
__global__ void
find_core_electrons_PBC_kernel(T **R, int numElec,
			       T *I, int firstIon, int lastIon,
			       T rcut, T *L_global, T *Linv_global,
			       T *quadPoints, int numQuadPoints,
			       int **elecs, T **ratioPos, 
			       T **dist_list, T **cosTheta_list, int *numPairs)
{
  int tid = threadIdx.x;
  __shared__ T *myR, *myRatioPos, *myDist, *myCosTheta;
  __shared__ int *myElecs;
  __shared__ T qp[BS][3];

  for (int i=0; i<3; i++)
    if (i*BS + tid < 3*numQuadPoints)
      qp[0][i*BS+tid] = quadPoints[i*BS+tid];
  if (tid == 0) {
    myR        =        R[blockIdx.x];
    myElecs    =    elecs[blockIdx.x];
    myRatioPos = ratioPos[blockIdx.x];
    myDist     = dist_list[blockIdx.x];
    myCosTheta = cosTheta_list[blockIdx.x];
  }
  __shared__ T L[3][3], Linv[3][3];
  if (tid < 9) {
    L[0][tid] = L_global[tid];
    Linv[0][tid] = Linv_global[tid];
  }

  __syncthreads();

  int numIon = lastIon - firstIon + 1;
  int numElecBlocks = numElec/BS + ((numElec % BS) ? 1 : 0);
  int numIonBlocks  = numIon /BS + ((numIon  % BS) ? 1 : 0);

  __shared__ T r[BS][3];
  __shared__ T i[BS][3];
  __shared__ int blockElecs[BS];
  __shared__ T blockPos[BS][3];
  __shared__ T dist[BS], disp[BS][3];
  __shared__ T blockDist[BS], blockCosTheta[BS];
  int posIndex=0, posBlockNum = 0;
  int npairs=0, index=0, blockNum=0;


  for (int iBlock=0; iBlock<numIonBlocks; iBlock++) {
    for (int dim=0; dim<3; dim++) 
      if ((3*iBlock+dim)*BS+tid < 3*numIon)
	i[0][dim*BS+tid] = I[3*BS*iBlock + 3*firstIon + dim*BS+tid];
    int ionEnd = ((iBlock+1)*BS < numIon) ? BS : (numIon - iBlock*BS);

    for (int eBlock=0; eBlock<numElecBlocks; eBlock++) {
      int elecEnd = ((eBlock+1)*BS < numElec) ? BS : (numElec - eBlock*BS);
      for (int dim=0; dim<3; dim++) 
	if ((3*eBlock+dim)*BS+tid < 3*numElec)
	  r[0][dim*BS+tid] = myR[3*BS*eBlock + dim*BS+tid];
      __syncthreads();
      for (int ion=0; ion<ionEnd; ion++) {
	disp[tid][0] = r[tid][0]-i[ion][0];
	disp[tid][1] = r[tid][1]-i[ion][1];
	disp[tid][2] = r[tid][2]-i[ion][2];
	dist[tid] =  min_dist<T>(disp[tid][0], disp[tid][1], disp[tid][2], L, Linv);
	for (int elec=0; elec<elecEnd; elec++) {
	  __syncthreads();
	  if (dist[elec] < rcut) {
	    // First, write quadrature points
	    if (numQuadPoints + posIndex <= BS) {
	      if (tid < numQuadPoints) {
	    	blockPos[posIndex+tid][0] = r[elec][0] - disp[elec][0] /*i[ion][0]*/ + dist[elec]*qp[tid][0];
		blockPos[posIndex+tid][1] = r[elec][1] - disp[elec][1] /*i[ion][1]*/ + dist[elec]*qp[tid][1];
	    	blockPos[posIndex+tid][2] = r[elec][2] - disp[elec][2] /*i[ion][2]*/ + dist[elec]*qp[tid][2];
		blockCosTheta[posIndex+tid] = 
		  (disp[elec][0]*qp[tid][0] +
		   disp[elec][1]*qp[tid][1] +
		   disp[elec][2]*qp[tid][2]) / dist[elec];
	      }
	      posIndex += numQuadPoints;
	    }
	    else {
	      // Write whatever will fit in the shared buffer
	      int numWrite = BS - posIndex;
	      if (tid < numWrite) {
	    	blockPos[posIndex+tid][0] = r[elec][0] - disp[elec][0] /*i[ion][0]*/ + dist[elec]*qp[tid][0];
	    	blockPos[posIndex+tid][1] = r[elec][1] - disp[elec][1] /*i[ion][1]*/ + dist[elec]*qp[tid][1];
	    	blockPos[posIndex+tid][2] = r[elec][2] - disp[elec][2] /*i[ion][2]*/ + dist[elec]*qp[tid][2];
		blockCosTheta[posIndex+tid] = 
		  (disp[elec][0]*qp[tid][0] +
		   disp[elec][1]*qp[tid][1] +
		   disp[elec][2]*qp[tid][2]) / dist[elec];
	      }
	      __syncthreads();
	      // dump the full buffer to global memory
	      for (int j=0; j<3; j++)
	    	myRatioPos[(posBlockNum*3+j)*BS+tid] = 
	    	  blockPos[0][j*BS+tid];
	      myCosTheta[posBlockNum*BS+tid] = blockCosTheta[tid];
	      posBlockNum++;
	      __syncthreads();
	      // Write the remainder into shared memory
	      if (tid < (numQuadPoints - numWrite)) {
	    	blockPos[tid][0] = r[elec][0] - disp[elec][0] /*i[ion][0]*/ + dist[elec]*qp[tid+numWrite][0];
	    	blockPos[tid][1] = r[elec][1] - disp[elec][1] /*i[ion][1]*/ + dist[elec]*qp[tid+numWrite][1];
	    	blockPos[tid][2] = r[elec][2] - disp[elec][2] /*i[ion][2]*/ + dist[elec]*qp[tid+numWrite][2];
		blockCosTheta[tid] = 
		  (disp[elec][0]*qp[tid+numWrite][0] +
		   disp[elec][1]*qp[tid+numWrite][1] +
		   disp[elec][2]*qp[tid+numWrite][2]) / dist[elec];
	      }
	      posIndex = numQuadPoints - numWrite;
	    }
	    
	    // Now, write electron IDs
	    if (index == BS) {
	      myElecs[blockNum*BS+tid] = blockElecs[tid];
	      myDist[blockNum*BS+tid]  = blockDist[tid];
	      blockNum++;
	      index = 0;
	    }

	    if (tid == 0) {
 	      blockElecs[index] = eBlock*BS+elec;
	      blockDist[index]  = dist[elec];
	    }

	    index++;
	    npairs++;
	    __syncthreads();
	  }
	}
      }
    }
  }
  for (int j=0; j<3; j++)
    if (j*BS + tid < 3*posIndex)
      myRatioPos[(posBlockNum*3+j)*BS + tid] =
	blockPos[0][j*BS+tid];
  if (tid < posIndex)
    myCosTheta[posBlockNum*BS+tid] = blockCosTheta[tid];

  // Write pairs and distances remaining the final block
  if (tid < index) {
    myElecs[blockNum*BS+tid] = blockElecs[tid];
    myDist[blockNum*BS+tid]  = blockDist[tid];
  }
  if (tid == 0)
    numPairs[blockIdx.x] = npairs;
}



void
find_core_electrons_PBC (float *R[], int numElec, 
			 float I[], int firstIon, int lastIon,
			 float rcut, float L[], float Linv[], 
			 float quadPoints[], int numQuadPoints,
			 int *elecs[], float *ratioPos[], 
			 float *dist[], float *cosTheta[],
			 int numPairs[], int numWalkers)
{
 const int BS = 32;
  
  dim3 dimBlock(BS);
  dim3 dimGrid(numWalkers);
  
  find_core_electrons_PBC_kernel<float,BS><<<dimGrid,dimBlock>>> 
    (R, numElec, I, firstIon, lastIon, rcut, L, Linv, 
     quadPoints, numQuadPoints, elecs, ratioPos, dist, cosTheta, numPairs);


}


void
find_core_electrons_PBC (double *R[], int numElec, 
			 double I[], int firstIon, int lastIon,
			 double rcut, double L[], double Linv[], 
			 double quadPoints[], int numQuadPoints,
			 int *elecs[], double *ratioPos[], 
			 double *dist[], double *cosTheta[],
			 int numPairs[], int numWalkers)
{
 const int BS = 32;
  
  dim3 dimBlock(BS);
  dim3 dimGrid(numWalkers);
  
  find_core_electrons_PBC_kernel<double,BS><<<dimGrid,dimBlock>>> 
    (R, numElec, I, firstIon, lastIon, rcut, L, Linv, 
     quadPoints, numQuadPoints, elecs, ratioPos, dist, cosTheta, numPairs);
}





//////////////////////
// Non-PBC versions //
//////////////////////

template<typename T, int BS>
__global__ void
find_core_electrons_kernel(T **R, int numElec,
			   T *I, int firstIon, int lastIon,
			   T rcut, int2 **pairs, T **dist, int *numPairs)
{
  int tid = threadIdx.x;
  __shared__ T *myR, *mydist;
  __shared__ int2 *mypairs;
  if (tid == 0) {
    myR     =     R[blockIdx.x];
    mydist  =  dist[blockIdx.x];
    mypairs = pairs[blockIdx.x];
  }

  __syncthreads();

//   int i0 = tid / 9;
//   int i1 = (tid - 9*i0)/3;
//   int i2 = (tid - 9*i0 - 3*i1);
//   __syncthreads();


  int numIon = lastIon - firstIon + 1;
  int numElecBlocks = numElec/BS + ((numElec % BS) ? 1 : 0);
  int numIonBlocks  = numIon /BS + ((numIon  % BS) ? 1 : 0);

  __shared__ T r[BS][3];
  __shared__ T i[BS][3];
  __shared__ T d[BS];
  __shared__ int2 blockpairs[BS];
  __shared__ T blockdist[BS];
  int npairs=0, index=0, blockNum=0;


  for (int iBlock=0; iBlock<numIonBlocks; iBlock++) {
    for (int dim=0; dim<3; dim++) 
      if (dim*BS+tid < 3*numIon)
	i[0][dim*BS+tid] = I[3*BS*iBlock + 3*firstIon + dim*BS+tid];
    int ionEnd = ((iBlock+1)*BS < numIon) ? BS : (numIon - iBlock*BS);

    for (int eBlock=0; eBlock<numElecBlocks; eBlock++) {
      int elecEnd = ((eBlock+1)*BS < numElec) ? BS : (numElec - eBlock*BS);
      for (int dim=0; dim<3; dim++) 
	if (dim*BS+tid < 3*numElec)
	  r[0][dim*BS+tid] = myR[3*BS*eBlock + dim*BS+tid];
      for (int ion=0; ion<ionEnd; ion++) {
	d[tid] = sqrtf((r[tid][0]-i[ion][0])*(r[tid][0]-i[ion][0]) +
		       (r[tid][0]-i[ion][1])*(r[tid][1]-i[ion][1]) +
		       (r[tid][0]-i[ion][2])*(r[tid][2]-i[ion][2]));
	for (int elec=0; elec<elecEnd; elec++) {
	  if (d[elec] < rcut) {
	    if (index == BS) {
	      mypairs[blockNum*BS+tid] = blockpairs[tid];
	      mydist[blockNum*BS+tid]  = blockdist[tid];
	      blockNum++;
	      index = 0;
	    }

	    if (tid == 0) {
	      blockpairs[index].x = iBlock*BS+ion;
	      blockpairs[index].y = eBlock*BS+elec;
	      blockdist[index]    = d[tid];
	    }
	    index++;
	    npairs++;
	  }
	}
      }
      __syncthreads();
    }
    __syncthreads();
  }
  // Write pairs and distances remaining the final block
  if (tid < index) {
    mypairs[blockNum*BS+tid] = blockpairs[tid];
    mydist[blockNum*BS+tid]  = blockdist[tid];
  }
  if (tid == 0)
    numPairs[blockIdx.x] = npairs;
}



void
find_core_electrons (float *R[], int numElec, 
		     float I[], int firstIon, int lastIon,
		     float rcut, int2 *pairs[], float *dist[], 
		     int numPairs[], int numWalkers)
{
  const int BS = 32;
  
  dim3 dimBlock(BS);
  dim3 dimGrid(numWalkers);
  
  find_core_electrons_kernel<float,BS><<<dimGrid,dimBlock>>> 
    (R, numElec, I, firstIon, lastIon, rcut, pairs, dist, numPairs);
}

void
find_core_electrons (double *R[], int numElec, 
		     double I[], int firstIon, int lastIon,
		     double rcut, int2 *pairs[], double *dist[], 
		     int numPairs[], int numWalkers)
{
  const int BS = 32;
  
  dim3 dimBlock(BS);
  dim3 dimGrid(numWalkers);
  
  find_core_electrons_kernel<double,BS><<<dimGrid,dimBlock>>> 
    (R, numElec, I, firstIon, lastIon, rcut, pairs, dist, numPairs);
}



template<typename T, int BS>
__global__ void
find_core_electrons_kernel(T **R, int numElec,
			   T *I, int firstIon, int lastIon,
			   T rcut, T *quadPoints, int numQuadPoints,
			   int **elecs, T **ratioPos, 
			   T **dist_list, T **cosTheta_list, int *numPairs)
{
  int tid = threadIdx.x;
  __shared__ T *myR, *myRatioPos, *myDist, *myCosTheta;
  __shared__ int *myElecs;
  __shared__ T qp[BS][3];

  for (int i=0; i<3; i++)
    if (i*BS + tid < 3*numQuadPoints)
      qp[0][i*BS+tid] = quadPoints[i*BS+tid];
  if (tid == 0) {
    myR        =        R[blockIdx.x];
    myElecs    =    elecs[blockIdx.x];
    myRatioPos = ratioPos[blockIdx.x];
    myDist     = dist_list[blockIdx.x];
    myCosTheta = cosTheta_list[blockIdx.x];
  }

  __syncthreads();

  int numIon = lastIon - firstIon + 1;
  int numElecBlocks = numElec/BS + ((numElec % BS) ? 1 : 0);
  int numIonBlocks  = numIon /BS + ((numIon  % BS) ? 1 : 0);

  __shared__ T r[BS][3];
  __shared__ T i[BS][3];
  __shared__ int blockElecs[BS];
  __shared__ T blockPos[BS][3];
  __shared__ T dist[BS], disp[BS][3];
  __shared__ T blockDist[BS], blockCosTheta[BS];
  int posIndex=0, posBlockNum = 0;
  int npairs=0, index=0, blockNum=0;


  for (int iBlock=0; iBlock<numIonBlocks; iBlock++) {
    for (int dim=0; dim<3; dim++) 
      if ((3*iBlock+dim)*BS+tid < 3*numIon)
	i[0][dim*BS+tid] = I[3*BS*iBlock + 3*firstIon + dim*BS+tid];
    int ionEnd = ((iBlock+1)*BS < numIon) ? BS : (numIon - iBlock*BS);

    for (int eBlock=0; eBlock<numElecBlocks; eBlock++) {
      int elecEnd = ((eBlock+1)*BS < numElec) ? BS : (numElec - eBlock*BS);
      for (int dim=0; dim<3; dim++) 
	if ((3*eBlock+dim)*BS+tid < 3*numElec)
	  r[0][dim*BS+tid] = myR[3*BS*eBlock + dim*BS+tid];
      __syncthreads();
      for (int ion=0; ion<ionEnd; ion++) {
	disp[tid][0] = r[tid][0]-i[ion][0];
	disp[tid][1] = r[tid][1]-i[ion][1];
	disp[tid][2] = r[tid][2]-i[ion][2];
	dist[tid] =  sqrtf(disp[tid][0]*disp[tid][0] +
			   disp[tid][1]*disp[tid][1] +
			   disp[tid][2]*disp[tid][2]);
	for (int elec=0; elec<elecEnd; elec++) {
	  __syncthreads();
	  if (dist[elec] < rcut) {
	    // First, write quadrature points
	    if (numQuadPoints + posIndex <= BS) {
	      if (tid < numQuadPoints) {
	    	blockPos[posIndex+tid][0] = r[elec][0] - disp[elec][0] /*i[ion][0]*/ + dist[elec]*qp[tid][0];
		blockPos[posIndex+tid][1] = r[elec][1] - disp[elec][1] /*i[ion][1]*/ + dist[elec]*qp[tid][1];
	    	blockPos[posIndex+tid][2] = r[elec][2] - disp[elec][2] /*i[ion][2]*/ + dist[elec]*qp[tid][2];
		blockCosTheta[posIndex+tid] = 
		  (disp[elec][0]*qp[tid][0] +
		   disp[elec][1]*qp[tid][1] +
		   disp[elec][2]*qp[tid][2]) / dist[elec];
	      }
	      posIndex += numQuadPoints;
	    }
	    else {
	      // Write whatever will fit in the shared buffer
	      int numWrite = BS - posIndex;
	      if (tid < numWrite) {
	    	blockPos[posIndex+tid][0] = r[elec][0] - disp[elec][0] /*i[ion][0]*/ + dist[elec]*qp[tid][0];
	    	blockPos[posIndex+tid][1] = r[elec][1] - disp[elec][1] /*i[ion][1]*/ + dist[elec]*qp[tid][1];
	    	blockPos[posIndex+tid][2] = r[elec][2] - disp[elec][2] /*i[ion][2]*/ + dist[elec]*qp[tid][2];
		blockCosTheta[posIndex+tid] = 
		  (disp[elec][0]*qp[tid][0] +
		   disp[elec][1]*qp[tid][1] +
		   disp[elec][2]*qp[tid][2]) / dist[elec];
	      }
	      __syncthreads();
	      // dump the full buffer to global memory
	      for (int j=0; j<3; j++)
	    	myRatioPos[(posBlockNum*3+j)*BS+tid] = 
	    	  blockPos[0][j*BS+tid];
	      myCosTheta[posBlockNum*BS+tid] = blockCosTheta[tid];
	      posBlockNum++;
	      __syncthreads();
	      // Write the remainder into shared memory
	      if (tid < (numQuadPoints - numWrite)) {
	    	blockPos[tid][0] = r[elec][0] - disp[elec][0] /*i[ion][0]*/ + dist[elec]*qp[tid+numWrite][0];
	    	blockPos[tid][1] = r[elec][1] - disp[elec][1] /*i[ion][1]*/ + dist[elec]*qp[tid+numWrite][1];
	    	blockPos[tid][2] = r[elec][2] - disp[elec][2] /*i[ion][2]*/ + dist[elec]*qp[tid+numWrite][2];
		blockCosTheta[tid] = 
		  (disp[elec][0]*qp[tid+numWrite][0] +
		   disp[elec][1]*qp[tid+numWrite][1] +
		   disp[elec][2]*qp[tid+numWrite][2]) / dist[elec];
	      }
	      posIndex = numQuadPoints - numWrite;
	    }
	    
	    // Now, write electron IDs
	    if (index == BS) {
	      myElecs[blockNum*BS+tid] = blockElecs[tid];
	      myDist[blockNum*BS+tid]  = blockDist[tid];
	      blockNum++;
	      index = 0;
	    }

	    if (tid == 0) {
 	      blockElecs[index] = eBlock*BS+elec;
	      blockDist[index]  = dist[elec];
	    }

	    index++;
	    npairs++;
	    __syncthreads();
	  }
	}
      }
    }
  }
  for (int j=0; j<3; j++)
    if (j*BS + tid < 3*posIndex)
      myRatioPos[(posBlockNum*3+j)*BS + tid] =
	blockPos[0][j*BS+tid];
  if (tid < posIndex)
    myCosTheta[posBlockNum*BS+tid] = blockCosTheta[tid];

  // Write pairs and distances remaining the final block
  if (tid < index) {
    myElecs[blockNum*BS+tid] = blockElecs[tid];
    myDist[blockNum*BS+tid]  = blockDist[tid];
  }
  if (tid == 0)
    numPairs[blockIdx.x] = npairs;
}



void
find_core_electrons (float *R[], int numElec, 
		     float I[], int firstIon, int lastIon,
		     float rcut, float quadPoints[], int numQuadPoints,
		     int *elecs[], float *ratioPos[], 
		     float *dist[], float *cosTheta[],
		     int numPairs[], int numWalkers)
{
 const int BS = 32;
  
  dim3 dimBlock(BS);
  dim3 dimGrid(numWalkers);
  
  find_core_electrons_kernel<float,BS><<<dimGrid,dimBlock>>> 
    (R, numElec, I, firstIon, lastIon, rcut, quadPoints, numQuadPoints, 
     elecs, ratioPos, dist, cosTheta, numPairs);


}


void
find_core_electrons (double *R[], int numElec, 
		     double I[], int firstIon, int lastIon,
		     double rcut, 
		     double quadPoints[], int numQuadPoints,
		     int *elecs[], double *ratioPos[], 
		     double *dist[], double *cosTheta[],
		     int numPairs[], int numWalkers)
{
 const int BS = 32;
  
  dim3 dimBlock(BS);
  dim3 dimGrid(numWalkers);
  
  find_core_electrons_kernel<double,BS><<<dimGrid,dimBlock>>> 
    (R, numElec, I, firstIon, lastIon, rcut, 
     quadPoints, numQuadPoints, elecs, ratioPos, dist, cosTheta, numPairs);
}
















// Maximum quadrature points of 32;


// This kernel assumes that the pair are sorted according to ion
// number
template<typename T, int BS>
__global__ void
make_work_list_kernel (int2 **pairs, T **dist, int *numPairs,
		       T *I, int numIons, T *quadPoints, int numQuadPoints,
		       T **ratio_pos)
{
  __shared__ T qp[BS][3];
  __shared__ T *myPairs, *myDist, *myRatioPos;
  __shared__ int np;
  int tid = threadIdx.x;
  if (tid == 0) {
    myPairs = pairs[blockIdx.x];
    myDist =   dist[blockIdx.x];
    myRatioPos = ratio_pos[blockIdx.x];
    np = numPairs[blockIdx.x];
  }

  for (int i=0; i<3; i++)
    if (i*BS+tid < 3*numQuadPoints)
      qp[0][i*BS + tid] = quadPoints[i*BS + tid];
  __syncthreads();

  __shared__ int2 sharedPairs[BS];
  __shared__ T i[BS][3];

  int iBlock = 0;
  int numPairBlocks = np/BS + ((np % BS) ? 1 : 0);
  for (int pairBlock=0; pairBlock<numPairBlocks; pairBlock++) {
    if (pairBlock*BS + tid < np)
      sharedPairs[tid] = myPairs[pairBlock*BS+tid];
    int end = ((pairBlock+1)*BS < np) ? BS : (np - pairBlock*BS);
    for (int ip=0; ip<end; ip++) {
      if (iBlock*BS < sharedPairs[ip].x) {
	while ((iBlock++*BS) < sharedPairs[ip].x);
	for (int dim=0; dim<3; dim++) 
	  if (dim*BS+tid < 3*numIons)
	    i[0][dim*BS+tid] = I[3*BS*iBlock + dim*BS+tid];
      }
      
      
    }
  }
}

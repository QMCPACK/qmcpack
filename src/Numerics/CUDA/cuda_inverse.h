#ifndef CUDA_INVERSE_H
#define CUDA_INVERSE_H

//////////////////////
// Single precision //
//////////////////////
size_t
cuda_inverse_many_worksize(int N);

void
cuda_inverse_many (float *Alist_d[], float *worklist_d[],
		   int N, int num_mats);



//////////////////////
// Double precision //
//////////////////////
size_t
cuda_inverse_many_double_worksize(int N);

void
cuda_inverse_many_double (float *Alist_d[], float *worklist_d[],
			  int N, int num_mats);
void
cuda_inverse_many_double (float *Alist_d[], float *worklist_d[],
			  int N, int row_stride, int num_mats);

void
cuda_inverse_many_double (double *Alist_d[], double *worklist_d[],
			  int N, int num_mats);




#endif

#ifndef ACCEPT_KERNEL_H
#define ACCEPT_KERNEL_H

void
accept_move_GPU_cuda (float* Rlist[], float new_pos[],
                      int toAccept[], int iat, int N);
void
accept_move_GPU_cuda (double* Rlist[], double new_pos[],
                      int toAccept[], int iat, int N);

void NL_move_cuda ( float* Rlist[],  float new_pos[], int N);
void NL_move_cuda (double* Rlist[], double new_pos[], int N);

#endif

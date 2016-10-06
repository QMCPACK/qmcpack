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

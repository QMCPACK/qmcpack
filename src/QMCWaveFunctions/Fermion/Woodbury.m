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
    
    
N = 512
M = 16

A = rand(N,N);

V = rand(M,N);
U = [eye(M,M);zeros(N-M,M)];

Ainv = A^(-1);
Anew = A + U*V;
Ainv_new = Anew^(-1);

Ainv_WF = Ainv - Ainv*U*(eye(M,M) + V*Ainv*U)^(-1)*V*Ainv;

max(max(abs(Ainv_WF - Ainv_new)))
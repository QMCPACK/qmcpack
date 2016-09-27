//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, StoneRidge Inc.
//                    Jeremy McMinnis, jmcminis@gmail.com, Navar Inc.
//
// File created by: Ken Esler, kpesler@gmail.com, StoneRidge Inc.
//////////////////////////////////////////////////////////////////////////////////////
    
    
#ifndef NLJOB_GPU_H
#define NLJOB_GPU_H

template <typename S>
struct NLjobGPU
{
  int Elec, NumQuadPoints;
  S *R, *QuadPoints, *Ratios;
};

#endif

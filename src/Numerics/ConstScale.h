//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    




#ifndef OHMMS_CONSTSCALINGSKMODEL_H
#define OHMMS_CONSTSCALINGSKMODEL_H

// SK<T,L1,L2,MID> where MID = 0 for no scaling
// dum scaling function which returns a constant function
struct ConstScale
{
  double C;
  ConstScale(double c=1.0):C(c) { }
  ~ConstScale() { }
  inline double operator()(double r)
  {
    return C;
  }
  inline double operator()(double r, double& vr)
  {
    vr = C/r;
    return 0.0e0;
  }
  inline double operator()(double r, double& vr, double& dvr)
  {
    vr = C/r/r;
    dvr = 0.0e0;
    return 0.0e0;
  }
};

#endif

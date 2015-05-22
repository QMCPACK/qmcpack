//////////////////////////////////////////////////////////////////
// (c) Copyright 2009-  by Jeongnim Kim and Ken Esler           //
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &          //
//   Materials Computation Center                               //
//   University of Illinois, Urbana-Champaign                   //
//   Urbana, IL 61801                                           //
//   e-mail: jnkim@ncsa.uiuc.edu                                //
//                                                              //
// Supported by                                                 //
//   National Center for Supercomputing Applications, UIUC      //
//   Materials Computation Center, UIUC                         //
//////////////////////////////////////////////////////////////////

#ifndef EINSPLINE_INTERFACES_H
#define EINSPLINE_INTERFACES_H
#include "OhmmsPETE/Tensor.h"
#include "OhmmsPETE/TinyVector.h"	
#include "OhmmsPETE/OhmmsArray.h"	
#include <einspline/bspline.h>
#include <einspline/nubspline.h>
#include <einspline/multi_bspline_structs.h>

#if   OHMMS_DIM==3

#define EINSPLINE_3D_INTERFACE(SPT,RT,DT)                                   \
inline void create_bspline(SPT *Spline, Array<DT,3> &data                   \
    , TinyVector<RT,3>& start, TinyVector<RT,3>& end                        \
    , BCtype_d xBC, BCtype_d yBC, BCtype_d zBC=0)                           \
{                                                                           \
  if (Spline != NULL) destroy_Bspline (Spline);                             \
  Ugrid xgrid, ygrid, zgrid;                                                \
  xgrid.start = start[0];   xgrid.end = end[0];  xgrid.num =data.size(0);   \
  ygrid.start = start[1];   ygrid.end = end[1];  ygrid.num =data.size(1);   \
  zgrid.start = start[2];   zgrid.end = end[2];  zgrid.num =data.size(2);   \
  Spline = create_SPT(xgrid, ygrid, zgrid, xBC, yBC, zBC, data.data());     \
}                                                                           \
                                                                            \
inline DT evaluate(SPT *Spline, const TinyVector<RT,3>& r)                  \
{                                                                           \
  DT val;                                                                   \
  eval_SPT(Spline, r[0], r[1], r[2], &val);                                 \
  return val;                                                               \
}                                                                           \
                                                                            \
inline void evaluate(SPT *Spline, const TinyVector<RT,3>& r                 \
    , DT &val, TinyVector<DT,3> &grad)                                      \
{                                                                           \
  eval_SPT_vg (Spline, r[0], r[1], r[2], &val, &grad[0]);                   \
}                                                                           \
                                                                            \
inline void evaluate(SPT *Spline, const TinyVector<RT,3>& r                 \
    , DT& val, TinyVector<DT,3> &grad, DT &lapl)                            \
{                                                                           \
  eval_SPT_vgl (Spline, r[0], r[1], r[2], &val, &grad[0], &lapl);           \
}                                                                           \
                                                                            \
inline void evaluate(SPT *Spline, const TinyVector<RT,3>& r                 \
    , DT &val, TinyVector<DT,3> &grad, Tensor<DT,3> &hess)                  \
{                                                                           \
  eval_SPT_vgh (Spline, r[0], r[1], r[2], &val, &grad[0],hess.begin());     \
}                                                                           \
                                                                            \
inline void evaluate(multi_SPT *spline                                      \
    , const TinyVector<RT,3>& r, Vector<DT> &psi)                           \
{                                                                           \
  eval_multi_SPT (spline, r[0], r[1], r[2], psi.data());                    \
}                                                                           \
                                                                            \
inline void evaluate(multi_SPT *spline, const TinyVector<RT,3>& r           \
    , Vector<DT> &psi, Vector<TinyVector<DT,3> > &grad)                     \
{                                                                           \
  eval_multi_SPT_vg (spline, r[0], r[1], r[2], psi.data(), &(grad[0][0]));  \
}                                                                           \
                                                                            \
inline void evaluate(multi_SPT *spline, const TinyVector<RT,3>& r           \
    , Vector<DT> &psi, Vector<TinyVector<DT,3> > &grad                      \
    , Vector<Tensor<DT,3> > &hess)                                          \
{                                                                           \
  eval_multi_SPT_vgh (spline, r[0], r[1], r[2], psi.data(), &(grad[0][0]),&(hess[0][0]));\
}                                                                           \

//declare inline functions 
namespace qmcplusplus { namespace einspline {

EINSPLINE_3D_INTERFACE(UBspline_3d_d,double,double)
EINSPLINE_3D_INTERFACE(UBspline_3d_z,double,complex<double>)
EINSPLINE_3D_INTERFACE(NUBspline_3d_d,double,double)
EINSPLINE_3D_INTERFACE(NUBspline_3d_z,double,complex<double>)
}}

#elif OHMMS_DIM==2

#define EINSPLINE_2D_INTERFACE(SPT,RT,DT)                                   \
  inline void create_bspline(SPT *Spline, Array<DT,2> &data                 \
      , TinyVector<RT,2>& start, TinyVector<RT,2>& end                      \
      , BCtype_d xBC, BCtype_d yBC, BCtype_d zBC=0)                         \
{                                                                           \
  if (Spline != NULL) destroy_Bspline (Spline);                             \
  Ugrid xgrid, ygrid;                                                       \
  xgrid.start = start[0];   xgrid.end = end[0];  xgrid.num =data.size(0);   \
  ygrid.start = start[1];   ygrid.end = end[1];  ygrid.num =data.size(1);   \
  Spline = create_SPT_2d_d (xgrid, ygrid, xBC, yBC, data.data());           \
}                                                                           \
                                                                            \
inline DT evaluate(SPT *Spline, const TinyVector<RT,2>& r)                  \
{                                                                           \
  DT val;                                                                   \
  eval_SPT (Spline, r[0], r[1], &val);                                      \
  return val;                                                               \
}                                                                           \
                                                                            \
inline void evaluate(SPT *Spline, const TinyVector<RT,2>& r                 \
    , DT &val, TinyVector<DT,2> &grad)                                      \
{                                                                           \
  eval_SPT_vg (Spline, r[0], r[1], &val, &grad[0]);                         \
}                                                                           \
                                                                            \
inline void evaluate(SPT *Spline, const TinyVector<RT,2>& r                 \
    , DT& val, TinyVector<DT,2> &grad, DT &lapl)                            \
{                                                                           \
  eval_SPT_vgl (Spline, r[0], r[1], &val, &grad[0], &lapl);                 \
}                                                                           \
                                                                            \
inline void evaluate(SPT *Spline, const TinyVector<RT,2>& r                 \  
    , DT &val, TinyVector<DT,2> &grad, Tensor<DT,2> &hess)                  \
{                                                                           \
  eval_SPT_vgh (Spline, r[0], r[1], &val, &grad[0], &hess(0,0));            \
}                                                                           \

namespace qmcplusplus { namespace einspline{
EINSPLINE_2D_INTERFACE(UBspline_2d_d,double,double)
EINSPLINE_2D_INTERFACE(UBspline_2d_z,double,complex<double>)
EINSPLINE_2D_INTERFACE(NUBspline_2d_d,double,double)
EINSPLINE_2D_INTERFACE(NUBspline_2d_z,double,complex<double>)
}}
#endif

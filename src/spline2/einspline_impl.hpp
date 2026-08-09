//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2026 QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


/** @file einspline_impl.hpp
 * @brief interface to hide types in einspline library
 *
 */
#ifndef QMCPLUSPLUS_EINSPLINE_IMPL_H
#define QMCPLUSPLUS_EINSPLINE_IMPL_H

#include "CPU/SIMD/algorithm.hpp"
#include "einspline/multi_bspline.h"
#include "einspline/multi_bspline_copy.h"

namespace qmcplusplus
{
/** einspline name space to define functions to handle functions
   *
   * functions to handle einspline calls transparently
   * For datatype (double,complex<double>,float,complex<float>)
   *  - create(spline,start,end,bc,num_splines)
   *  - set(spline,i,data)
   *  - evaluate(spline,r,psi)
   *  - evaluate(spline,r,psi,grad)
   *  - evaluate(spline,r,psi,grad,lap)
   *  - evaluate(spline,r,psi,grad,hess)
   * are defined to wrap einspline calls. A similar pattern is used for BLAS/LAPACK.
   * The template parameters of the functions  are
   * \tparam PT position type, e.g. TinyVector<T,D>
   * \tparam VT array of values, e.g. Vector<T>
   * \tparam GT array of gradients, e.g. Vector<TinyVector<T,D> >
   * \tparam HT array of hessian tensors, e.g. Vector<Tensor<T,D> >
   */
namespace einspline
{
/** create spline for double */
template<typename GT, typename BCT>
multi_UBspline_3d_d* create(multi_UBspline_3d_d* s, GT& grid, BCT& bc, int num_splines)
{
  return create_multi_UBspline_3d_d(grid[0], grid[1], grid[2], bc[0], bc[1], bc[2], num_splines);
}

/** create spline for float */
template<typename GT, typename BCT>
multi_UBspline_3d_s* create(multi_UBspline_3d_s* s, GT& grid, BCT& bc, int num_splines)
{
  return create_multi_UBspline_3d_s(grid[0], grid[1], grid[2], bc[0], bc[1], bc[2], num_splines);
}

/////another creation functions
/** create spline and initialized it */
template<typename VT, typename IT>
multi_UBspline_3d_s* create(multi_UBspline_3d_s* s, VT& start, VT& end, IT& ng, bc_code bc, int num_splines)
{
  Ugrid x_grid, y_grid, z_grid;
  BCtype_s xBC, yBC, zBC;
  x_grid.start = start[0];
  x_grid.end   = end[0];
  x_grid.num   = ng[0];
  y_grid.start = start[1];
  y_grid.end   = end[1];
  y_grid.num   = ng[1];
  z_grid.start = start[2];
  z_grid.end   = end[2];
  z_grid.num   = ng[2];
  xBC.lCode = xBC.rCode = bc;
  yBC.lCode = yBC.rCode = bc;
  zBC.lCode = zBC.rCode = bc;
  return create_multi_UBspline_3d_s(x_grid, y_grid, z_grid, xBC, yBC, zBC, num_splines);
}

/** convert double to single precision */
inline void convert(multi_UBspline_3d_d* in, multi_UBspline_3d_s* out)
{
  BCtype_s xbc, ybc, zbc;
  xbc.lCode = in->xBC.lCode;
  xbc.rCode = in->xBC.rCode;
  ybc.lCode = in->yBC.lCode;
  ybc.rCode = in->yBC.rCode;
  zbc.lCode = in->zBC.lCode;
  zbc.rCode = in->zBC.rCode;

  xbc.lVal = (float)(in->xBC.lVal);
  xbc.rVal = (float)(in->xBC.rVal);
  ybc.lVal = (float)(in->yBC.lVal);
  ybc.rVal = (float)(in->yBC.rVal);
  zbc.lVal = (float)(in->zBC.lVal);
  zbc.rVal = (float)(in->zBC.rVal);

  out = create_multi_UBspline_3d_s(in->x_grid, in->y_grid, in->z_grid, xbc, ybc, zbc, in->num_splines);
  simd::copy(out->coefs, in->coefs, in->coefs_size);
}
/** create  multi_UBspline_3d_d*
     * @param s dummy multi_UBspline_3d_d* 
     * @param start starting grid values
     * @param end ending grid values
     * @param ng number of grids for [start,end)
     * @param bc boundary condition
     * @param num_splines number of splines to do 
     */
template<typename VT, typename IT>
multi_UBspline_3d_d* create(multi_UBspline_3d_d* s, VT& start, VT& end, IT& ng, bc_code bc, int num_splines)
{
  Ugrid x_grid, y_grid, z_grid;
  BCtype_d xBC, yBC, zBC;
  x_grid.start = start[0];
  x_grid.end   = end[0];
  x_grid.num   = ng[0];
  y_grid.start = start[1];
  y_grid.end   = end[1];
  y_grid.num   = ng[1];
  z_grid.start = start[2];
  z_grid.end   = end[2];
  z_grid.num   = ng[2];
  xBC.lCode = xBC.rCode = bc;
  yBC.lCode = yBC.rCode = bc;
  zBC.lCode = zBC.rCode = bc;
  return create_multi_UBspline_3d_d(x_grid, y_grid, z_grid, xBC, yBC, zBC, num_splines);
}

template<typename VT, typename IT>
multi_UBspline_3d_d* create(multi_UBspline_3d_d* s,
                            VT& start,
                            VT& end,
                            IT& ng,
                            bc_code xbc,
                            bc_code ybc,
                            bc_code zbc,
                            int num_splines)
{
  Ugrid x_grid, y_grid, z_grid;
  BCtype_d xBC, yBC, zBC;
  x_grid.start = start[0];
  x_grid.end   = end[0];
  x_grid.num   = ng[0];
  y_grid.start = start[1];
  y_grid.end   = end[1];
  y_grid.num   = ng[1];
  z_grid.start = start[2];
  z_grid.end   = end[2];
  z_grid.num   = ng[2];
  xBC.lCode = xBC.rCode = xbc;
  yBC.lCode = yBC.rCode = ybc;
  zBC.lCode = zBC.rCode = zbc;
  return create_multi_UBspline_3d_d(x_grid, y_grid, z_grid, xBC, yBC, zBC, num_splines);
}

/** interfaces to use UBspline_3d_X 
     *
     * - create
     * - set
     * - evaluate
     */
template<typename VT, typename IT>
UBspline_3d_d* create(UBspline_3d_d* s, VT& start, VT& end, IT& ng, bc_code bc, int n = 1)
{
  Ugrid x_grid, y_grid, z_grid;
  BCtype_d xBC, yBC, zBC;
  x_grid.start = start[0];
  x_grid.end   = end[0];
  x_grid.num   = ng[0];
  y_grid.start = start[1];
  y_grid.end   = end[1];
  y_grid.num   = ng[1];
  z_grid.start = start[2];
  z_grid.end   = end[2];
  z_grid.num   = ng[2];
  xBC.lCode = xBC.rCode = bc;
  yBC.lCode = yBC.rCode = bc;
  zBC.lCode = zBC.rCode = bc;
  return create_UBspline_3d_d(x_grid, y_grid, z_grid, xBC, yBC, zBC, NULL);
}

template<typename VT, typename IT>
UBspline_3d_d* create(UBspline_3d_d* s, VT& start, VT& end, IT& ng, IT& halfg, int n = 1)
{
  Ugrid x_grid, y_grid, z_grid;
  BCtype_d xBC, yBC, zBC;
  x_grid.start = start[0];
  x_grid.end   = end[0];
  x_grid.num   = ng[0];
  y_grid.start = start[1];
  y_grid.end   = end[1];
  y_grid.num   = ng[1];
  z_grid.start = start[2];
  z_grid.end   = end[2];
  z_grid.num   = ng[2];
  xBC.lCode = xBC.rCode = (halfg[0]) ? ANTIPERIODIC : PERIODIC;
  yBC.lCode = yBC.rCode = (halfg[1]) ? ANTIPERIODIC : PERIODIC;
  zBC.lCode = zBC.rCode = (halfg[2]) ? ANTIPERIODIC : PERIODIC;
  ;
  return create_UBspline_3d_d(x_grid, y_grid, z_grid, xBC, yBC, zBC, NULL);
}

inline void set(UBspline_3d_d* s, double* restrict data) { recompute_UBspline_3d_d(s, data); }

inline void set(multi_UBspline_3d_d* spline, int i, UBspline_3d_d* spline_in, const int* offset, const int* N)
{
  copy_UBspline_3d_d(spline, i, spline_in, offset, N);
}

inline void set(multi_UBspline_3d_s* spline, int i, UBspline_3d_d* spline_in, const int* offset, const int* N)
{
  copy_UBspline_3d_d_s(spline, i, spline_in, offset, N);
}



// 1D spline
/** create a single spline for double */
template<typename VT>
UBspline_1d_d* create(UBspline_1d_d* s,
                      const VT& start,
                      const VT& end,
                      const int spline_npoints,
                      double* indata,
                      bool lFlat)
{
  BCtype_d bc;
  if (lFlat)
    bc.lCode = FLAT;
  else
    bc.lCode = NATURAL;
  bc.rCode = NATURAL;
  Ugrid grid;
  grid.start = start;
  grid.end   = end;
  grid.num   = spline_npoints;
  return create_UBspline_1d_d(grid, bc, indata);
}



/** spline destroy */
template<typename SplineType>
inline void destroy(SplineType* restrict& spline)
{
  if (spline != NULL)
  {
    if (spline->coefs != NULL)
      free(spline->coefs);
    free(spline);
    spline = NULL;
  }
}

// 1D spline interface to einspline routines.
inline void set(multi_UBspline_1d_d* spline, int i, UBspline_1d_d* spline_in, const int offset, const int N)
{
  copy_UBspline_1d_d(spline, i, spline_in, offset, N);
}

inline void set(multi_UBspline_1d_s* spline, int i, UBspline_1d_d* spline_in, const int offset, const int N)
{
  copy_UBspline_1d_d_s(spline, i, spline_in, offset, N);
}

inline multi_UBspline_1d_d* create(multi_UBspline_1d_d* s, Ugrid& grid, BCtype_d& bc, int num_splines)
{
  multi_UBspline_1d_d* newspline = create_multi_UBspline_1d_d(grid, bc, num_splines);
  std::fill(newspline->coefs, newspline->coefs + newspline->coefs_size, 0.0);
  return newspline;
}

inline multi_UBspline_1d_s* create(multi_UBspline_1d_s* s, Ugrid& grid, BCtype_s& bc, int num_splines)
{
  multi_UBspline_1d_s* newspline = create_multi_UBspline_1d_s(grid, bc, num_splines);
  std::fill(newspline->coefs, newspline->coefs + newspline->coefs_size, 0.0f);
  return newspline;
}

} // namespace einspline
} // namespace qmcplusplus
#endif

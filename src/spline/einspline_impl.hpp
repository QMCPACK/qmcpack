//////////////////////////////////////////////////////////////////
// (c) Copyright 2006-  by Jeongnim Kim and Ken Esler           //
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
/** @file einspline_impl.hpp
 * @brief define free functions in qmcplusplus::einspline namespace to handle einspline calls transparently
 *
 * The template parameters of the functions  are
 * \tparam PT position type, e.g. TinyVector<T,D>
 * \tparam VT array of values, e.g. Vector<T>
 * \tparam GT array of gradients, e.g. Vector<TinyVector<T,D> >
 * \tparam HT array of hessian tensors, e.g. Vector<Tensor<T,D> >
 */
#ifndef QMCPLUSPLUS_EINSPLINE_IMPL_H
#define QMCPLUSPLUS_EINSPLINE_IMPL_H

#include <bspline/bspline_traits.hpp>
namespace qmcplusplus
{
  namespace einspline
  {
    /** evaluate values only using multi_UBspline_3d_d 
     */
    template<typename PT, typename VT, typename GT>
      inline void evaluate(multi_UBspline_3d_d *restrict spline, const PT& r, VT &psi)
      { eval_multi_UBspline_3d_d (spline, r[0], r[1], r[2], psi.data()); }

    /** evaluate values and gradients using multi_UBspline_3d_d 
     */
    template<typename PT, typename VT, typename GT>
      inline void evaluate(multi_UBspline_3d_d *restrict spline, const PT& r, VT &psi, GT &grad)
      { eval_multi_UBspline_3d_d_vg (spline, r[0], r[1], r[2], psi.data(), grad[0].data()); }

    /** evaluate values, gradients and laplacians using multi_UBspline_3d_d 
     */
    template<typename PT, typename VT, typename GT>
      inline void evaluate(multi_UBspline_3d_d *restrict spline, const PT& r, VT &psi, GT &grad, VT& lap)
      { eval_multi_UBspline_3d_d_vgl (spline, r[0], r[1], r[2], psi.data(), grad[0].data(), lap.data()); }

    /** evaluate values, gradients and hessians using multi_UBspline_3d_d 
     */
    template<typename PT, typename VT, typename GT, typename HT>
      inline void evaluate(multi_UBspline_3d_d *restrict spline, const PT& r, VT &psi, GT &grad, HT& hess)
      { eval_multi_UBspline_3d_d_vgh (spline, r[0], r[1], r[2], psi.data(), grad[0].data(),hess[0].data()); }

    /** evaluate values only using multi_UBspline_3d_z 
     */
    template<typename PT, typename VT, typename GT>
      inline void evaluate(multi_UBspline_3d_z *restrict spline, const PT& r, VT &psi)
      { eval_multi_UBspline_3d_z (spline, r[0], r[1], r[2], psi.data()); }

    /** evaluate values and gradients using multi_UBspline_3d_z 
     */
    template<typename PT, typename VT, typename GT>
      inline void evaluate(multi_UBspline_3d_z *restrict spline, const PT& r, VT &psi, GT &grad)
      { eval_multi_UBspline_3d_z_vg (spline, r[0], r[1], r[2], psi.data(), grad[0].data()); }

    /** evaluate values, gradients and laplacians using multi_UBspline_3d_z 
     */
    template<typename PT, typename VT, typename GT>
      inline void evaluate(multi_UBspline_3d_z *restrict spline, const PT& r, VT &psi, GT &grad, VT& lap)
      { eval_multi_UBspline_3d_z_vgl (spline, r[0], r[1], r[2], psi.data(), grad[0].data(), lap.data()); }

    /** evaluate values, gradients and hessians using multi_UBspline_3d_z 
     */
    template<typename PT, typename VT, typename GT, typename HT>
      inline void evaluate(multi_UBspline_3d_z *restrict spline, const PT& r, VT &psi, GT &grad, HT& hess)
      { eval_multi_UBspline_3d_z_vgh (spline, r[0], r[1], r[2], psi.data(), grad[0].data(),hess[0].data());}

    /** evaluate values only using multi_UBspline_3d_s 
     */
    template<typename PT, typename VT, typename GT>
      inline void evaluate(multi_UBspline_3d_s *restrict spline, const PT& r, VT &psi)
      { eval_multi_UBspline_3d_s (spline, r[0], r[1], r[2], psi.data()); }

    /** evaluate values and gradients using multi_UBspline_3d_s 
     */
    template<typename PT, typename VT, typename GT>
      inline void evaluate(multi_UBspline_3d_s *restrict spline, const PT& r, VT &psi, GT &grad)
      { eval_multi_UBspline_3d_s_vg (spline, r[0], r[1], r[2], psi.data(), grad[0].data()); }

    /** evaluate values, gradients and laplacians using multi_UBspline_3d_s 
     */
    template<typename PT, typename VT, typename GT>
      inline void evaluate(multi_UBspline_3d_s *restrict spline, const PT& r, VT &psi, GT &grad, VT& lap)
      { eval_multi_UBspline_3d_s_vgl (spline, r[0], r[1], r[2], psi.data(), grad[0].data(), lap.data()); }

    /** evaluate values, gradients and hessians using multi_UBspline_3d_s 
     */
    template<typename PT, typename VT, typename GT, typename HT>
      inline void evaluate(multi_UBspline_3d_s *restrict spline, const PT& r, VT &psi, GT &grad, HT& hess)
      { eval_multi_UBspline_3d_s_vgh (spline, r[0], r[1], r[2], psi.data(), grad[0].data(),hess[0].data()); }

    /** evaluate values only using multi_UBspline_3d_c 
     */
    template<typename PT, typename VT, typename GT>
      inline void evaluate(multi_UBspline_3d_c *restrict spline, const PT& r, VT &psi)
      { eval_multi_UBspline_3d_c (spline, r[0], r[1], r[2], psi.data()); }

    /** evaluate values and gradients using multi_UBspline_3d_c 
     */
    template<typename PT, typename VT, typename GT>
      inline void evaluate(multi_UBspline_3d_c *restrict spline, const PT& r, VT &psi, GT &grad)
      { eval_multi_UBspline_3d_c_vg (spline, r[0], r[1], r[2], psi.data(), grad[0].data()); }

    /** evaluate values, gradients and laplacians using multi_UBspline_3d_c 
     */
    template<typename PT, typename VT, typename GT>
      inline void evaluate(multi_UBspline_3d_c *restrict spline, const PT& r, VT &psi, GT &grad, VT& lap)
      { eval_multi_UBspline_3d_c_vgl (spline, r[0], r[1], r[2], psi.data(), grad[0].data(), lap.data()); }

    /** evaluate values, gradients and hessians using multi_UBspline_3d_c 
     */
    template<typename PT, typename VT, typename GT, typename HT>
      inline void evaluate(multi_UBspline_3d_c *restrict spline, const PT& r, VT &psi, GT &grad, HT& hess)
      { eval_multi_UBspline_3d_c_vgh (spline, r[0], r[1], r[2], psi.data(), grad[0].data(),hess[0].data());}
  }
}
#endif

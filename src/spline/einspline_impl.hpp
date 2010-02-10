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
 * @brief define einspline class with static functions
 *
 */
#ifndef QMCPLUSPLUS_EINSPLINE_IMPL_H
#define QMCPLUSPLUS_EINSPLINE_IMPL_H

#ifndef QMCPLUSPLUS_EINSPLINE_ENGINE_HPP
#error "einspline_impl.hpp is used only by einspline_eingine.hpp"
#else
namespace qmcplusplus
{
  /** class to handle functions
   *
   * static functions to handle einspline calls transparently
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
  struct einspline
  {
    /** create  multi_UBspline_3d_d*
     * @param s dummy multi_UBspline_3d_d* 
     * @param start starting grid values
     * @param end ending grid values
     * @param ng number of grids for [start,end)
     * @param bc boundary condition
     * @param num_splines number of splines to do 
     */
    template<typename VT, typename IT>
    static multi_UBspline_3d_d*  create(multi_UBspline_3d_d* s
        , VT& start , VT& end, IT& ng , bc_code bc, int num_splines)
    { 
      Ugrid x_grid, y_grid, z_grid;
      BCtype_d xBC,yBC,zBC;
      x_grid.start=start[0]; x_grid.end=end[0]; x_grid.num=ng[0];
      y_grid.start=start[1]; y_grid.end=end[1]; y_grid.num=ng[1];
      z_grid.start=start[2]; z_grid.end=end[2]; z_grid.num=ng[2];
      xBC.lCode=xBC.rCode=bc;
      yBC.lCode=yBC.rCode=bc;
      zBC.lCode=zBC.rCode=bc;
      return create_multi_UBspline_3d_d(x_grid,y_grid,z_grid, xBC, yBC, zBC, num_splines);
    }

    template<typename VAT>                                                                 
      static inline void  set(multi_UBspline_3d_d* spline, int i, VAT& indata)
      { set_multi_UBspline_3d_d(spline, i, indata.data()); }                                                            

    /** evaluate values only using multi_UBspline_3d_d 
    */
    template<typename PT, typename VT>
      static inline void  evaluate(multi_UBspline_3d_d *restrict spline, const PT& r, VT &psi)
      { eval_multi_UBspline_3d_d (spline, r[0], r[1], r[2], psi.data()); }

    /** evaluate values and gradients using multi_UBspline_3d_d 
    */
    template<typename PT, typename VT, typename GT>
      static inline void  evaluate(multi_UBspline_3d_d *restrict spline, const PT& r, VT &psi, GT &grad)
      { eval_multi_UBspline_3d_d_vg (spline, r[0], r[1], r[2], psi.data(), grad[0].data()); }

    /** evaluate values, gradients and laplacians using multi_UBspline_3d_d 
    */
    template<typename PT, typename VT, typename GT>
      static inline void  evaluate(multi_UBspline_3d_d *restrict spline, const PT& r, VT &psi, GT &grad, VT& lap)
      { eval_multi_UBspline_3d_d_vgl (spline, r[0], r[1], r[2], psi.data(), grad[0].data(), lap.data()); }

    /** evaluate values, gradients and hessians using multi_UBspline_3d_d 
    */
    template<typename PT, typename VT, typename GT, typename HT>
      static inline void  evaluate(multi_UBspline_3d_d *restrict spline, const PT& r, VT &psi, GT &grad, HT& hess)
      { eval_multi_UBspline_3d_d_vgh (spline, r[0], r[1], r[2], psi.data(), grad[0].data(),hess[0].data()); }

    /** create spline for complex<double> */
    template<typename VT, typename IT>
    static multi_UBspline_3d_z*  create(multi_UBspline_3d_z* s
        , VT& start , VT& end, IT& ng , bc_code bc, int num_splines)
    { 
      Ugrid x_grid, y_grid, z_grid;
      BCtype_z xBC,yBC,zBC;
      x_grid.start=start[0]; x_grid.end=end[0]; x_grid.num=ng[0];
      y_grid.start=start[1]; y_grid.end=end[1]; y_grid.num=ng[1];
      z_grid.start=start[2]; z_grid.end=end[2]; z_grid.num=ng[2];
      xBC.lCode=xBC.rCode=bc;
      yBC.lCode=yBC.rCode=bc;
      zBC.lCode=zBC.rCode=bc;
      return create_multi_UBspline_3d_z(x_grid,y_grid,z_grid, xBC, yBC, zBC, num_splines);
    }

    template<typename VAT>                                                                 
      static inline void  set(multi_UBspline_3d_z* spline, int i, VAT& indata)
      { set_multi_UBspline_3d_z(spline, i, indata.data()); }                                                            

    /** evaluate values only using multi_UBspline_3d_z 
    */
    template<typename PT, typename VT>
      static inline void  evaluate(multi_UBspline_3d_z *restrict spline, const PT& r, VT &psi)
      { eval_multi_UBspline_3d_z (spline, r[0], r[1], r[2], psi.data()); }

    /** evaluate values and gradients using multi_UBspline_3d_z 
    */
    template<typename PT, typename VT, typename GT>
      static inline void  evaluate(multi_UBspline_3d_z *restrict spline, const PT& r, VT &psi, GT &grad)
      { eval_multi_UBspline_3d_z_vg (spline, r[0], r[1], r[2], psi.data(), grad[0].data()); }

    /** evaluate values, gradients and laplacians using multi_UBspline_3d_z 
    */
    template<typename PT, typename VT, typename GT>
      static inline void  evaluate(multi_UBspline_3d_z *restrict spline, const PT& r, VT &psi, GT &grad, VT& lap)
      { eval_multi_UBspline_3d_z_vgl (spline, r[0], r[1], r[2], psi.data(), grad[0].data(), lap.data()); }

    /** evaluate values, gradients and hessians using multi_UBspline_3d_z 
    */
    template<typename PT, typename VT, typename GT, typename HT>
      static inline void  evaluate(multi_UBspline_3d_z *restrict spline, const PT& r, VT &psi, GT &grad, HT& hess)
      { eval_multi_UBspline_3d_z_vgh (spline, r[0], r[1], r[2], psi.data(), grad[0].data(),hess[0].data());}

    /** create spline and initialized it */
    template<typename VT, typename IT>
    static multi_UBspline_3d_s*  create(multi_UBspline_3d_s* s
        , VT& start , VT& end, IT& ng , bc_code bc, int num_splines)
    { 
      Ugrid x_grid, y_grid, z_grid;
      BCtype_s xBC,yBC,zBC;
      x_grid.start=start[0]; x_grid.end=end[0]; x_grid.num=ng[0];
      y_grid.start=start[1]; y_grid.end=end[1]; y_grid.num=ng[1];
      z_grid.start=start[2]; z_grid.end=end[2]; z_grid.num=ng[2];
      xBC.lCode=xBC.rCode=bc;
      yBC.lCode=yBC.rCode=bc;
      zBC.lCode=zBC.rCode=bc;
      return create_multi_UBspline_3d_s(x_grid,y_grid,z_grid, xBC, yBC, zBC, num_splines);
    }

    template<typename VAT>                                                                 
      static inline void  set(multi_UBspline_3d_s* spline, int i, VAT& indata)
      { set_multi_UBspline_3d_s(spline, i, indata.data()); }                                                            

    /** evaluate values only using multi_UBspline_3d_s 
    */
    template<typename PT, typename VT>
      static inline void  evaluate(multi_UBspline_3d_s *restrict spline, const PT& r, VT &psi)
      { eval_multi_UBspline_3d_s (spline, r[0], r[1], r[2], psi.data()); }

    /** evaluate values and gradients using multi_UBspline_3d_s 
    */
    template<typename PT, typename VT, typename GT>
      static inline void  evaluate(multi_UBspline_3d_s *restrict spline, const PT& r, VT &psi, GT &grad)
      { eval_multi_UBspline_3d_s_vg (spline, r[0], r[1], r[2], psi.data(), grad[0].data()); }

    /** evaluate values, gradients and laplacians using multi_UBspline_3d_s 
    */
    template<typename PT, typename VT, typename GT>
      static inline void  evaluate(multi_UBspline_3d_s *restrict spline, const PT& r, VT &psi, GT &grad, VT& lap)
      { eval_multi_UBspline_3d_s_vgl (spline, r[0], r[1], r[2], psi.data(), grad[0].data(), lap.data()); }

    /** evaluate values, gradients and hessians using multi_UBspline_3d_s 
    */
    template<typename PT, typename VT, typename GT, typename HT>
      static inline void  evaluate(multi_UBspline_3d_s *restrict spline, const PT& r, VT &psi, GT &grad, HT& hess)
      { eval_multi_UBspline_3d_s_vgh (spline, r[0], r[1], r[2], psi.data(), grad[0].data(),hess[0].data()); }

    /** create spline for complex<double> */
    template<typename VT, typename IT>
    static multi_UBspline_3d_c*  create(multi_UBspline_3d_c* s
        , VT& start , VT& end, IT& ng , bc_code bc, int num_splines)
    { 
      Ugrid x_grid, y_grid, z_grid;
      BCtype_c xBC,yBC,zBC;
      x_grid.start=start[0]; x_grid.end=end[0]; x_grid.num=ng[0];
      y_grid.start=start[1]; y_grid.end=end[1]; y_grid.num=ng[1];
      z_grid.start=start[2]; z_grid.end=end[2]; z_grid.num=ng[2];
      xBC.lCode=xBC.rCode=bc;
      yBC.lCode=yBC.rCode=bc;
      zBC.lCode=zBC.rCode=bc;
      return create_multi_UBspline_3d_c(x_grid,y_grid,z_grid, xBC, yBC, zBC, num_splines);
    }

    template<typename VAT>                                                                 
      static inline void  set(multi_UBspline_3d_c* spline, int i, VAT& indata)
      { set_multi_UBspline_3d_c(spline, i, indata.data()); }                                                            

    /** evaluate values only using multi_UBspline_3d_c 
    */
    template<typename PT, typename VT>
      static inline void  evaluate(multi_UBspline_3d_c *restrict spline, const PT& r, VT &psi)
      { eval_multi_UBspline_3d_c (spline, r[0], r[1], r[2], psi.data()); }

    /** evaluate values and gradients using multi_UBspline_3d_c 
    */
    template<typename PT, typename VT, typename GT>
      static inline void  evaluate(multi_UBspline_3d_c *restrict spline, const PT& r, VT &psi, GT &grad)
      { eval_multi_UBspline_3d_c_vg (spline, r[0], r[1], r[2], psi.data(), grad[0].data()); }

    /** evaluate values, gradients and laplacians using multi_UBspline_3d_c 
    */
    template<typename PT, typename VT, typename GT>
      static inline void  evaluate(multi_UBspline_3d_c *restrict spline, const PT& r, VT &psi, GT &grad, VT& lap)
      { eval_multi_UBspline_3d_c_vgl (spline, r[0], r[1], r[2], psi.data(), grad[0].data(), lap.data()); }

    /** evaluate values, gradients and hessians using multi_UBspline_3d_c 
    */
    template<typename PT, typename VT, typename GT, typename HT>
      static inline void  evaluate(multi_UBspline_3d_c *restrict spline, const PT& r, VT &psi, GT &grad, HT& hess)
      { eval_multi_UBspline_3d_c_vgh (spline, r[0], r[1], r[2], psi.data(), grad[0].data(),hess[0].data());}
  };
}
#endif
#endif

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
    
    



#ifndef QMCPLUSPLUS_CREATE_GRID_FUNCTION_H
#define QMCPLUSPLUS_CREATE_GRID_FUNCTION_H

#include <cmath>
/**Abstract class to perform a transformation from a function (analytic or grid)
 * to a grid function
 */
template<class FnOut>
struct Transform2GridFunctorBase
{
  typedef typename FnOut::point_type point_type;
  virtual void generate(point_type ri, point_type rf, int ng, int np=0) = 0;
  virtual void generate(int np=0) = 0;
};

/**Class to transform a function (analytic or grid) to a new grid function.
 *  For an input function \f$ F_{in} \f$
 * \f[ \frac{F_{in}(r)}{r^{np}} \longrightarrow F'_{out}(r'), \f]
 * where \f$ r' \f$ is the new grid and \f$ np \f$ is an integer
 * power.
 */
template<class FnIn, class FnOut>
struct Transform2GridFunctor:
  public Transform2GridFunctorBase<FnOut>
{

  typedef typename FnIn::value_type  result_t;
  typedef typename FnOut::point_type point_type;

  FnIn& in_;
  FnOut& out_;

  /*! \fn Transform2GridFunctor(FnIn& in, FnOut& out)
   * \param in the input function
   * \param out returns the output function
   * \brief The constructor.
  */

  Transform2GridFunctor(FnIn& in, FnOut& out):in_(in), out_(out) {}

  /*! \fn void generate(point_type ri, point_type rf, int ng, int np=0)
   * \param ri the initial grid point
   * \param rf the final grid point
   * \param ng number of grid points
   * \param np the inverse power
   * \brief Initializes the grid of the output function, assigns
   values to the grid points and performs the 1D-Cubic Spline
   for interpolation.
  */
  void generate(point_type ri, point_type rf, int ng, int np=0)
  {
    //reference to the output functions grid
    typename FnOut::grid_type& grid = out_.grid();
    //set the output functions grid
    grid.set(ri,rf,ng);
    out_.resize(ng);
    //assign values to the output function
    if(np == 0)
    {
      for(int i=0; i<grid.size(); i++)
      {
        out_(i) = in_.f(grid(i));
      }
    }
    else
    {
      for(int i=0; i<grid.size(); i++)
      {
        out_(i) = in_.f(grid(i))*std::pow(grid(i),np);
      }
    }
    //boundary conditions
    result_t deriv= in_.df(grid(0));
    if(np>0)
    {
      point_type r = grid(0);
      deriv = (deriv*pow(r,np)+static_cast<result_t>(np)*out_(0)*pow(r,np-1));
    }
    //spline the output function
    out_.spline(0,deriv,grid.size()-1,0.0);
  }

  inline void generate(point_type rf)
  {
    //reference to the output functions grid
    typename FnOut::grid_type& grid = out_.grid();
    grid.locate(rf);
    int npts(grid.currentIndex()+1);
    //int npts = grid.index(rf)+1;
    out_.resize(grid.size());
    for(int i=0; i<npts; i++)
      out_(i) = in_.f(grid(i));
    //boundary conditions
    result_t deriv= in_.df(grid(0));
    //spline the output function
    out_.spline(0,deriv,npts-1,0.0);
  }


  /*! \fn void generate(int np=0)
    * \param np the inverse power
    * \brief Assigns values to the grid points and performs the 1D-Cubic Spline
    for interpolation.  Assumes that the grid of the output function has
    already been initialized.
   */

  void generate(int np=0)
  {
    typename FnOut::grid_type& grid = out_.grid();
    out_.resize(grid.size());
    if(np == 0)
    {
      for(int i=0; i<grid.size(); i++)
      {
        out_(i) = in_.f(grid(i));
      }
    }
    else
    {
      for(int i=0; i<grid.size(); i++)
      {
        out_(i) = in_.f(grid(i))*std::pow(grid(i),np);
      }
    }
    result_t deriv= in_.df(grid(0));
    if(np>0)
    {
      point_type r = grid(0);
      deriv = (deriv*std::pow(r,np)+static_cast<result_t>(np)*out_(0)*std::pow(r,np-1));
    }
    out_.spline(0,deriv,grid.size()-1,0.0);
  }
};

/**create a transform engine with a given input function **/
template<class FnIn, class FnOut>
Transform2GridFunctorBase<FnOut>*
createTransform2GridFunctor(FnIn& in, FnOut& out)
{
  return new Transform2GridFunctor<FnIn,FnOut>(in,out);
}


template<class FnIn, class FnOut>
struct TestTransform
{

  typedef typename FnOut::value_type value_type;
  typedef typename FnOut::point_type point_type;

  static void check(FnIn& af, FnOut& nf, point_type ri, point_type rf, int n)
  {
    value_type yp0, ypn;
    value_type ypp0, yppn;
    value_type y,yp,y2p;
    value_type z,zp,z2p;
    ///create a transform engine to convert a function, e.g., RadialSTO to a grid functor
    Transform2GridFunctorBase<FnOut> *transformer
    = createTransform2GridFunctor(af,nf);
    /// In = out/r^np
    transformer->generate(ri,rf,n);
    for(int i=1; i<nf.size()-1; i++)
    {
      value_type x0 = nf.r(i);
      nf.setgrid(x0);
      y = nf.evaluate(x0,1.0/x0,yp,y2p);
      z = af.evaluate(x0,1.0/x0,zp,z2p);
      std::cout << std::setw(10) << x0
                << std::setw(20) << std::setprecision(10) << y
                << std::setw(20) << (y-z)/z
                << std::setw(20) << yp
                << std::setw(20) << (yp-zp)/zp
                << std::setw(20) << y2p
                << std::setw(20) << (y2p-z2p)/z2p << std::endl;
    }
    delete transformer;
  }

};


#endif

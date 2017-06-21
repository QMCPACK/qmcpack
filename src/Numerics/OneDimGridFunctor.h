//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



#ifndef QMCPLUSPLUS_GRID_FUNCTOR_H
#define QMCPLUSPLUS_GRID_FUNCTOR_H

#include "Numerics/OneDimGridBase.h"
#include "Optimize/VarList.h"

namespace qmcplusplus
{
template<class T, unsigned D>
struct FunctorBase { };

/** Implement One-Dimensional function on a radial grid.
 *
 * template parameters
 * - Td return type
 * - Tg grid value type
 * - CTd container type associated with the values and derivatives
 * - CTg container type associated with the grid data
 *
 * Store the values of the function for the
 * cooresponding grid points, \f$ y_i = y(x_i) \f$.
 */
template <class Td,
         class Tg = Td,
         class CTd= Vector<Td>,
         class CTg= Vector<Tg>  >
struct OneDimGridFunctor//: public FunctorBase<Td,1> {
{

  /// the type of the value on a grid
  typedef Td  value_type;
  /// the type of the grid value
  typedef Tg  point_type;
  /// the type of the containers Y, dY and d2Y
  typedef CTd data_type;
  /// the grid type
  typedef OneDimGridBase<Tg,CTg> grid_type;
  /// the type of this class
  typedef OneDimGridFunctor<Td,Tg,CTd,CTg>  this_type;

  /** constructor
   *@param gt a radial grid
   */
  OneDimGridFunctor(grid_type* gt = 0)
    : GridManager(true), OwnGrid(false), m_grid(gt)
  {
    if(m_grid)
      resize(m_grid->size());
    //FirstAddress.resize(3,0);
  }

  /** virtual destructor */
  inline virtual ~OneDimGridFunctor()
  {
    if(OwnGrid&&m_grid)
      delete m_grid;
  }

  /////copy constructor
  //OneDimGridFunctor(const this_type& a): GridManager(true), m_grid(a.m_grid){
  //  if(m_grid) resize(m_grid->size());
  //  FirstAddress.resize(3,0);
  //}


  OneDimGridFunctor<Td,Tg,CTd,CTg>(const OneDimGridFunctor<Td,Tg,CTd,CTg>& a)
  {
    GridManager = a.GridManager;
    OwnGrid=true;
    m_grid = a.m_grid->makeClone();
    Y = a.Y;
    dY = a.dY;
    d2Y = a.d2Y;
    m_Y.resize(a.m_Y.size());
    m_Y = a.m_Y;
    NumNodes = a.NumNodes;
  }

  /// assignment operator
  const this_type& operator=(const this_type& a)
  {
    //This object does not manage the grid
    GridManager=false;
    OwnGrid=false;
    m_grid = a.m_grid;
    m_Y = a.m_Y;
    //m_Y2 = a.m_Y2;
    return *this;
  }

  template<class T1>
  const this_type& operator=(const T1& x)
  {
    Y = x;
    return *this;
  }

  template<typename TT>
  inline void resetParameters(const TT& active)
  {
  }

  ///set the number of nodes
  inline void setNumOfNodes(int n)
  {
    NumNodes = n;
  }

  ///return the number of nodes
  inline int getNumOfNodes() const
  {
    return NumNodes;
  }

  ///return the grid data
  inline value_type* data()
  {
    return &(m_Y[0]);
  }
  ///assign the grid data
  inline const value_type* data() const
  {
    return &(m_Y[0]);
  }
  ///return the number of data points
  inline int size() const
  {
    return m_Y.size();
  }
  ///resize the number of data points
  inline void resize(int n)
  {
    m_Y.resize(n);
  }
  ///return the radial grid
  inline const grid_type& grid() const
  {
    return *m_grid;
  }
  ///assign a radial grid
  inline grid_type& grid()
  {
    return *m_grid;
  }
  ///set the status of GridManager
  inline void setGridManager(bool willmanage)
  {
    GridManager = willmanage;
  }

  /**returns a value
   * @param i grid index
   * @return the value at i
   */
  inline value_type operator()(int i) const
  {
    return m_Y[i];
  }

  /**asign a value at i
   * @param i grid index
   * @return the value at i
   */
  inline value_type& operator()(int i)
  {
    return m_Y[i];
  }

  ///** return the address of the values
  // * @param i index, i=0 value, i=1 first derivative, i=2 second
  // */
  //inline value_type* data(int i) {
  //  return 0;
  //  //return FirstAddress[i];
  //}

  /**return the differntial spacing for the grid
   *@warning only for LinearGrid and LogGrid
  */
  inline point_type dh() const
  {
    return m_grid->dh();
  }

  ///return \f$r(i)\f$ the grid point at index i
  inline point_type r(int i) const
  {
    return m_grid->r(i);
  }
  ///return \f$r(i+1)-r(i)\f$
  inline point_type dr(int i) const
  {
    return m_grid->dr(i);
  }

  /** Evaluate the function and its derivatives, store the derivatives.
   *@param r radial distance
   *@return the value of the function
   */
  inline value_type f(point_type r)
  {
    //setgrid(r);
    return Y=splint(r);
  }

  /** Evaluate the function and its derivatives, store the derivatives.
   *@param r radial distance
   *@return the derivative of the function
   */
  inline value_type df(point_type r)
  {
    //setgrid(r);
    Y=splint(r,dY,d2Y);
    return dY;
  }

  /** Evaluate the function value only
   * @param r value on a grid
   * @param rinv inverse of r
   * @return value at r
   */
  inline value_type evaluate(point_type r, point_type rinv)
  {
    return Y = splint(r);
  }

  /** Evaluate the function and its derivatives.
   * @param r value on a grid
   * @param rinv inverse of r
   * @return value at r
   *
   * Derivatives are storged.
   */
  inline value_type evaluateAll(point_type r, point_type rinv)
  {
    return Y = splint(r,dY,d2Y);
  }

  /////reset the values: do nothing
  //virtual void reset() { }
  /////reset the values from the pool
  //virtual void resetParameters(VarRegistry<point_type>& vlist)
  //{
  //  ///DO NOTHING
  //}

  virtual
  value_type
  splint(point_type r, value_type& du, value_type& d2u)
  {
    return 0.0;
  }

  virtual
  value_type splint(point_type r)
  {
    return 0.0;
  }

  virtual
  void spline(int imin, value_type yp1, int imax, value_type ypn) {  }

  virtual void spline() {  }

  /**
   *@param r radial distance
   *@param rinv inverse of radial distance
   *@param du return derivative
   *@param d2u return 2nd derivative
   *@return the value of the function
   *@brief Evaluate the function and its derivatives.
   */
  inline value_type
  evaluate(point_type r, point_type rinv, value_type& du, value_type& d2u)
  {
    return splint(r,du,d2u);
  }

  ///true, if this object manages the grid
  bool GridManager;
  ///true, if owns the grid to clean up
  bool OwnGrid;
  ///pointer to the radial grid
  grid_type* m_grid;

  ///store the value of the function
  value_type Y;
  ///store the derivative of the function
  value_type dY;
  ///store the second derivative of the function
  value_type d2Y;

  ///data for the function on the grid
  data_type m_Y;

  ///the number of nodes
  int NumNodes;

  ///address of coefficients Y and dY or d2Y
  //std::vector<value_type*> FirstAddress;
};

/** One-dimensional grid functor that returns a constant
 */
template <class Td,
         class Tg = Td,
         class CTd= Vector<Td>,
         class CTg= Vector<Tg> >
class OneDimConstFunctor: public OneDimGridFunctor<Td,Tg,CTd,CTg>
{

public:

  Td ConstValue;
  typedef OneDimGridFunctor<Td,Tg,CTd,CTg> base_type;
  typedef typename base_type::value_type  value_type;
  typedef typename base_type::point_type  point_type;
  typedef typename base_type::data_type data_type;
  typedef typename base_type::grid_type grid_type;


  OneDimConstFunctor(grid_type* gt = 0): base_type(gt), ConstValue(0.0) { }

  inline value_type splint(point_type r)
  {
    return ConstValue;
  }

  inline value_type
  splint(point_type r, value_type& du, value_type& d2u)
  {
    du=0.0;
    d2u=0.0;
    return ConstValue;
  }

  inline
  void spline(int imin, value_type yp1, int imax, value_type ypn)
  {
  }

  inline void spline() { }
};

}
#endif

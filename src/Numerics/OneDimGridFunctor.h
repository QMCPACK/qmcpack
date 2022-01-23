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

#include <memory>
#include "Numerics/OneDimGridBase.h"

namespace qmcplusplus
{
/** Implement One-Dimensional function on a radial grid.
 *
 * template parameters
 * - Td return type
 * - Tg grid value type
 * - CTd container type associated with the values and derivatives
 * - CTg container type associated with the grid data
 *
 * Store the values of the function for the
 * corresponding grid points, \f$ y_i = y(x_i) \f$.
 */
template<class Td, class Tg = Td, class CTd = Vector<Td>, class CTg = Vector<Tg>>
struct OneDimGridFunctor
{
  /// the type of the value on a grid
  using value_type = Td;
  /// the type of the grid value
  using point_type = Tg;
  /// the type of the containers Y, dY and d2Y
  using data_type = CTd;
  /// the grid type
  using grid_type = OneDimGridBase<Tg, CTg>;
  /// the type of this class
  using this_type = OneDimGridFunctor<Td, Tg, CTd, CTg>;

  /** constructor
   *@param gt a radial grid. The pointer is treated as a reference
   */
  OneDimGridFunctor(std::unique_ptr<grid_type> gt = std::unique_ptr<grid_type>()) : m_grid(std::move(gt))
  {
    if (m_grid)
      resize(m_grid->size());
  }

  OneDimGridFunctor(const OneDimGridFunctor& a)
  {
    if (a.m_grid)
      m_grid = a.m_grid->makeClone();
    Y   = a.Y;
    dY  = a.dY;
    d2Y = a.d2Y;
    m_Y.resize(a.m_Y.size());
    m_Y      = a.m_Y;
    NumNodes = a.NumNodes;
  }

  virtual ~OneDimGridFunctor() = default;

  template<typename TT>
  inline void resetParameters(const TT& active)
  {}

  ///set the number of nodes
  inline void setNumOfNodes(int n) { NumNodes = n; }

  ///return the number of nodes
  inline int getNumOfNodes() const { return NumNodes; }

  ///return the grid data
  inline value_type* data() { return &(m_Y[0]); }
  ///assign the grid data
  inline const value_type* data() const { return &(m_Y[0]); }
  ///return the number of data points
  inline int size() const { return m_Y.size(); }
  ///resize the number of data points
  inline void resize(int n) { m_Y.resize(n); }
  ///return the radial grid
  inline const grid_type& grid() const { return *m_grid; }
  ///assign a radial grid
  inline grid_type& grid() { return *m_grid; }

  /**returns a value
   * @param i grid index
   * @return the value at i
   */
  inline value_type operator()(int i) const { return m_Y[i]; }

  /**asign a value at i
   * @param i grid index
   * @return the value at i
   */
  inline value_type& operator()(int i) { return m_Y[i]; }

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
  inline point_type dh() const { return m_grid->dh(); }

  ///return \f$r(i)\f$ the grid point at index i
  inline point_type r(int i) const { return m_grid->r(i); }
  ///return \f$r(i+1)-r(i)\f$
  inline point_type dr(int i) const { return m_grid->dr(i); }

  /** Evaluate the function and its derivatives, store the derivatives.
   *@param r radial distance
   *@return the value of the function
   */
  inline value_type f(point_type r)
  {
    //setgrid(r);
    return Y = splint(r);
  }

  /** Evaluate the function and its derivatives, store the derivatives.
   *@param r radial distance
   *@return the derivative of the function
   */
  inline value_type df(point_type r)
  {
    //setgrid(r);
    Y = splint(r, dY, d2Y);
    return dY;
  }

  /** Evaluate the function value only
   * @param r value on a grid
   * @param rinv inverse of r
   * @return value at r
   */
  inline value_type evaluate(point_type r, point_type rinv) { return Y = splint(r); }

  /** Evaluate the function and its derivatives.
   * @param r value on a grid
   * @param rinv inverse of r
   * @return value at r
   *
   * Derivatives are storged.
   */
  inline value_type evaluateAll(point_type r, point_type rinv) { return Y = splint(r, dY, d2Y); }

  virtual value_type splint(point_type r, value_type& du, value_type& d2u) const { return 0.0; }

  virtual value_type splint(point_type r) const { return 0.0; }

  virtual void spline(int imin, value_type yp1, int imax, value_type ypn) {}

  virtual void spline() {}

  /**
   *@param r radial distance
   *@param rinv inverse of radial distance
   *@param du return derivative
   *@param d2u return 2nd derivative
   *@return the value of the function
   *@brief Evaluate the function and its derivatives.
   */
  inline value_type evaluate(point_type r, point_type rinv, value_type& du, value_type& d2u)
  {
    return splint(r, du, d2u);
  }

  ///pointer to the radial grid
  std::unique_ptr<grid_type> m_grid;

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
template<class Td, class Tg = Td, class CTd = Vector<Td>, class CTg = Vector<Tg>>
class OneDimConstFunctor : public OneDimGridFunctor<Td, Tg, CTd, CTg>
{
public:
  Td ConstValue;
  using base_type  = OneDimGridFunctor<Td, Tg, CTd, CTg>;
  using value_type = typename base_type::value_type;
  using point_type = typename base_type::point_type;
  using data_type  = typename base_type::data_type;
  using grid_type  = typename base_type::grid_type;


  OneDimConstFunctor(grid_type* gt = 0) : base_type(gt), ConstValue(0.0) {}

  inline value_type splint(point_type r) const override { return ConstValue; }

  inline value_type splint(point_type r, value_type& du, value_type& d2u) const override
  {
    du  = 0.0;
    d2u = 0.0;
    return ConstValue;
  }

  inline void spline(int imin, value_type yp1, int imax, value_type ypn) {}

  inline void spline() {}
};

} // namespace qmcplusplus
#endif

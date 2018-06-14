//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



#ifndef QMCPLUSPLUS_ONEDIMGRID_BASE
#define QMCPLUSPLUS_ONEDIMGRID_BASE

/**@file OneDimGridBase.h
 *@brief Decalaration of One-Dimesional grids
 */

#include <algorithm>
#include <cmath>
#include "OhmmsPETE/OhmmsVector.h"
#include "Numerics/GridTraits.h"
#include "einspline/bspline_base.h"

namespace qmcplusplus
{


/** An abstract base class to implement a One-Dimensional grid
 */
template <class T, class CT=Vector<T> >
struct OneDimGridBase
{

  typedef T value_type;
  typedef CT Array_t;


  ///the current index of the grid
  int Loc;
  int GridTag;
  int num_points;
  value_type lower_bound;
  value_type upper_bound;
  ///differential spacing of the grid
  value_type Delta;
  double     DeltaInv;
  value_type dL;
  value_type dLinv;
  value_type cL;
  value_type cR;
  value_type p1,p2,q1,q2;
  value_type dp1,dq1,dq2;
  value_type d2p1,d2q1,d2q2;

  ///array to store the radial grid data
  Array_t X;


  inline OneDimGridBase():num_points(0)
  {}

  virtual OneDimGridBase<T,CT>* makeClone() const =0;

  virtual ~OneDimGridBase() {}

  inline int getGridTag() const
  {
    return GridTag;
  }

  ///return the current index
  inline int currentIndex() const
  {
    return Loc;
  }

  inline int getIndex(T r) const
  {
    int k;
    int klo=0;
    int khi=this->size();
    while(khi-klo > 1)
    {
      k=(khi+klo) >> 1;
      if(X[k] > r)
        khi=k;
      else
        klo=k;
    }
    return klo;
  }

  ///assign a value
  inline T& operator[](int i)
  {
    return X[i];
  }
  ///assign a value
  inline T& operator()(int i)
  {
    return X[i];
  }
  ///return a value
  inline T operator[](int i) const
  {
    return X[i];
  }
  ///return a value
  inline T operator()(int i) const
  {
    return X[i];
  }

  inline const T* data() const
  {
    return &(X[0]);
  }
  inline T* data()
  {
    return &(X[0]);
  }

  ///return the differential spacing of the grid
  inline T dh() const
  {
    return Delta;
  }
  ///returns \f$ r(i)\f$
  inline T r(int i) const
  {
    return X[i];
  }
  ///returns \f$ r(i+1)-r(i)\f$
  inline T dr(int i) const
  {
    return X[i+1]-X[i];
  }
  ///returns the size of the grid
  inline int size() const
  {
    return num_points;
  }
  ///return the first grid point
  inline T rmin() const
  {
    return lower_bound;
  }
  ///return the last grid point
  inline T rmax() const
  {
    return upper_bound;
  }

  ///update the variables for interpolations
  inline void updateFirstOrder(T r, bool all)
  {
    //int khi(Loc+1);
    //h=X[khi]-X[Loc];
    ////hinv=1.0/h;
    //value_type hinv(1.0e0/h);
    //value_type t((r-X[Loc])*hinv);
    //value_type tm(t-1.0);
    //replcate t -> cL tm ->cR
    //p1=tm*tm*(1.0+2.0*t);
    //p2=t*t*(3.0-2.0*t);
    //q1=t*tm*tm;
    //q2=t*t*tm;
    //
    //dp1=6.0*t*tm*hinv;
    //dq1=(1.0-4.0*t+3.0*t*t);
    //dq2=t*(3.0*t-2.0);
    //
    //d2p1=(12.0*t-6.0)*hinv*hinv;
    //d2q1=(6.0*t-4.0)*hinv;
    //d2q2=(6.0*t-2.0)*hinv;
    //locate is called separately (PBC)
    //locate(r);
    dL = X[Loc+1]-X[Loc];
    dLinv = 1.0e0/dL;
    cL = (r-X[Loc])*dLinv;
    cR = cL-1.0;
    p1=cR*cR*(1.0+2.0*cL);
    p2=cL*cL*(3.0-2.0*cL);
    q1=cL*cR*cR;
    q2=cL*cL*cR;
    if(all)
    {
      dp1=6.0*cL*cR*dLinv;
      dq1=1.0-4.0*cL+3.0*cL*cL;
      dq2=cL*(3.0*cL-2.0);
      d2p1=(12.0*cL-6.0)*dLinv*dLinv;
      d2q1=(6.0*cL-4.0)*dLinv;
      d2q2=(6.0*cL-2.0)*dLinv;
    }
  }

  template <typename T1>
  inline T1 cubicInterpolateFirst(T1 a, T1 b, T1 a1, T1 b1)
  {
    return p1*a+p2*b+dL*(q1*a1+q2*b1);
  }

  template <typename T1>
  inline T1 cubicInterpolateFirst(T1 a, T1 b, T1 a1, T1 b1, T1& du, T1& d2u)
  {
    du = dp1*(a-b)+dq1*a1+dq2*b1;
    d2u = d2p1*(a-b)+d2q1*a1+d2q2*b1;
    return p1*a+p2*b+dL*(q1*a1+q2*b1);
  }

  inline void updateSecondOrder(T r, bool all)
  {
    //Find Loc
    locate(r);
    dL = X[Loc+1]-X[Loc];
    dLinv = 1.0e0/dL;
    cL = (r-X[Loc])*dLinv; //B
    cR = (X[Loc+1]-r)*dLinv;//A
    const T onesixth = 1.0/6.0;
    T h6(dL*onesixth);
    q1 = cR*(cR*cR-1.0)*h6*dL; //C
    q2 = cL*(cL*cL-1.0)*h6*dL; //D
    if(all)
    {
      dq1 = h6*(1.0-3.0*cR*cR);
      dq2 = h6*(3.0*cL*cL-1.0);
    }
  }

  inline void updateForQuintic(T r, bool all)
  {
    //Find Loc
    locate(r);
    cL = (r-X[Loc]);
  }

  template <typename T1>
  inline T1 cubicInterpolateSecond(T1 y1, T1 y2, T1 d2y1, T1 d2y2)
  {
    return cR*y1+cL*y2+q1*d2y1+q2*d2y2;
  }

  template <typename T1>
  inline T1 cubicInterpolateSecond(T1 y1, T1 y2, T1 d2y1, T1 d2y2,
                                   T& du, T& d2u)
  {
    du = dLinv*(y2-y1)+dq1*d2y1+dq2*d2y2;
    d2u = cR*d2y1+cL*d2y2;
    return cR*y1+cL*y2+q1*d2y1+q2*d2y2;
  }

  template <typename T1>
  inline T1 quinticInterpolate(T1 a, T1 b, T1 c, T1 d, T1 e, T1 f)
  {
    return a+cL*(b+cL*(c+cL*(d+cL*(e+cL*f))));
  }

  template <typename T1>
  inline T1 quinticInterpolate(T1 a, T1 b, T1 c, T1 d, T1 e, T1 f, T1 &du, T1 &d2u)
  {
    du = b+cL*(2.0*c+cL*(3.0*d+cL*(4.0*e+cL*f*5.0)));
    d2u = 2.0*c+cL*(6.0*d+cL*(12.0*e+cL*f*20.0));
    return a+cL*(b+cL*(c+cL*(d+cL*(e+cL*f))));
  }

  template <typename T1>
  inline T1 quinticInterpolate(T1 a, T1 b, T1 c, T1 d, T1 e, T1 f, T1 &du, T1 &d2u, T1 &d3u)
  {
    du = b+cL*(2.0*c+cL*(3.0*d+cL*(4.0*e+cL*f*5.0)));
    d2u = 2.0*c+cL*(6.0*d+cL*(12.0*e+cL*f*20.0));
    d3u = 6.0*d+cL*(24.0*e+cL*f*60.0);
    return a+cL*(b+cL*(c+cL*(d+cL*(e+cL*f))));
  }

  /** evaluate the index of r
   * @param r current position
   *
   * The grid index satifies \f$ X[Loc] \ge r < X[Loc+1]\f$.
   */
  virtual void locate(T r)=0;

  /** Set the grid given the parameters.
   *@param ri initial grid point
   *@param rf final grid point
   *@param n number of grid points
   */
  virtual void set(T ri, T rf, int n)=0;

};

/** One-Dimensional linear-grid.
 *
 * The analytic form \f$ r_i = r_0 +
 * i\left( \frac{r_f - r_0}{N-1} \right), \f$
 * where \f$ N \f$ is the number of points and the index
 * \f$ i \f$ runs from 0 to \f$ N-1 \f$.
 */
template <class T, class CT=Vector<T> >
struct LinearGrid: public OneDimGridBase<T,CT>
{

  using OneDimGridBase<T,CT>::Loc;
  using OneDimGridBase<T,CT>::GridTag;
  using OneDimGridBase<T,CT>::num_points;
  using OneDimGridBase<T,CT>::lower_bound;
  using OneDimGridBase<T,CT>::upper_bound;
  using OneDimGridBase<T,CT>::Delta;
  using OneDimGridBase<T,CT>::DeltaInv;
  using OneDimGridBase<T,CT>::X;

  OneDimGridBase<T,CT>* makeClone() const
  {
    return new LinearGrid<T,CT>(*this);
  }

  inline void locate(T r)
  {
    Loc = static_cast<int>((static_cast<double>(r)-X[0])*DeltaInv);
  }

  inline void set(T ri, T rf, int n)
  {
    GridTag = LINEAR_1DGRID;
    lower_bound=ri;
    upper_bound=rf;
    num_points=n;
    // Delta is the differential spacing
    X.resize(n);
    double Delta_DP = (static_cast<double>(rf)-static_cast<double>(ri))/static_cast<double>(n-1);
    Delta = Delta_DP;
    DeltaInv = 1.0/Delta_DP;
    for(int i=0; i<n; i++)
      X[i] = ri+Delta_DP*i;
  }
  
  inline Ugrid einspline_grid()
  {
    Ugrid grid;
    grid.start = lower_bound;
    grid.end   = upper_bound;
    grid.num   = num_points;
    grid.delta = Delta;
    grid.delta_inv= DeltaInv;

    return grid;   
  }
};

/** One-Dimensional logarithmic-grid.
 *
 * The analytic form \f$ r_i = r_0
 * \left( \frac{r_f}{r_0} \right) ^{\frac{i}{N-1}}, \f$
 * where \f$ N \f$ is the number of points and the index
 * \f$ i \f$ runs from 0 to \f$ N-1 \f$
 */
template <class T, class CT=Vector<T> >
struct LogGrid: public OneDimGridBase<T,CT>
{

  using OneDimGridBase<T,CT>::Loc;
  using OneDimGridBase<T,CT>::GridTag;
  using OneDimGridBase<T,CT>::num_points;
  using OneDimGridBase<T,CT>::lower_bound;
  using OneDimGridBase<T,CT>::upper_bound;
  using OneDimGridBase<T,CT>::Delta;
  using OneDimGridBase<T,CT>::X;
  // T Delta;
  T OneOverLogDelta;

  OneDimGridBase<T,CT>* makeClone() const
  {
    return new LogGrid<T,CT>(*this);
  }

  inline void locate(T r)
  {
    Loc = static_cast<int>(std::log(r/X[0])*OneOverLogDelta);
  }

  inline void set(T ri, T rf, int n)
  {
    GridTag = LOG_1DGRID;
    lower_bound=ri;
    upper_bound=rf;
    num_points=n;
    // r(i) = ri*(rf/ri)^(i/(n-1))
    // this expression is equal to:
    // r(i) = ri*exp(dlog_ratio*i)
    // where dlog_ratio = (1/(n-1))*log(rf/ri)
    // dlog_ratio is the differential spacing
    double ratio = rf/ri;
    double log_ratio = std::log(ratio);
    double dlog_ratio = log_ratio/static_cast<double>(n-1);
    X.resize(n);
    X[0] = ri;
    for(int i=1; i < n; i++)
      X[i] = ri*std::exp(i*dlog_ratio);
    Delta = dlog_ratio;
    OneOverLogDelta = 1.0/dlog_ratio;
  }
};

/** One-Dimensional logarithmic-grid starting at the origin (Used in Siesta).
 *
 * The analytic form \f$ r_i = B
 * \left[ \exp(Ai)-1 \right] , \f$
 * where the number of points is \f$ N \f$ and the index
 * \f$ i \f$ runs from 0 to \f$ N-1 \f$
 */
template <class T, class CT=Vector<T> >
struct LogGridZero: public OneDimGridBase<T,CT>
{

  using OneDimGridBase<T,CT>::Loc;
  using OneDimGridBase<T,CT>::GridTag;
  using OneDimGridBase<T,CT>::num_points;
  using OneDimGridBase<T,CT>::lower_bound;
  using OneDimGridBase<T,CT>::upper_bound;
  using OneDimGridBase<T,CT>::Delta;
  using OneDimGridBase<T,CT>::X;
  T OneOverA;
  T OneOverB;

  OneDimGridBase<T,CT>* makeClone() const
  {
    return new LogGridZero<T,CT>(*this);
  }

  inline void locate(T r)
  {
    Loc= static_cast<int>(std::log(r*OneOverB+1.0)*OneOverA);
  }

  /** the meaing of ri/rf are different from the convetions of other classes
   * @param ri ratio
   * @param rf norm
   * @param n number of grid, [0,n)
   */
  inline void set(T ri, T rf, int n)
  {
    GridTag = LOGZERO_1DGRID;
    lower_bound=ri;
    upper_bound=rf;
    num_points=n;
    OneOverA = 1.0/ri;
    OneOverB = 1.0/rf;
    X.resize(n);
    for(int i=0; i<n; i++)
      X[i] = rf*(std::exp(ri*i)-1.0);
    Delta = 0.0;
  }
};

/** One-Dimensional numerical grid with arbitrary grid spacings.
 *
 * Non-Analytic grid, uses an array of values
 * (typically read in from a file).
 */
template <class T, class CT=Vector<T> >
struct NumericalGrid: public OneDimGridBase<T,CT>
{

  using OneDimGridBase<T,CT>::Loc;
  using OneDimGridBase<T,CT>::GridTag;
  using OneDimGridBase<T,CT>::num_points;
  using OneDimGridBase<T,CT>::lower_bound;
  using OneDimGridBase<T,CT>::upper_bound;
  using OneDimGridBase<T,CT>::Delta;
  using OneDimGridBase<T,CT>::X;

  NumericalGrid()
  {
    GridTag = CUSTOM_1DGRID;
  }

  template<class VA>
  NumericalGrid(const VA& nv)
  {
    GridTag = CUSTOM_1DGRID;
    assign(nv.begin(),nv.end());
  }

  OneDimGridBase<T,CT>* makeClone() const
  {
    return new NumericalGrid<T,CT>(*this);
  }


  template<class IT>
  void assign(IT g_first, IT g_last)
  {
    num_points=g_last-g_first;
    X.resize(num_points);
    copy(g_first,g_last,X.begin());
    lower_bound=X[0];
    upper_bound=X[num_points-1];
    int nf=num_points/2;
    Delta=X[nf]-X[nf-1];
  }

  inline void resize(int n)
  {
    if(X.size() != n)
      X.resize(n);
  }

  inline void locate(T r)
  {
    int k;
    int klo=0;
    int khi=num_points;
    //int khi=this->size()-1;
    while(khi-klo > 1)
    {
      k=(khi+klo) >> 1;
      if(X[k] > r)
        khi=k;
      else
        klo=k;
    }
    Loc= klo;
  }

  inline void set(T ri, T rf, int n)
  {
    lower_bound=ri;
    upper_bound=rf;
    num_points=n;
    //X.resize(n);
    //Delta = 0.0;
  }
};

template<class T>
std::ostream& operator<<(std::ostream& out, const OneDimGridBase<T>& rhs)
{
  for(int i=0; i<rhs.size(); i++)
    out << i << " " << rhs.r(i) << " " << rhs.dr(i)<< std::endl;
  return out;
}

}
#endif

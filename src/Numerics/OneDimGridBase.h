//////////////////////////////////////////////////////////////////
// (c) Copyright 2003  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef OHMMS_ONEDIMGRID_BASE
#define OHMMS_ONEDIMGRID_BASE

/**@file OneDimGridBase.h
 *@brief Decalaration of One-Dimesional grids
 */

#include <algorithm>
#include "OhmmsPETE/OhmmsVector.h"

/** An abstract base class to implement a One-Dimensional grid 
 */
template <class T, class CT=Vector<T> >
struct OneDimGridBase {

  typedef T value_type;
  typedef CT Array_t;
  ///the current index of the grid
  int Loc;
  ///differential spacing of the grid
  T Delta;
  ///temporary data for interpolations
  T h,p1,p2,q1,q2,dp1,dq1,dq2,d2p1,d2q1,d2q2;

  ///array to store the radial grid data
  Array_t X;
  ///assign a value
  inline T& operator[](int i) { return X[i];}
  ///assign a value
  inline T& operator()(int i) { return X[i];}
  ///return a value
  inline T operator[](int i) const { return X[i];}
  ///return a value
  inline T operator()(int i) const { return X[i];}

  inline const T* data() const { return &(X[0]);}
  inline T* data() { return &(X[0]);}

  ///return the differential spacing of the grid
  inline T dh() const { return Delta;}
  ///returns \f$r(i)\f$
  inline T r(int i) const {return X[i];}
  ///returns \f$r(i+1)-r(i)\f$
  inline T dr(int i) const { return X[i+1]-X[i];}  
  ///returns the size of the grid
  inline int size() const { return X.size();}
  ///return the first grid point
  inline T rmin() const { return X[0];}
  ///return the last grid point
  inline T rmax() const { return X[X.size()-1];}

  ///update the variables for interpolations
  inline void update(T r) {
    int khi(Loc+1);
    h=X[khi]-X[Loc]; 
    //hinv=1.0/h;
    value_type hinv(1.0e0/h); 
    value_type t((r-X[Loc])*hinv); 
    value_type tm(t-1.0);
    p1=tm*tm*(1.0+2.0*t);
    p2=t*t*(3.0-2.0*t);
    q1=t*tm*tm;
    q2=t*t*tm;

    dp1=6.0*t*tm*hinv;
    dq1=(1.0-4.0*t+3.0*t*t);
    dq2=t*(3.0*t-2.0);

    d2p1=(12.0*t-6.0)*hinv*hinv;
    d2q1=(6.0*t-4.0)*hinv;
    d2q2=(6.0*t-2.0)*hinv;
  }

  ///assign and return the index for radial point r
  virtual int index(T r) = 0;

  inline T cubicInterpolate(T a, T b, T a1, T b1, T& du, T& d2u) {
    du = dp1*(a-b)+dq1*a1+dq2*b1;
    d2u = d2p1*(a-b)+d2q1*a1+d2q2*b1;
    return p1*a+p2*b+h*(q1*a1+q2*b1);
  }
  /**
   *@param ri initial grid point
   *@param rf final grid point
   *@param n number of grid points
   *@brief Set the grid given the parameters.
   */
  virtual void set(T ri, T rf, int n) = 0;
};

/** One-Dimensional linear-grid.
 *
 * The analytic form \f[ r_i = r_0 + 
 * i\left( \frac{r_f - r_0}{N-1} \right), \f]
 * where \f$ N \f$ is the number of points and the index
 * \f$ i \f$ runs from 0 to \f$ N-1 \f$
 */
template <class T, class CT=Vector<T> >
struct LinearGrid: public OneDimGridBase<T,CT> {

  using OneDimGridBase<T,CT>::X;
  using OneDimGridBase<T,CT>::Loc;
  using OneDimGridBase<T,CT>::Delta;

  // T Delta;
  T DeltaInv;

  inline int index(T r) {
    return Loc = static_cast<int>((r-X[0])*DeltaInv);
  }

  inline void set(T ri, T rf, int n) {
    // Delta is the differential spacing
    X.resize(n);
    Delta = (rf-ri)/static_cast<T>(n-1);
    DeltaInv = 1.0/Delta;
    X[0] = ri;
    for(int i=0; i<n-1; i++) X[i+1] = X[i]+Delta;
  }

};

/** One-Dimensional logarithmic-grid.
 *
 * The analytic form \f[ r_i = r_0 
 * \left( \frac{r_f}{r_0} \right) ^{\frac{i}{N-1}}, \f]
 * where \f$ N \f$ is the number of points and the index
 * \f$ i \f$ runs from 0 to \f$ N-1 \f$
 */
template <class T, class CT=Vector<T> >
struct LogGrid: public OneDimGridBase<T,CT> {

  using OneDimGridBase<T,CT>::X;
  using OneDimGridBase<T,CT>::Loc;
  using OneDimGridBase<T,CT>::Delta;
  // T Delta;
  T OneOverLogDelta; 
  
  inline int index(T r){
    return Loc = static_cast<int>(log(r/X[0])*OneOverLogDelta);
  }

  inline void set(T ri, T rf, int n) {
    // r(i) = ri*(rf/ri)^(i/(n-1))
    // this expression is equal to:
    // r(i) = ri*exp(dlog_ratio*i)
    // where dlog_ratio = (1/(n-1))*log(rf/ri) 
    // dlog_ratio is the differential spacing
    T ratio = rf/ri;
    T log_ratio = log(ratio);
    T dlog_ratio = log_ratio/static_cast<T>(n-1);
    T expdr = exp(dlog_ratio);
    X.resize(n);

    X[0] = ri;
    for(int i=0; i < n-1; i++) {
      X[i+1] = X[i]*expdr;
    }
    Delta = dlog_ratio;
    OneOverLogDelta = 1.0/Delta;
  }
};

/**One-Dimensional logarithmic-grid starting at the
 *origin (Used in Siesta).
 *
 * The analytic form \f[ r_i = B 
 * \left[ \exp(Ai)-1 \right] , \f]
 * where the number of points is \f$ N \f$ and the index
 * \f$ i \f$ runs from 0 to \f$ N-1 \f$
 */
template <class T, class CT=Vector<T> >
struct LogGridZero: public OneDimGridBase<T,CT> {

  using OneDimGridBase<T,CT>::X;
  using OneDimGridBase<T,CT>::Loc;
  using OneDimGridBase<T,CT>::Delta;
  T OneOverA; 
  T OneOverB;

  inline int index(T r){
    return Loc = static_cast<int>(log(r*OneOverB+1.0)*OneOverA);
  }

  inline void set(T a, T b, int n) {
    OneOverA = 1.0/a;
    OneOverB = 1.0/b;
    X.resize(n);
    for(int i=0; i<n; i++) X[i] = b*(exp(a*i)-1.0);

    Delta = 0.0;
  }
};

/** One-Dimensional numerical grid with arbitrary grid spacings. 
 *
 * Non-Analytic grid, uses an array of values
 * (typically read in from a file).
 */
template <class T, class CT=Vector<T> >
struct NumericalGrid: public OneDimGridBase<T,CT> {

  using OneDimGridBase<T,CT>::X;
  using OneDimGridBase<T,CT>::Loc;
  using OneDimGridBase<T,CT>::Delta;
 
  template<class VA>
  NumericalGrid(const VA& nv) {
    // NumericalGrid(const std::vector<T>& nv) {
    X.resize(nv.size());
    std::copy(nv.begin(), nv.end(), X.data());
  }

  inline int index(T r){
    int k;
    int klo=0;
    int khi=this->size()-1;
    while(khi-klo > 1){
      k=(khi+klo) >> 1;
      if(X[k] > r) khi=k;
      else klo=k;
    }
    return Loc = klo;
  }
  
  inline void set(T ri, T rf, int n) {
    X.resize(n);
    Delta = 0.0;
  }
};



template<class T>
ostream& operator<<(ostream& out, const OneDimGridBase<T>& rhs)
{
  for(int i=0; i<rhs.size(); i++)
    out << i << " " << rhs.r(i) << " " << rhs.dr(i)<< endl;
  return out;
}
#endif

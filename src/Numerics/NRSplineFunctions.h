//////////////////////////////////////////////////////////////////
// (c) Copyright 2003  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef QMCPLUSPLUS_NR_CUBICSPLINE_H
#define QMCPLUSPLUS_NR_CUBICSPLINE_H

/**template function: converted from Numerical Recipe spline.c
 *note that the range of data is [0,n) instead of [1,n]
 */
template<typename Tg, typename T>
inline void 
NRCubicSpline(const Tg* x, const T* y, int n, T yp1, T ypn, T* y2) {

  int i,k;
  T p,qn,sig,un;
  std::vector<T> u(n);
  if (yp1 > 0.99e30)
    y2[0]=u[0]=0.0;
  else {
    y2[0] = -0.5;
    u[0]=(3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
  }
  for (i=1;i<n-1;i++) {
    sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
    p=sig*y2[i-1]+2.0;
    y2[i]=(sig-1.0)/p;
    u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
    u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
  }
  if (ypn > 0.99e30)
    qn=un=0.0;
  else {
    qn=0.5;
    un=(3.0/(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
  }
  y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.0);
  for (k=n-2;k>=0;k--)
    y2[k]=y2[k]*y2[k+1]+u[k];
}

/**template function: converted from Numerical Recipe spline.c and with PBC
 * @param x first address of the x-axis
 * @param y first address of the y-axis
 * @param n number of data point
 * @param d1 first address of the first-derivatives coefficients
 * @param d2 first address of the second-derivatives coefficients
 *
 *Note that the range of data is [0,n) instead of [1,n]
 */
template<typename Tg, typename T>
inline void 
NRCubicSplinePBC(const Tg* x, const T* y, int n, T* d1, T* d2 ) {

  //store the last valid index
  int last(n-1),last1(n-2),last2(n-3);

  //temporary array for the d1(frist) and d3(third) derivatives
  std::vector<T> d3(n);

  //clear the contents
  d1[0]=d2[0]=d3[0]=0.0;
  d1[last]=d2[last]=d3[last]=0.0;

  if(n == 2) {
    d1[0]=(y[1]-y[0])/(x[1]-x[0]);
    d1[1]=d1[0];
  } else {
    d3[0]=x[1]-x[0];
    d2[1]=(y[1]-y[0])/d3[0];
    for(int cur=1, next=2; cur<last; cur++, next++) {
      d3[cur]=x[next]-x[cur];
      d1[cur]=2.0*(d3[cur-1]+d3[cur]);
      d2[next]=(y[next]-y[cur])/d3[cur];
      d2[cur]=d2[next]-d2[cur];
    }
  }

  //Time to take care of BC
  std::vector<T> wk(n,0.0);
  d1[0]=2.0*(d3[0]+d3[last1]);
  d2[0]=(y[1]-y[0])/d3[0]-(y[last]-y[last1])/d3[last1];
  wk[0]=d3[last1];
  wk[last2]=d3[last2];
  wk[last1]=d3[last1];

  //forward elimination
  for(int cur=1, prev=0; cur<last1; cur++,prev++) {
    T t(d3[prev]/d1[prev]);
    d1[cur] -= t*d3[prev];
    d2[cur] -= t*d2[prev];
    wk[cur] -= t*wk[prev];

    T q(wk[last1]/d1[prev]);
    wk[last1]=-q*d3[prev];
    d1[last1]-= q*wk[prev];
    d2[last1]-= q*d2[prev];
  }
  //correct the last-1
  wk[last1]+=d3[last2];

  T t(wk[last1]/d1[last2]);
  d1[last1] -= t*wk[last2];
  d2[last1] -= t*d2[last2];

  //back substitution
  d2[last1] /= d1[last1];
  d2[last2] = (d2[last2]-wk[last2]*d2[last1])/d1[last2];
  for(int cur=last2-1,next=last2; next>0; cur--,next--) {
    d2[cur]=(d2[cur]-d3[cur]*d2[next]-wk[cur]*d2[last1])/d1[cur];
  }
  d2[last]=d2[0];

  for(int cur=0,next=1; cur<last; cur++,next++) {
    d1[cur]=(y[next]-y[cur])/d3[cur]-d3[cur]*(d2[next]+2.0*d2[cur]);
    d3[cur]=6.0*(d2[next]-d2[cur])/d3[cur];
    d2[cur]*=6.0;
  }
  d1[last]=d1[0];
  d2[last]=d2[0];
}

template<typename Tg, typename T>
inline void 
NRCubicSplineFirst(const Tg* x, const T* y, int n, T* d1, T* d2 ) {

  int last(n-1),last1(n-2),last2(n-3);

  //store the last valid index
  T a1=d1[0];
  T an=d1[last];

  //temporary array for the d1(frist) and d3(third) derivatives
  std::vector<T> d3(n);

  //clear the contents
  d2[0]=d3[0]=0.0;
  d1[last]=d2[last]=d3[last]=0.0;

  if(n == 2) {
    d1[0]=(y[1]-y[0])/(x[1]-x[0]);
    d1[1]=d1[0];
  } else {
    d3[0]=x[1]-x[0];
    d2[1]=(y[1]-y[0])/d3[0];
    for(int cur=1, next=2; cur<last; cur++, next++) {
      d3[cur]=x[next]-x[cur];
      d1[cur]=2.0*(d3[cur-1]+d3[cur]);
      d2[next]=(y[next]-y[cur])/d3[cur];
      d2[cur]=d2[next]-d2[cur];
    }
  }

  T elem21(d3[0]), elemnn1(d3[last1]);

  d1[0]=2.0*d3[0];
  d2[0]=(y[1]-y[0])/d3[0]-a1;

  d1[last]=2.0*d3[last1];
  d2[last]=-(y[last]-y[last1])/d3[last1]+an;

  //forward elimination
  T t1(elem21/d1[0]);
  d1[1] -= t1*d3[0];
  d2[1] -= t1*d2[0];
  for(int cur=2, prev=1; cur<last; cur++, prev++) {
    T t(d3[prev]/d1[prev]);
    d1[cur] -= t*d3[prev];
    d2[cur] -= t*d2[prev];
  }
  t1=elemnn1/d1[last1];
  d1[last] -= t1*d3[last1];
  d2[last] -= t1*d2[last1];

  //back substitution
  d2[last] /= d1[last];

  for(int cur=last1; cur>=0; cur--) {
    d2[cur]=(d2[cur]-d3[cur]*d2[cur+1])/d1[cur];
  }
  //d2[0]=(d2[0]-d3[0]*d2[1])/d1[0];

  d3[0]=x[1]-x[0];
  d3[last1]=x[last]-x[last1];
  for(int cur=0,next=1; cur<last; cur++,next++) {
    d1[cur]=(y[next]-y[cur])/d3[cur]-d3[cur]*(d2[next]+2.0*d2[cur]);
    d3[cur] =6.0*(d2[next]-d2[cur])/d3[cur];
    d2[cur]*=6.0;
  }
  T hn(x[last]-x[last1]);
  d1[last]=an;
  d2[last]=d2[last1]+hn*d3[last1];
}

/** Solve Tridiagonal problem
 */
template<class CT>
struct TriDiagSolver
{
  static void solve(const CT& a, const CT& b, const CT& c, const CT& d, CT& p, int n)
  {
    typedef typename CT::value_type value_type;

    value_type bet;
    CT gamma(n);

    bet=b[0];
    p[0]=d[0]/bet;

    for(int j=1; j<n; j++)
    {
      gamma[j]=c[j-1]/bet;
      bet=b[j]-a[j]*gamma[j];
      p[j]=(d[j]-a[j]*p[j-1])/bet;
    }
    for(int j=n-2; j>=0; j--)
      p[j]=p[j]-gamma[j+1]*p[j+1];
  }

  //static void check(const CT& a, const CT& b, const CT& c, const CT& d, CT& p, int n)
  //{
  //  typedef typename CT::value_type value_type;

  //  qmcplusplus::Matrix<value_type> M(n,n);
  //  for(int i=0; i<n; i++) M(i,i)=b[i];
  //  for(int i=0; i<n-1; i++) M(i+1,i)=a[i];
  //  for(int i=0; i<n-1; i++) M(i,i+1)=c[i];

  //  CT res(n,0.0);
  //  for(int i=0; i<n; i++)
  //  {
  //    value_type s=0;
  //    for(int j=0; j<n; j++) s += M(i,j)*p[j];
  //    res[i]=s;
  //  }

  //  for(int i=0; i<n; i++)
  //    std::cout << d[i] << " " << res[i] << std::endl;
  //}
};
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#ifndef QMCPLUSPLUS_NRC_OPTIMIZATION_H
#define QMCPLUSPLUS_NRC_OPTIMIZATION_H

#include "Numerics/MatrixOperators.h"
#include "Numerics/DeterminantOperators.h"
#include "Numerics/LinearFit.h"
#include "Configuration.h"
#include <math.h>
#if (__GNUC__ == 2)
#include <algo.h>
#else
#include <algorithm>
#endif
#include <limits>
//template<class T>
//inline void shift (T& a, T& b, T& c, T d) {
//  a = b; b= c; c = d;
//}
//double sign(double a, double b) {
//  return (b > 0.0)? std::abs(a): -std::abs(a);
//}

template<class T>
struct sign2 { };

template<>
struct sign2<double>
{
  inline static double apply(double a, double b)
  {
    return (b > 0.0)? std::abs(a): -std::abs(a);
  }
};

template<>
struct sign2<float>
{
  inline static float apply(float a, float b)
  {
    return (b > 0.0)? std::abs(a): -std::abs(a);
  }
};

template<class T>
struct NRCOptimization
{

  typedef T Return_t;

  /** number of line iteration for Brent method */
  int ITMAX;

  /** maximum CG mixture for line minimization
   *
   * y'=y + Lambda*cg where Lambda < LambdaMax
   */
  Return_t LambdaMax;
  Return_t quadstep;
  Return_t quadoffset;
  Return_t largeQuarticStep;
  bool validFuncVal;
  //if tol is absolute, not a percent
  bool AbsFuncTol;

  int current_step;

  NRCOptimization()
  {
    ITMAX=100;
    ZEPS = 1.0e-10;
    CGOLD = 0.3819660e0;
    GOLD = 1.618034e0;
    TOL = 2.0e-4;
    GLIMIT = 100.0;
    TINY = std::numeric_limits<T>::epsilon();
    LambdaMax = 0.02;
    current_step = 0;
    quadstep=0.01;
    quadoffset=0.0;
    largeQuarticStep=2.0;
    validFuncVal=true;
    AbsFuncTol=false;
  }

  virtual ~NRCOptimization() { }

  /** evaluate the value for y+dl*x
   *
   * Lineminimization uses this function to find the minimum along the x direction
   */
  virtual Return_t Func(Return_t dl) = 0;

  Return_t Lambda;
  Return_t ZEPS, CGOLD, TOL, GLIMIT, TINY, GOLD;

  // Returns the number of real roots
  inline int CubicFormula (double a, double b, double c, double d,
                           double &x1, double &x2, double &x3)
  {
    double A = b/a;
    double B = c/a;
    double C = d/a;
    double Q = (A*A - 3.0*B)/9.0;
    double R = (2.0*A*A*A - 9.0*A*B + 27.0*C)/54.0;
    //cerr << "Q = " << Q << " R = " << R << "\n";
    if ((R*R) < (Q*Q*Q))
    {
      double theta = acos(R/sqrt(Q*Q*Q));
      double twosqrtQ = 2.0*sqrt(Q);
      double third = 1.0/3.0;
      double thirdA = third * A;
      x1 = -twosqrtQ*cos(third*theta) - thirdA;
      x2 = -twosqrtQ*cos(third*(theta + 2.0*M_PI)) - thirdA;
      x3 = -twosqrtQ*cos(third*(theta - 2.0*M_PI)) - thirdA;
      return 3;
    }
    else
    {
      double D = -Q*Q*Q + R*R;
      double u = cbrt(-R + sqrt(D));
      double v = cbrt(-R - sqrt(D));
      double y1 = u+v;
      x1 = y1 - A/3.0;
      return 1;
    }
  }

  inline Return_t QuarticMinimum (std::vector<Return_t> &coefs)
  {
    double a, b, c, d;
    a = 4.0*coefs[4];
    b = 3.0*coefs[3];
    c = 2.0*coefs[2];
    d = coefs[1];
    double x1, x2, x3;
    int numroots = CubicFormula (a, b, c, d, x1, x2, x3);
    if (numroots == 1)
      return x1;
    else
    {
      double v1 = coefs[0] + coefs[1]*x1 + coefs[2]*x1*x1 + coefs[3]*x1*x1*x1
                  + coefs[4]*x1*x1*x1*x1;
      double v2 = coefs[0] + coefs[1]*x2 + coefs[2]*x2*x2 + coefs[3]*x2*x2*x2
                  + coefs[4]*x2*x2*x2*x2;
      double v3 = coefs[0] + coefs[1]*x3 + coefs[2]*x3*x3 + coefs[3]*x3*x3*x3
                  + coefs[4]*x3*x3*x3*x3;
      if (v1 < v2 && v1 < v3)
        return x1;
      if (v2 < v1 && v2 < v3)
        return x2;
      if (v3 < v1 && v3 < v2)
        return x3;
      return x1;
    }
  }



  bool lineoptimization()
  {
    std::vector<Return_t> x(5), y(5), coefs(5), deriv(4);
    qmcplusplus::Matrix<Return_t> S(5,5);
    x[0]=-2*quadstep;
    x[1]=-quadstep;
    x[2]=0.0;
    x[3]=quadstep;
    x[4]=2*quadstep;
    Return_t start_cost, cost;
    validFuncVal=true;
    for (int i=0; i<5; i++)
    {
      y[i] = Func(x[i]);
      for (int j=0; j<5; j++)
        S(i,j) = std::pow(x[i],j);
    }
    start_cost = y[2];
    if(validFuncVal)
    {
      qmcplusplus::invert_matrix(S, false);
      qmcplusplus::MatrixOperators::product(S, &(y[0]), &(coefs[0]));
      Lambda = QuarticMinimum (coefs);
#if (__cplusplus >= 201103L)
      if (std::abs(Lambda) > largeQuarticStep || std::isnan(Lambda))
        return lineoptimization2();
      cost = Func(Lambda);
      if (std::isnan(cost) || cost > start_cost)
        return lineoptimization2();
#else
      if (std::abs(Lambda) > largeQuarticStep || isnan(Lambda))
        return lineoptimization2();
      cost = Func(Lambda);
      if (isnan(cost) || cost > start_cost)
        return lineoptimization2();
#endif
    }
    else
    {
      return lineoptimization2();
    }
    //fprintf (stderr, "Minimum found at %1.8f\n", Lambda);
    current_step++;
    return true;
    // HACK HACK HACK
//     if (Lambda < 0.0) {
    // char fname[50];
    // snprintf (fname, 50, "line_opt_%d.dat", current_step);
    // FILE *fout = fopen (fname, "w");
    // for (double lam=-0.01; lam<=0.01; lam+=0.0001) {
    //   double val = 0.0;
    //   for (int j=0; j<5; j++)
    // 	val += coefs[j] * std::pow(lam, j);
    //   fprintf (fout, "%1.8f %1.12e %1.12e\n", lam, Func(lam), val);
    // }
    // fclose(fout);
//     }
    // END HACK HACK HACK
  }

  bool lineoptimization3(int points, Return_t& zeroCost)
  {
    // quartic fit with variable number of points for input.
    //  least squares solver
    std::vector<Return_t> x(points), y(points), coefs(5);
    qmcplusplus::Matrix<Return_t> S(points,5);
    for(int i=0; i<points; i++)
      x[i]=Return_t(i-1)*quadstep;
    Return_t start_cost, cost;
    std::vector<bool> cFailed(points,false);
    int nFailed(0);
    validFuncVal=true;
    y[0] = Func(x[0]);
    if (!validFuncVal)
    {
      cFailed[0]=true;
      nFailed++;
    }
    y[1]=zeroCost;
    y[2] = Func(x[2]);
    if (!validFuncVal)
    {
      cFailed[2]=true;
      nFailed++;
    }
    for (int j=0; j<5; j++)
    {
      S(2,j) = std::pow(x[2],j);
      S(1,j) = 0.0;
      S(0,j) = std::pow(x[0],j);
    }
    //if our guess is in the wrong direction, flip the line search
    if (y[0]<y[1])
      for (int i=3; i<points; i++)
        x[i]*=-1.0;
    else
      if ((y[0]>y[1])&&(y[2]>y[1]))
      {
        Return_t stp=x[2]-x[0];
        for (int i=1; i<points-2; i++)
          x[i]=x[0]+1.0*i*stp;
      }
    for (int i=3; i<points; i++)
    {
      y[i] = Func(x[i]);
      if (!validFuncVal)
      {
        cFailed[i]=true;
        nFailed++;
      }
      for (int j=0; j<5; j++)
        S(i,j) = std::pow(x[i],j);
    }
    start_cost = y[1];
//     for (int i=0; i<points; i++) std::cout <<x[i]<<": "<<y[i]<< std::endl;
    if (nFailed>0)
    {
      int ok_pts(points-nFailed);
      if (ok_pts>=5)
      {
        //some failed but we still have enough to do a fit
        std::vector<Return_t> xp(ok_pts), yp(ok_pts);
        qmcplusplus::Matrix<Return_t> Sp(ok_pts,5);
        for (int i=0,ip=0; i<points; i++)
        {
          if (!cFailed[i])
          {
            yp[ip]=y[i];
            xp[ip]=x[i];
            for (int j=0; j<5; j++)
              Sp(ip,j) = S(i,j);
            ip++;
          }
        }
        x=xp;
        y=yp;
        S=Sp;
        validFuncVal=true;
      }
      else
        validFuncVal=false;
    }
    if(validFuncVal)
    {
      qmcplusplus::LinearFit(y,S,coefs);
      Lambda = QuarticMinimum (coefs);
#if (__cplusplus >= 201103L)
      if (std::abs(Lambda) > largeQuarticStep || std::isnan(Lambda) || (Lambda==0.0))
        return lineoptimization2(largeQuarticStep);
      zeroCost = Func(Lambda);
//       std::cout <<"Start Cost:"<< start_cost<<" Lambda:"<<Lambda<<" FinalCost:"<<cost<< std::endl;
      if (std::isnan(zeroCost) || zeroCost > start_cost)
        return lineoptimization2(largeQuarticStep);
#else
      if (std::abs(Lambda) > largeQuarticStep || isnan(Lambda) || (Lambda==0.0))
        return lineoptimization2(largeQuarticStep);
      zeroCost = Func(Lambda);
//       std::cout <<"Start Cost:"<< start_cost<<" Lambda:"<<Lambda<<" FinalCost:"<<cost<< std::endl;
      if (isnan(zeroCost) || zeroCost > start_cost)
        return lineoptimization2(largeQuarticStep);
#endif
    }
    else
    {
      return false;
    }
    //fprintf (stderr, "Minimum found at %1.8f\n", Lambda);
    current_step++;
    return true;
    // HACK HACK HACK
//     if (Lambda < 0.0) {
    // char fname[50];
    // snprintf (fname, 50, "line_opt_%d.dat", current_step);
    // FILE *fout = fopen (fname, "w");
    // for (double lam=-0.01; lam<=0.01; lam+=0.0001) {
    //   double val = 0.0;
    //   for (int j=0; j<5; j++)
    //    val += coefs[j] * std::pow(lam, j);
    //   fprintf (fout, "%1.8f %1.12e %1.12e\n", lam, Func(lam), val);
    // }
    // fclose(fout);
//     }
    // END HACK HACK HACK
  }

  bool lineoptimization2(Return_t maxStep=1e9)
  {
    Return_t ax = 0;
    Return_t bx(0), fa, fx, fb;
    Return_t xx = LambdaMax;
    // HACK HACK HACK
//     char fname[50];
//     snprintf (fname, 50, "line_opt_%d.dat", current_step++);
//     FILE *fout = fopen (fname, "w");
//     for (double lam=-0.01; lam<=0.01; lam+=0.0001)
//       fprintf (fout, "%1.8f %1.12e\n", lam, Func(lam));
//     fclose(fout);
    // END HACK HACK HACK
    bool success=true;
    validFuncVal=true;
    qmcplusplus::app_log()<<"Before:  ax = "<<ax<<"  bx="<<xx<<"  cx="<<bx<< std::endl;
    success=mnbrakNRC(ax,xx,bx,fa,fx,fb,maxStep);
    if((!success && !validFuncVal) || (success && !validFuncVal))
    {
      Lambda = 0.0;
      qmcplusplus::app_log()<<"Problems bracketing minimum.\n";
      return false;
    }
    else if(!success && validFuncVal)
      {
// in case the function is unable to bracket the minimum but
// still finds a point with lower cost and validFuncVal, take this point
        Lambda=xx;
        qmcplusplus::app_log()<<"Problems bracketing minimum. Lower Value returned.\n";
        return false;
      }
    qmcplusplus::app_log()<<"After:  ax = "<<ax<<"  bx="<<xx<<"  cx="<<bx<< std::endl;
    Lambda = 0.0e0;
    Return_t ep = brentNRC(ax,xx,bx,Lambda);
    if(validFuncVal)
    {
      qmcplusplus::app_log()<<"Minimum found at lambda = "<<Lambda<< std::endl;
      if(std::abs(Lambda)<TINY)
        return false;
      else
        return true;
    }
    else
    {
      Lambda = 0.0;
      qmcplusplus::app_log()<<"Problems bracketing minimum.\n";
      return false;
    }
  }

  inline void shift(Return_t& a, Return_t& b, Return_t& c, Return_t d)
  {
    a = b;
    b= c;
    c = d;
  }

  T brentNRC(Return_t ax, Return_t bx, Return_t cx, Return_t & xmin);

  bool mnbrakNRC(Return_t& ax,Return_t& bx,Return_t& cx,
                 Return_t& fa,Return_t& fb,Return_t& fc, Return_t maxStep=1e9);

};

template<class T>
T NRCOptimization<T>::brentNRC(Return_t ax, Return_t bx,  Return_t cx, Return_t& xmin)
{
  Return_t a,b,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
  Return_t e=0.0, d=0.0;
  validFuncVal=true;
  a=((ax < cx) ? ax : cx);
  b=((ax > cx) ? ax : cx);
  x=w=v=bx;
  fw=fv=fx=Func(x);
  if(!validFuncVal)
    return 0.0;
  for(int iter=1; iter<=ITMAX; iter++)
  {
    xm=0.5*(a+b);
    if (AbsFuncTol)
    {
      tol2=2.0*(tol1=TOL);
      if (std::abs(x-xm) <= tol2)
      {
        xmin=x;
        return fx;
      }
    }
    else
    {
      tol2=2.0*(tol1=TOL*std::abs(x)+ZEPS);
      if (std::abs(x-xm) <= (tol2-0.5*(b-a)))
      {
        xmin=x;
        return fx;
      }
    }
    if (std::abs(e) > tol1)
    {
      r=(x-w)*(fx-fv);
      q=(x-v)*(fx-fw);
      p=(x-v)*q-(x-w)*r;
      q=2.0*(q-r);
      if (q > 0.0)
        p = -p;
      q=std::abs(q);
      etemp=e;
      e=d;
      if (std::abs(p) >= std::abs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
        d=CGOLD*(e=(x >= xm ? a-x : b-x));
      else
      {
        d=p/q;
        u=x+d;
        //if (u-a < tol2 || b-u < tol2) d=sign(tol1,xm-x);
        if (u-a < tol2 || b-u < tol2)
          d=sign2<T>::apply(tol1,xm-x);
      }
    }
    else
    {
      d=CGOLD*(e=(x >= xm ? a-x : b-x));
    }
    //u=(std::abs(d) >= tol1 ? x+d : x+sign(tol1,d));
    u=(std::abs(d) >= tol1 ? x+d : x+sign2<T>::apply(tol1,d));
    fu = Func(u); // fu=(*f)(u);
    if(!validFuncVal)
      return 0.0;
    if (fu <= fx)
    {
      if (u >= x)
        a=x;
      else
        b=x;
      shift(v,w,x,u);
      shift(fv,fw,fx,fu);
    }
    else
    {
      if (u < x)
        a=u;
      else
        b=u;
      if (fu <= fw || w == x)
      {
        v=w;
        w=u;
        fv=fw;
        fw=fu;
      }
      else
        if (fu <= fv || v == x || v == w)
        {
          v=u;
          fv=fu;
        }
    }
  }
  xmin=x;
  return fx;
}

template<class T>
bool NRCOptimization<T>::mnbrakNRC(Return_t& ax, Return_t& bx, Return_t& cx,
                              Return_t& fa, Return_t& fb, Return_t& fc, Return_t maxStep)
{
  Return_t ulim,u,r,q,fu,dum = 0.0e0;
  validFuncVal=true;
  fa = Func(ax); // *fa=(*func)(*ax);
  if(!validFuncVal)
    return false;
  fb = Func(bx); // *fb=(*func)(*bx);
  if(!validFuncVal)
  {
    validFuncVal=true;
    bx = ax;
    return false;
  }
  if (fb > fa)
  {
    shift(dum,ax,bx,dum);
    shift(dum,fb,fa,dum);
  }
  cx=bx+GOLD*(bx-ax);
  fc = Func(cx); // *fc=(*func)(*cx);
  if(!validFuncVal)
  {
    validFuncVal=true;
    return false;
  }
  while (fb > fc)
  {
    r=(bx-ax)*(fb-fc);
    q=(bx-cx)*(fb-fa);
    //u=(bx)-((bx-cx)*q-(bx-ax)*r)/(2.0*sign(std::max<T>(std::abs(q-r),TINY),q-r));
    u=(bx)-((bx-cx)*q-(bx-ax)*r)/(2.0*sign2<T>::apply(std::max<T>(std::abs(q-r),TINY),q-r));
    ulim=(bx)+GLIMIT*(cx-bx);
    if ((bx-u)*(u-cx) > 0.0)
    {
      fu = Func(u); // fu=(*func)(u);
// this is a problematic case, since both bx,cx is good,
// but u, which is in between {bx,cx} is bad.
      if(!validFuncVal)
      {
        validFuncVal=true;
        return false;
      }
      if (fu < fc)
      {
        ax=bx;
        bx=u;
        fa=fb;
        fb=fu;
        return true;
      }
      else if (fu > fb)
        {
          cx=u;
          fc=fu;
          return true;
        }
      u=cx+GOLD*(cx-bx);
      fu = Func(u);//fu=(*func)(u);
      if(!validFuncVal)
      {
        bx=cx;
        validFuncVal=true;
        return false;
      }
    }
    else if ((cx-u)*(u-ulim) > 0.0)
      {
        fu = Func(u);//fu=(*func)(u);
        if(!validFuncVal)
        {
          bx=cx;
          validFuncVal=true;
          return false;
        }
        if (fu < fc)
        {
          shift(bx,cx,u,cx + GOLD*(cx-bx));
          shift(fb,fc,fu,Func(u)) ;
          if(!validFuncVal)
          {
            bx=cx;
            validFuncVal=true;
            return false;
          }
        }
      }
      else if ((u-ulim)*(ulim-cx) >= 0.0)
        {
          u=ulim;
          fu = Func(u);//fu=(*func)(u);
          if(!validFuncVal)
          {
            bx=cx;
            validFuncVal=true;
            return false;
          }
        }
        else
        {
          u=cx+GOLD*(cx-bx);
          fu = Func(u);//fu=(*func)(u);
          if(!validFuncVal)
          {
            bx=cx;
            validFuncVal=true;
            return false;
          }
        }
    shift(ax,bx,cx,u);
    shift(fa,fb,fc,fu);
//     if we are out of bounds totally then return false
    if((std::abs(ax)>maxStep) and (std::abs(bx)>maxStep) and (std::abs(cx)>maxStep))
      return false;
  }
  return true;
}

#endif

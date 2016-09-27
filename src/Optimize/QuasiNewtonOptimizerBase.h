//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#ifndef QN_H
#define QN_H
#include <iostream>
#include <fstream>
#include <vector>
//#define debugmsg(msg) {std::cout << msg << std::endl;}
#define debugmsg(msg) { }


template<class T>
class QuasiNewtonOptimizerBase
{

protected:

  T DeltaFn;
  T HH; //hh
  T Epsilon; //eps
  int OptMode; //mode
  int MaxFn; //maxfn
  int PrintMode; //iprint

  virtual T Funct(int nparam, T* data) = 0;

public:

  QuasiNewtonOptimizerBase()
  {
    DeltaFn = -0.5;
    HH = 1e-6;
    Epsilon = 1e-6;
    OptMode = 1;
    MaxFn = 1;//MaxFn = 1000;
    PrintMode = 1;
  }

  void search (int nparam, T* val, std::ofstream& OutStream)
  {
    typedef std::vector<T> Array_t;
    Array_t h;
    Array_t w;
    Array_t xm;
    Array_t g;
    Array_t x;
    T aeps, alpha, tot, ff;
    int n;
    int link, itn;
    int ifn;
    T df;
    T dmin;
    T f, f1, f2;
    T z, zz;
    int intt;
    int np, n1, nn, is, iu;
    int iv, ib, idiff, iexit;
    int ij, i1;
    T sig, gys, gs0, dgs;
    n = nparam;
    OutStream << "Entry into va10a" << std::endl;
    np = n+1;
    n1 = n-1;
    nn = (n*np)/2;
    is = n;
    iu = n;
    iv = n + n;
    ib = iv + n;
    idiff = 1;
    iexit = 0;
    h.resize(n*np/2+1,0.0);///n*np/2
    w.resize(4*n+1,0.0);
    xm.resize(n+1,0.0);
    g.resize(n+1,0.0);
    x.resize(n+1,0.0);
    for(int i=1; i<=n*np/2+1; i++)
      h[i] = 0.0;
    for(int i=1; i<=4*n+1; i++)
      w[i] = 0.0;
    for(int i=1; i<=n+1; i++)
    {
      xm[i] = 0.0;
      g[i] = 0.0;
      x[i] = 0.0;
    }
    OutStream << "x = ";
    for(int i=1; i<=n; i++)
    {
      x[i] = val[i-1];
      xm[i] = std::abs(x[i]) + 0.01;
      OutStream << x[i] << " ";
    }
    OutStream << std::endl;
    f = Funct(n,&x[1]);
    OutStream << "f = " << f << std::endl;
    ///Section A
    if(OptMode  != 3)
    {
      if(OptMode == 2)
      {
        ij = 1;
        for(int i=2; i <= n; i++)
        {
          z = h[ij];
          if (z <= 0)
            return;
          ij++;
          i1 = ij;
          for(int j=i; j <= n; j++)
          {
            zz = h[ij];
            h[ij] /= z;
            int jk = ij;
            int ik = i1;
            for(int k=i; k <= j; k++)
            {
              jk += np -k;
              h[jk] -= h[ik] * zz;
              ik++;
            }
          }
        }
        if(h[ij] <= 0)
          return;
      }
      else
        /// OptMode == 1
      {
        ij = nn + 1;
        for(int i=1; i <= n; i++)
        {
          for(int j=1; j <= i; j++)
          {
            ij--;
            h[ij] = 0.0;
          }
          h[ij] = 1;
        }
      }
    }
    ij = np;
    dmin = h[1];
    int i=2;
    while(i<= n)
    {
      if(h[ij] < dmin)
        dmin = h[ij];
      ij = ij+np-i;
      i++;
    }
    ///Section C
    if(dmin <= 0)
    {
      OutStream << "dmin <= 0" << std::endl;
      return;
    }
    z = f;
    itn = 0;
    f = Funct(n,&x[1]);
    OutStream << "f = " << f << std::endl;
    ifn = 1;
    df = DeltaFn;
    if(DeltaFn == 0)
      df = f - z;
    if(DeltaFn < 0)
      df = std::abs(df * f);
    if(df <= 0)
      df = 1;
funct17:
    ;
    debugmsg("inside 17");
    for(int i=1; i <= n; i++)
      w[i] = x[i];
    link = 1;
    if((idiff - 1) < 0)
      goto funct100;
    else
      if((idiff -1) == 0)
        goto funct100;
      else
        goto funct110;
funct18:
    ;
    debugmsg("inside 18");
    if(ifn >= MaxFn)
      goto funct90;
    goto funct20;
funct20:
    ;
    debugmsg("inside 20");
    if(PrintMode == 0)
      goto funct21;
    if((itn%PrintMode) != 0)
      goto funct21;
    OutStream << "itn = " << itn << " ifn = " << ifn << std::endl;
    OutStream << "f = " << f << std::endl;
    if(PrintMode < 0)
      goto funct21;
    OutStream << "x = ";
    for(int i=1; i<=n; i++)
      OutStream << x[i] << " ";
    OutStream << std::endl;
    OutStream << "g = ";
    for(int i=1; i<=n; i++)
      OutStream << g[i] << " ";
    OutStream << std::endl;
    goto funct21;
funct21:
    ;
    debugmsg( "inside 21");
    itn++;
    w[1] = -1*g[1];
    for(int i=2; i<=n; i++)
    {
      ij = i;
      i1 = i -1;
      z = -1*g[i];
      for(int j=1; j<=i1; j++)
      {
        z -= h[ij] * w[j];
        ij += n - j;
      }
      w[i] = z;
    }
    w[is+n] = w[n] / h[nn];
    ij = nn;
    for(int i=1; i<=n1; i++)
    {
      ij--;
      z = 0;
      for(int j=1; j<=i; j++)
      {
        z += h[ij] * w[is+np-j];
        ij--;
      }
      w[is+n-i] = w[n-i] / h[ij] - z;
    }
    z = 0;
    gs0 = 0;
    for(int i=1; i<=n; i++)
    {
      if(z * xm[i] < std::abs(w[is+i]))
        z = std::abs(w[is+i]) / xm[i];
      gs0 += g[i] * w[is+i];
    }
    aeps = Epsilon / z;
    iexit = 2;
    if(gs0 >= 0)
      goto funct92;
    alpha = -2.0 * df / gs0;
    if(alpha > 1)
      alpha = 1;
    ff = f;
    tot = 0;
    intt = 0;
    iexit = 1;
    goto funct30;
    /////////////////////////
funct30:
    ;
    debugmsg("inside 30");
    if(ifn >= MaxFn)
      goto funct90;
    for(int i=1; i<=n; i++)
    {
      w[i] = x[i] + alpha * w[is+1];
    }
    f1 = Funct(n,&w[1]);
    debugmsg("f1 = " << f1);
    ifn++;
    if(f1 >= f)
      goto funct40;
    f2 = f;
    tot += alpha;
    goto funct32;
funct32:
    ;
    debugmsg("inside 32");
    for(int i=1; i<=n; i++)
      x[i] = w[i];
    f = f1;
    if((intt -1) < 0)
      goto funct35;
    else
      if((intt-1)==0)
        goto funct49;
      else
        goto funct50;
funct35:
    ;
    debugmsg("inside 35");
    if(ifn >= MaxFn)
      goto funct90;
    for(int i=1; i<=n; i++)
    {
      w[i] = x[i] + alpha * w[is+i];
    }
    f1 = Funct(n,&w[1]);
    debugmsg("f1 = " << f1);
    debugmsg("f2 = " << f2);
    ifn++;
    if(f1 >= f)
      goto funct50;
    if(((f1 + f2) >= (f + f)) && (7*f1+5*f2 > 12*f))
      intt = 2;
    tot += alpha;
    alpha *= 2;
    goto funct32;
funct40:
    ;
    debugmsg("inside 40");
    if(alpha < aeps)
      goto funct92;
    if (ifn >= MaxFn)
      goto funct90;
    alpha *= 0.5;
    for(int i=1; i <= n; i++)
    {
      w[i] = x[i] + alpha * w[is+i];
    }
    f2 = Funct(n,&w[1]);
    ifn++;
    if(f2 >= f)
      goto funct45;
    tot += alpha;
    f = f2;
    for(int i=1; i<=n; i++)
      x[i] = w[i];
    goto funct49;
funct45:
    ;
    debugmsg("inside 45");
    z = 0.1;
    if((f1 + f) > (f2 + f2))
      z = 1.0 + 0.5*(f - f1)/(f + f1 - f2 - f2);
    if(z < 0.1)
      z = 0.1;
    alpha *= z;
    intt = 1;
    goto funct30;
funct49:
    ;
    debugmsg("inside 49");
    if(tot < aeps)
      goto funct92;
    goto funct50;
funct50:
    ;
    debugmsg("inside 50");
    alpha = tot;
    for(int i=1; i<=n; i++)
    {
      x[i] = w[i];
      w[ib+i] = g[i];
    }
    link = 2;
    if((idiff - 1) < 0)
      goto funct100;
    else
      if((idiff -1) == 0)
        goto funct100;
      else
        goto funct110;
funct52:
    ;
    debugmsg("inside 52");
    zz = alpha / (dgs - alpha*gs0);
    z = dgs * zz - 1;
    for(int i=1; i<=n; i++)
    {
      w[iu+i] = z * w[i] + g[i];
    }
    sig = 1.0 / (zz*dgs*dgs);
    goto funct70;
funct54:
    ;
    debugmsg("inside 54");
    if(ifn >= MaxFn)
      goto funct90;
    gys = 0;
    for(int i=1; i <= n; i++)
    {
      w[i] = w[ib+i];
      gys += g[i] * w[is+i];
    }
    df = ff - f;
    dgs = gys - gs0;
    if(dgs <= 0)
      goto funct20;
    link = 1;
    if(dgs + alpha*gs0 > 0)
      goto funct52;
    for(int i=1; i <= n; i++)
    {
      w[iu+i] = g[i] - w[i];
    }
    sig = 1 / (alpha*dgs);
    goto funct70;
funct60:
    ;
    debugmsg("inside 60");
    link = 2;
    for(int i=1; i<=n; i++)
    {
      w[iu+i] = w[i];
    }
    if(dgs + alpha * gs0 > 0)
      goto funct62;
    sig = 1.0 / gs0;
    goto funct70;
funct62:
    ;
    debugmsg("inside 62");
    sig = -1*zz;
    goto funct70;
funct70:
    ;
    debugmsg("inside 70");
    w[iv+1] = w[iu+1];
    for(int i=2; i<=n; i++)
    {
      ij = i;
      i1 = i - 1;
      z = w[iu+i];
      for(int j=1; j<=i1; j++)
      {
        z -= h[ij] * w[iv+j];
        ij += n - j;
      }
      w[iv+i] = z;
    }
    ij = 1;
    for(int i=1; i <= n; i++)
    {
      z = h[ij] + sig * w[iv+i] * w[iv+i];
      if(z <= 0)
        z = dmin;
      if(z < dmin)
        dmin = z;
      h[ij] = z;
      w[ib+i] = w[iv+i] * sig / z;
      sig -= w[ib+i] * w[ib+i] * z;
      ij += np -i;
    }
    ij = 1;
    for(int i=1; i<=n1; i++)
    {
      ij++;
      i1 = i + 1;
      for(int j=i1; j<=n; j++)
      {
        w[iu+j] -= h[ij] * w[iv+i];
        h[ij] += w[ib+i] * w[iu+j];
        ij++;
      }
    }
    if(link == 1)
      goto funct60;
    else
      goto funct20;
funct90:
    ;
    debugmsg("inside 90");
    iexit = 3;
    goto funct94;
funct92:
    ;
    debugmsg("inside 92");
    if(idiff == 2)
      goto funct94;
    idiff = 2;
    goto funct17;
funct94:
    ;
    debugmsg("inside 94");
    /////////////
    for(int i=1; i<n+1; i++)
      val[i-1] = x[i];
    ////////////
    if(PrintMode == 0)
      return;
    OutStream << "itn, ifn, iexit " << itn << " " << ifn << " " << iexit << std::endl;
    OutStream << "f = " << f << std::endl;
    for(int i=1; i<=n; i++)
      OutStream << x[i] << " ";
    OutStream << std::endl;
    for(int i=1; i<=n; i++)
      OutStream << g[i] << " ";
    OutStream << std::endl;
    return;
funct100:
    ;
    debugmsg("inside 100");
    for(int i=1; i<= n; i++)
    {
      z = HH * xm[i];
      w[i] += z;
      f1 = Funct(n,&w[1]);
      g[i] = (f1 - f) / z;
      w[i] -= z;
    }
    ifn += n;
    if(link == 1)
    {
      goto funct18;
    }
    else
    {
      goto funct54;
    }
funct110:
    ;
    debugmsg("inside 110");
    for(int i=1; i <= n; i++)
    {
      z = HH * xm[i];
      w[i] += z;
      f1 = Funct(n,&w[1]);
      debugmsg("f1 = " << f1);
      w[i] -= 2*z;
      f2 = Funct(n,&w[1]);
      debugmsg("f2 = " << f2);
      g[i] = (f1 - f2) / (2.0 * z);
      w[i] += z;
    }
    ifn += 2*n;
    if(link == 1)
      goto funct18;
    else
      goto funct54;
  }
};

#endif

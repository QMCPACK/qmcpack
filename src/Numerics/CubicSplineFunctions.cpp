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
    
    




template<class T>
T
CubicSpline<T>::operator()(T r)
{
  const T onesixth = 0.16666666666666666667e0;
  int klo,khi,k;
  T h,b,a;
  if(r>X[1] && r < X[Npt])
    //NR part
  {
    klo=1;
    khi=Npt;
    while (khi-klo > 1)
    {
      k=(khi+klo) >> 1;
      if (X[k] > r)
        khi=k;
      else
        klo=k;
    }
    h=X[khi]-X[klo];
    if (h == 0.0)
    {
      // ERRORMSG("splint called with x = %lf\n",x);
    };
    a=(X[khi]-r)/h;
    b=(r-X[klo])/h;
    return a*Y[klo]+b*Y[khi]+
           (a*(a*a-1.e0)*Y2[klo]+b*(b*b-1.0e0)*Y2[khi])*(h*h)*onesixth;
  }
  else
    // Extrapolation past endpoints
  {
    if(r<=X[1])
    {
      return Y[1] + Yp1 * (r - X[1]);
    }
    else
    {
      return Y[Npt] + Ypn * (r - X[Npt]);
    }
  }
}

template<class T>
T
CubicSpline<T>::operator()(T r, T& val)
{
  const T onesixth = 0.16666666666666666667e0;
  int klo,khi,k;
  T h,b,a;
  if(r>X[1] && r < X[Npt])
    //NR part
  {
    klo=1;
    khi=Npt;
    while (khi-klo > 1)
    {
      k=(khi+klo) >> 1;
      if (X[k] > r)
        khi=k;
      else
        klo=k;
    }
    h=X[khi]-X[klo];
    if (h == 0.0)
    {
      //ERRORMSG("splint called with r = " << r << std::endl);
    };
    a=(X[khi]-r)/h;
    b=(r-X[klo])/h;
    val=a*Y[klo]+b*Y[khi]+
        (a*(a*a-1.e0)*Y2[klo]+b*(b*b-1.0e0)*Y2[khi])*(h*h)*onesixth;
    val /= r;
    return -(1/h)*(Y[klo]-Y[khi] + ((h*h)*onesixth)*((3*a*a-1)*Y2[klo] -
                   (3*b*b-1)*Y2[khi]) );
  }
  else
    // Extrapolation past endpoints
  {
    if(r<=X[1])
    {
      val =(Y[1] + Yp1 * (r - X[1]))/r;
      return Yp1;
    }
    else
    {
      val= (Y[Npt] + Ypn * (r - X[Npt]))/r;
      return  Ypn;
    }
  }
}

template<class T>
T
CubicSpline<T>::operator()(T r, T& val, T& yp)
{
  const T onesixth = 0.16666666666666666667e0;
  int klo,khi,k;
  T h,b,a;
  T rinv = 1.0/r;
  if(r>X[1] && r < X[Npt])
    //NR part
  {
    klo=1;
    khi=Npt;
    while (khi-klo > 1)
    {
      k=(khi+klo) >> 1;
      if (X[k] > r)
        khi=k;
      else
        klo=k;
    }
    h=X[khi]-X[klo];
    if (h == 0.0)
    {
      //ERRORMSG("splint called with r = " << r << std::endl);
    };
    T hinv = 1/h;
    a=(X[khi]-r)*hinv;
    b=(r-X[klo])*hinv;
    val=a*Y[klo]+b*Y[khi]+
        (a*(a*a-1.e0)*Y2[klo]+b*(b*b-1.0e0)*Y2[khi])*(h*h)*onesixth;
    val *= (rinv*rinv);
    yp = -rinv*hinv*(Y[klo]-Y[khi] +
                     (h*h*onesixth)*((3*a*a-1)*Y2[klo] -(3*b*b-1)*Y2[khi]) );
    return a*Y2[klo]+b*Y2[khi];
  }
  else
    // Extrapolation past endpoints
  {
    if(r<=X[1])
    {
      val =(Y[1] + Yp1 * (r - X[1]))*rinv*rinv;
      yp = Yp1*rinv;
      return 0.0e0;
    }
    else
    {
      val= (Y[Npt] + Ypn * (r - X[Npt]))*rinv*rinv;
      yp =  Ypn*rinv;
      return 0.0e0;
    }
  }
}


template<class T>
void CubicSpline<T>::spline()
{
  //This is just a modified version of the Numerical Recipes routine spline,
  //which generates spline coefficients given a set of spline parameters.
  int i,k;
  T p,qn,sig,un;
  std::vector<T> u(Npt+1);
  if (Yp1 > 0.99e30)
    Y2[1]=u[1]=0.0;
  else
  {
    Y2[1] = -0.5;
    u[1]=(3.0/(X[2]-X[1]))*((Y[2]-Y[1])/(X[2]-X[1])-Yp1);
  }
  for (i=2; i<=Npt-1; i++)
  {
    sig=(X[i]-X[i-1])/(X[i+1]-X[i-1]);
    p=sig*(Y2[i-1])+2.0;
    Y2[i]=(sig-1.0)/p;
    u[i]=(Y[i+1]-Y[i])/(X[i+1]-X[i]) - (Y[i]-Y[i-1])/(X[i]-X[i-1]);
    u[i]=(6.0*u[i]/(X[i+1]-X[i-1])-sig*u[i-1])/p;
  }
  if (Ypn > 0.99e30)
    qn=un=0.0;
  else
  {
    qn=0.5;
    un=(3.0/(X[Npt]-X[Npt-1]))*(Ypn-(Y[Npt]-Y[Npt-1])/(X[Npt]-X[Npt-1]));
  }
  Y2[Npt]=(un-qn*u[Npt-1])/(qn*(Y2[Npt-1])+1.0);
  for(k=Npt-1; k>=1; k--)
    Y2[k]= (Y2[k]*(Y2[k+1]))+u[k];
}

//
// Regular spacing
//
template<class T>
T
RegCubicSpline<T>::operator()(T r)
{
  if(r<Xmin)
  {
    return Y[1] + Yp1 * (r - Xmin);
  }
  else
    if(r<Xmax)
    {
      int i = static_cast<int>((r-Xmin)*DxInv+1);
      int ip = i+1;
      T rhi = static_cast<T>(i)*Dx+Xmin;
      T a=(rhi-r)*DxInv;
      T b=(r-rhi+Dx)*DxInv;
      return
        a*Y[i]+b*Y[ip]+(a*(a*a-1.e0)*Y2[i]+b*(b*b-1.0e0)*Y2[ip])*h2over6;
    }
    else
    {
      return Y[Npt] + Ypn * (r - Xmax);
    }
}

template<class T>
T
RegCubicSpline<T>::operator()(T r, T& val)
{
  if(r<Xmin)
  {
    val =(Y[1] + Yp1 * (r - Xmin))/r;
    return Yp1;
  }
  else
    if(r<Xmax)
    {
      int i = static_cast<int>((r-Xmin)*DxInv+1);
      int ip = i+1;
      T rhi = static_cast<T>(i)*Dx+Xmin;
      T a=(rhi-r)*DxInv;
      T b=(r-rhi+Dx)*DxInv;
      val = (a*Y[i]+b*Y[ip]+(a*(a*a-1.e0)*Y2[i]+b*(b*b-1.0e0)*Y2[ip])*h2over6)/r;
      return -DxInv*(Y[i] -Y[ip] + h2over6*((3*a*a-1)*Y2[i]-(3*b*b-1)*Y2[ip]));
    }
    else
    {
      val= (Y[Npt] + Ypn * (r - Xmax))/r;
      return  Ypn;
    }
}

template<class T>
T
RegCubicSpline<T>::operator()(T r, T& val, T& yp)
{
  T rinv = 1.0/r;
  if(r<Xmin)
  {
    val =(Y[1] + Yp1 * (r - Xmin))*rinv*rinv;
    yp = Yp1*rinv;
    return 0;
  }
  else
    if(r<Xmax)
    {
      int i = static_cast<int>((r-Xmin)*DxInv+1);
      int ip = i+1;
      T rhi = static_cast<T>(i)*Dx+Xmin;
      T a=(rhi-r)*DxInv;
      T b=(r-rhi+Dx)*DxInv;
      val = (a*Y[i]+b*Y[ip]
             +(a*(a*a-1.e0)*Y2[i]+b*(b*b-1.0e0)*Y2[ip])*h2over6)*rinv*rinv;
      yp = -rinv*DxInv*
           (Y[i] -Y[ip] + h2over6*((3*a*a-1)*Y2[i]-(3*b*b-1)*Y2[ip]));
      return a*Y2[i]+b*Y2[ip];
    }
    else
    {
      val= (Y[Npt] + Ypn * (r - Xmax))*rinv*rinv;
      yp = Ypn*rinv;
      return 0;
    }
}

template<class T>
void
RegCubicSpline<T>::spline()
{
  //This is just a modified version of the Numerical Recipes routine spline,
  //which generates spline coefficients given a set of spline parameters.
  int i,k;
  T p,qn,sig,un;
  std::vector<T> u(Npt+1);
  if (Yp1 > 0.99e30)
    Y2[1]=u[1]=0.0;
  else
  {
    Y2[1] = -0.5;
    u[1]=(3.0*DxInv)*((Y[2]-Y[1])*DxInv-Yp1);
  }
  sig=0.5;
  for (i=2; i<=Npt-1; i++)
  {
    p=sig*(Y2[i-1])+2.0;
    Y2[i]=(sig-1.0)/p;
    u[i]=(Y[i+1]+Y[i-1]-2*Y[i])*DxInv;
    u[i]=(3.0*u[i]*DxInv-sig*u[i-1])/p;
  }
  if (Ypn > 0.99e30)
    qn=un=0.0;
  else
  {
    qn=0.5;
    un=(3.0*DxInv)*(Ypn-(Y[Npt]-Y[Npt-1])*DxInv);
  }
  Y2[Npt]=(un-qn*u[Npt-1])/(qn*(Y2[Npt-1])+1.0);
  for(k=Npt-1; k>=1; k--)
    Y2[k]= (Y2[k]*(Y2[k+1]))+u[k];
}


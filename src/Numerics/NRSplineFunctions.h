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
    
    



#ifndef QMCPLUSPLUS_NR_CUBICSPLINE_H
#define QMCPLUSPLUS_NR_CUBICSPLINE_H

/**template function: converted from Numerical Recipe spline.c
 *note that the range of data is [0,n) instead of [1,n]
 */
template<typename Tg, typename T>
inline void
NRCubicSpline(const Tg* x, const T* y, int n, T yp1, T ypn, T* y2)
{
  int i,k;
  T p,qn,sig,un;
  std::vector<T> u(n);
  if (yp1 > 0.99e30)
    y2[0]=u[0]=0.0;
  else
  {
    y2[0] = -0.5;
    u[0]=(3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
  }
  for (i=1; i<n-1; i++)
  {
    sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
    p=sig*y2[i-1]+2.0;
    y2[i]=(sig-1.0)/p;
    u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
    u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
  }
  if (ypn > 0.99e30)
    qn=un=0.0;
  else
  {
    qn=0.5;
    un=(3.0/(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
  }
  y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.0);
  for (k=n-2; k>=0; k--)
    y2[k]=y2[k]*y2[k+1]+u[k];
}

#include <climits>

/**template function:
 *note that the range of data is [0,n) instead of [1,n]
 * solve the linear system for the first and second derivatives
 *  Copy of routine from: "One dimensional spline interpolation algorithms" by Helmuth Spath
 * QUINAT routine
 */
// morales:
// FIX FIX FIX: right now it uses zero derivatives at both ends...
template<typename Tg, typename T>
inline void
QuinticSplineSolve(int N, const Tg* X, T* Y, T* B, T* C, T* D, T* E, T* F)
{
  int M;
  // mmorales: had issues setting std::numeric_limits<T>::epsilon() with gcc, FIX later
  T eps = 1.0e-9;
  //eps= std::numeric_limits<double>::epsilon();
  T Q,R,Q2,Q3,R2,QR,P,P2,PQ,PR,PQQR,B1,V,TT,S,U,P3;
  M=N-2;
  Q=X[1]-X[0];
  R=X[2]-X[1];
  Q2=Q*Q;
  R2=R*R;
  QR=Q+R;
  D[0]=E[0]=D[1]=0.0;
  if(std::abs(Q) > eps)
    D[1]=6.0*Q*Q2/(QR*QR);
  for(int i=1; i<M; i++)
  {
    P=Q;
    Q=R;
    R=X[i+2]-X[i+1];
    P2=Q2;
    Q2=R2;
    R2=R*R;
    PQ=QR;
    QR=Q+R;
    if(std::abs(Q) <= eps)
    {
      D[i+1]=E[i]=F[i-1]=0.0;
    }
    else
    {
      Q3=Q2*Q;
      PR=P*R;
      PQQR=PQ*QR;
      D[i+1]=6.0*Q3/(QR*QR);
      D[i]  += 2.0*Q*( 15.0*PR*PR + (P+R)*Q*(20.0*PR+7.0*Q2) + Q2*(8.0*(P2+R2)
                       + 21.0*PR + 2.0*Q2) )/(PQQR*PQQR);
      D[i-1] += 6.0*Q3/(PQ*PQ);
      E[i] = Q2*(P*QR+3.0*PQ*(QR+2.0*R))/(PQQR*QR);
      E[i-1] += Q2*(R*PQ+3.0*QR*(PQ+2.0*P))/(PQQR*PQ);
      F[i-1] = Q3/PQQR;
    }
  }
  if(std::abs(R) > eps)
    D[N-3] += 6.0*R*R2/(QR*QR);
  for(int i=1; i<N; i++)
  {
    if(std::abs(X[i]-X[i-1]) < eps)
    {
      B[i]=Y[i];
      Y[i]=Y[i-1];
    }
    else
    {
      B[i] = (Y[i]-Y[i-1])/(X[i]-X[i-1]);
    }
  }
  for(int i=2; i<N; i++)
  {
    if(std::abs(X[i]-X[i-2]) < eps)
    {
      C[i]=B[i]*0.5;
      B[i]=B[i-1];
    }
    else
    {
      C[i] = (B[i]-B[i-1])/(X[i]-X[i-2]);
    }
  }
  // assume N > 5
  P=C[0]=E[N-3]=F[0]=F[N-4]=F[N-3]=0.0;
  C[1]=C[3]-C[2];
  D[1]=1.0/D[1];
  for(int i=2; i<M; i++)
  {
    Q=D[i-1]*E[i-1];
    D[i]=1.0/(D[i]-P*F[i-2]-Q*E[i-1]);
    E[i]-=Q*F[i-1];
    C[i]=C[i+2]-C[i+1]-P*C[i-2]-Q*C[i-1];
    P=D[i-1]*F[i-1];
  }
  C[N-2]=0.0;
  C[N-1]=0.0;
  for(int i=N-3; i>0; i--)
    C[i] = (C[i]-E[i]*C[i+1]-F[i]*C[i+2])*D[i];
  M=N-1;
  Q=X[1]-X[0];
  R=X[2]-X[1];
  B1=B[1];
  Q3=Q*Q*Q;
  QR=Q+R;
  if(std::abs(QR) < eps)
  {
    V=TT=0.0;
  }
  else
  {
    V=C[1]/QR;
    TT=V;
  }
  F[0]=0.0;
  if(std::abs(Q) > eps )
    F[0]=V/Q;
  for(int i=1; i<M; i++)
  {
    P=Q;
    Q=R;
    R=0.0;
    if(i != (M-1))
      R=X[i+2]-X[i+1];
    P3=Q3;
    Q3=Q*Q*Q;
    PQ=QR;
    QR=Q+R;
    S=TT;
    TT=0;
    if(std::abs(QR) > eps )
      TT=(C[i+1]-C[i])/QR;
    U=V;
    V=TT-S;
    if(std::abs(PQ) < eps )
    {
      C[i]=C[i-1];
      D[i]=E[i]=F[i]=0.0;
    }
    else
    {
      F[i] = F[i-1];
      if(std::abs(Q) > eps )
        F[i]=V/Q;
      E[i]=5.0*S;
      D[i]=10.0*(C[i]-Q*S);
      C[i]=D[i]*(P-Q)+(B[i+1]-B[i]+(U-E[i])*P3
                       -(V+E[i])*Q3)/PQ;
      B[i]=(P*(B[i+1]-V*Q3)+Q*(B[i]-U*P3))/PQ
           -P*Q*(D[i]+E[i]*(Q-P));
    }
  }
  P=X[1]-X[0];
  S=F[0]*P*P*P;
  E[0]=D[0]=E[N-1]=D[N-1]=0.0;
  C[0]=C[1]-10.0*S;
  B[0]=B1-(C[0]+S)*P;
  Q=X[N-1]-X[N-2];
  TT=F[N-2]*Q*Q*Q;
  C[N-1]=C[N-2]+10.0*TT;
  B[N-1]=B[N-1]+(C[N-1]-TT)*Q;
}


#endif

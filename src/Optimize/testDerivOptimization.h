//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#ifndef QMCPLUSPLUS_TESTDERIVOPTIMIZATION_H
#define QMCPLUSPLUS_TESTDERIVOPTIMIZATION_H

#include "Optimize/OptimizeBase.h"
#include "OhmmsData/ParameterSet.h"

#include <time.h>


template<class T>
class testDerivOptimization: public MinimizerBase<T>
{
public:

  typedef T Return_t;
  typedef typename MinimizerBase<T>::ObjectFuncType ObjectFuncType;
  using MinimizerBase<T>::msg_stream;


  ObjectFuncType* TargetFunc;
  int NumParams;
  Return_t Displacement;
  Return_t CostTol;
  Return_t GradTol;
  Return_t deltaG;
  int CurStep;

  int a_its ;
  int a_itmax ;
  int a_rich ;
  int a_verbose ;
  Return_t a_stepmax ;
  int a_linmin_maxits ;
  Return_t a_linmin_g1 ;
  Return_t a_linmin_g2 ;
  Return_t a_linmin_g3 ;
  Return_t a_lastx     ;
  Return_t a_lastx_default ;
  int CG_ortho;
  int xycleanup;
  Return_t a_gtyp ;
  std::vector<Return_t>  a_pt ,  a_gx ,  a_gy , a_gunused ;
  std::vector<Return_t>  a_xi , a_g , a_h, Parms;
  int a_restart ;
  bool Failed_Last;
  int xybisect;
  Return_t Gfactor;
  std::vector< std::vector<Return_t> > Past_a_xi, Past_Parms, Past_a_h;
  std::vector< std::vector<Return_t> > PastGradients;
  std::vector< std::vector<Return_t> > PastGradientsParameterPoints;

  testDerivOptimization(ObjectFuncType* atarget=0):  TargetFunc(atarget),
    Displacement(0),xybisect(0), Failed_Last(false),
    CostTol(1.e-3),GradTol(1.e-2), CG_ortho(1), a_linmin_maxits(5),
    a_lastx(1e-2), a_lastx_default(1e-2), a_linmin_g1(2), a_linmin_g2(1.0),
    a_linmin_g3(0.5), a_restart(0), xycleanup(4), Gfactor(3), deltaG(1e-5)
  {
    if (atarget)
      setTarget(atarget);
  }

  ~testDerivOptimization() {}

  void setTarget(ObjectFuncType* fn)
  {
    TargetFunc=fn;
    NumParams=TargetFunc->NumParams();
    resizeAllArray(NumParams);
    for (int i=0; i<NumParams; i++)
      Parms[i]= TargetFunc->Params(i);
    //max CG steps is the dimensionality of the problem
    CG_ortho = std::min(CG_ortho,NumParams);
    //if CG=0 we are using steepest descent. if <0 set to pure CG
    if (CG_ortho<0)
      CG_ortho=NumParams;
  }

  void resizeAllArray(int newSize)
  {
    a_g.resize(newSize,0.0);
    a_gunused.resize(newSize,0.0);
    a_h.resize(newSize,0.0);
    a_xi.resize(newSize,0.0);
    a_pt.resize(newSize,0.0);
    a_gx.resize(newSize,0.0);
    a_gy.resize(newSize,0.0);
    Parms.resize(newSize,0.0);
  }

  bool optimize(ObjectFuncType* fn)
  {
    setTarget(fn);
    return optimize();
  }

  bool optimize()
  {
    reset();
    bool successM = jmacOpt();
    //make sure wave function has the right parameters for next runs
    for (int i=0; i<NumParams; i++)
      TargetFunc->Params(i)=Parms[i];
    TargetFunc->Report();
    return successM;
  }

  bool get(std::ostream&) const;

  bool put(std::istream&);

  /**  Parse the xml file for parameters
   * @param cur current xmlNode
   * @param a_itmax maximum number of CG iterations
   * @param Displacement used for finite difference cost function
   * @param CG_ortho number of previous search directions we want to maintain orthogonality to.
   * @param a_rich 0=used approximate gradient at new location, 1=calculate new gradient at line minimum
   * @param xybisect Number of times to use bisection before linear interpolation for line minimization minimum finding routine.
   * @param Gfactor max allowed increase in Cost function gradient. Avoids falling in holes.
   * @param xycleanup max steps to tighten the line search
   * @param GradTol gradient to quit line minimization
   * @param a_verbose 0=quiet, 1=normal, 2=chatty, 3+=various degrees of debug (loud!)
   * @param a_lastx_default default step size. Only has a transitive effect.
   * @param a_linmin_maxits Maximum number of steps to try to find a minimum along a_h direction
   * @param deltaG if change in gradient is less than this relax the CG_ortho constraints
   *
   */

  bool put(xmlNodePtr cur)
  {
    ParameterSet p;
    p.add(a_itmax,"max_steps","none");
    p.add(a_itmax,"maxSteps","none");
    p.add(Displacement,"epsilon","none");
    p.add(CG_ortho,"length_CG","none");
    p.add(a_rich,"rich","none");
    p.add(xybisect,"xybisect","none");
    p.add(Gfactor,"Gfactor","none");
    p.add(xycleanup,"xypolish","none");
    p.add(CostTol,"tolerance","none");
    p.add(GradTol,"tolerance_g","none");
    p.add(GradTol,"toleranceG","none");
    p.add(deltaG,"tolerance_cg","none");
    p.add(deltaG,"toleranceCG","none");
    p.add(a_verbose,"verbose","none");
    p.add(a_lastx_default,"stepsize","none");
    p.add(a_lastx_default,"stepSize","none");
    p.add(a_linmin_maxits,"max_linemin","none");
    p.add(a_linmin_maxits,"maxLinemin","none");
    p.put(cur);
    a_lastx=a_lastx_default;
    return true;
  }

  void reset()
  {
    Past_a_xi.clear();
    Past_Parms.clear();
    Past_a_h.clear();
    PastGradients.clear();
    PastGradientsParameterPoints.clear();
    Failed_Last=false;
  }

  Return_t func(std::vector<Return_t> _p)
  {
    for (int i=0; i<NumParams; ++i)
      TargetFunc->Params(i) = _p[i] ;
    return TargetFunc->Cost();
  }


  void dfunc(std::vector<Return_t> RT, std::vector<Return_t> &FG)
  {
    ///To test we simply output the analytic and numeric gradients of the cost function. Make sure they agree.
    std::vector<Return_t> Dummy(FG);
    TargetFunc->GradCost( Dummy, RT, 1e-5);
    TargetFunc->GradCost(FG, RT, 0.0);
    PastGradients.push_back(FG);
    PastGradientsParameterPoints.push_back(RT);
    std::cout <<"Numeric            Analytic       Percent"<< std::endl;
    for(int k=0; k<NumParams; k++)
    {
      if (Dummy[k]!=0)
        std::cout <<Dummy[k]<<"  "<<FG[k]<<"  "<< 100*(Dummy[k]-FG[k])/Dummy[k]<< std::endl;
      else
        std::cout <<Dummy[k]<<"  "<<FG[k]<<"   inf"<< std::endl;
    }
    std::cout << std::endl;
  }

  /** Main Optimization routine
   *  adapted from the macoptII optimizer by David MacKay
   *  http://www.inference.phy.cam.ac.uk/mackay/c/macopt.html#license
   *
   * Several safeguards for stability.
   */
  bool jmacOpt()
  {
    Return_t gg(0.0) , gam(0.0) , dgg(0.0) ;
    Return_t step , tmpd ;
    dfunc(Parms, a_xi);
    Past_a_xi.push_back(a_xi);
    Past_Parms.push_back(Parms);
    restart(1);
    Past_a_h.push_back(a_h);
    gg = 0.0;
    for (int j = 0 ; j < NumParams ; j ++)
      gg += a_xi[j]*a_xi[j];
    a_gtyp = sqrt(gg / (1.0*NumParams));
    for (int a_its = 0 ; ((a_its <= a_itmax)&&(TargetFunc->IsValid)&(gg>CostTol)) ; ++a_its)
    {
      if (a_verbose > 0)
        printf("mac_it %d of %d : gg = %6.3g tol = %6.3g: \n", a_its , a_itmax , gg , CostTol) ;
      step = maclinminII(Parms);
      if (step<0.0)
        return false;
      Return_t GGnew = 0.0;
      for (int j = 0 ; j < NumParams ; j ++)
        GGnew += a_xi[j]*a_xi[j];
      if (GGnew>Gfactor*gg)
      {
        if (a_verbose>1)
        {
          std::cout <<"G grew too much "<<GGnew<<" > "<<Gfactor<<"*"<<gg<< std::endl;
          printf("mac_it %d of %d : gg = %6.3g tol = %6.3g: \n", a_its , a_itmax , gg , CostTol) ;
          if (a_verbose>3)
          {
            int PastLines = Past_a_h.size()-1;
            std::cout <<PastLines<< std::endl;
            for (int ih=PastLines; ((ih>-1) & (ih>(PastLines-CG_ortho))) ; ih--)
            {
              std::cout <<"Past_a_h ["<<ih<<"] ";
              for (int j = 0 ; j < NumParams ; j ++)
                std::cout <<Past_a_h[ih][j]<<" ";
              std::cout << std::endl;
            }
          }
        }
        Past_Parms.pop_back();
        Past_a_xi.pop_back();
        Past_a_h.pop_back();
        if (Past_Parms.size()>0)
        {
          Parms= Past_Parms[Past_Parms.size()-1];
          a_xi = Past_a_xi[Past_a_xi.size()-1];
          gg=0.0;
          for (int j = 0 ; j < NumParams ; j ++)
            gg += a_xi[j]*a_xi[j];
          a_gtyp = sqrt(gg / (1.0*NumParams));
          a_h = a_xi;
          for (int j = 0 ; j < NumParams ; j ++)
            a_h[j] *= -1.0;
          Past_a_h.push_back(a_h);
        }
        else
          return false;
        //dfunc( Parms , a_xi  ) ;
        if (Failed_Last)
          return false;
        else
          Failed_Last=true;
      }
      else
      {
        //here we relax the CG towards steepest descent. Is this the right thing to do?
        if ((deltaG > std::abs(gg-GGnew)) & (CG_ortho>0))
        {
          --CG_ortho;
          if (a_verbose>3)
          {
            std::cout <<"Relaxing Orthogonality condition. Rotating to Steepest Descent"<< std::endl;
          }
        }
        Past_Parms.push_back(Parms);
        Failed_Last=false;
        gg=GGnew;
      }
      if (a_restart)
      {
        restart(0) ;
      }
      if (a_verbose>3)
      {
        int PastLines = Past_a_h.size()-1;
        for (int ih=PastLines; ((ih>-1) & (ih>(PastLines-CG_ortho))) ; ih--)
        {
          std::cout <<"Past_a_h ["<<ih<<"] ";
          for (int j = 0 ; j < NumParams ; j ++)
            std::cout <<Past_a_h[ih][j]<<" ";
          std::cout << std::endl;
          std::cout <<"Past_a_xi ["<<ih<<"] ";
          for (int j = 0 ; j < NumParams ; j ++)
            std::cout <<Past_a_xi[ih][j]<<" ";
          std::cout << std::endl;
        }
      }
      if (!Failed_Last)
      {
        std::vector<Return_t> tmpA_H(a_xi);
        for (int j = 0 ; j < NumParams ; j ++)
          tmpA_H[j] *= -1.0;
        int PastLines = Past_a_h.size()-1;
        for (int ih=PastLines; ((ih>-1) & (ih>(PastLines-CG_ortho))) ; --ih)
        {
          Return_t tmpDGG(0.0), tempGG(0.0);
          for (int j = 0 ; j < NumParams ; j ++)
          {
            tmpDGG += Past_a_h[ih][j] * tmpA_H[j] ;
            tempGG += Past_a_h[ih][j] * Past_a_h[ih][j];
          }
          Return_t Gamm = (tmpDGG / tempGG);
          for (int j = 0 ; j < NumParams ; j ++)
            tmpA_H[j] -= Gamm * Past_a_h[ih][j];
        }
        a_h=tmpA_H;
        Past_a_h.push_back(tmpA_H);
        if (a_verbose>3)
        {
          PastLines = Past_a_h.size()-1;
          std::cout <<"Inner Products for line searches:"<< std::endl;
          for (int ih=PastLines; ((ih>-1) & (ih>(PastLines-CG_ortho-1))) ; --ih)
          {
            for (int ih2=ih; ((ih2>-1) & (ih2>(PastLines-CG_ortho-1))) ; --ih2)
            {
              std::cout <<"Between "<<ih<<"  "<<ih2<<"  ";
              Return_t tempGG(0.0);
              for (int j = 0 ; j < NumParams ; j ++)
              {
                tempGG += Past_a_h[ih][j] * Past_a_h[ih2][j];
              }
              std::cout <<tempGG<< std::endl;
            }
          }
        }
      }
    }
    if (a_verbose > 0)
      printf("FINAL : gg = %6.3g tol = %6.3g: \n", gg , CostTol) ;
    return ((TargetFunc->IsValid) && (gg <= CostTol));
  }


  Return_t maclinminII(std::vector<Return_t> &p)
  {
    Return_t y(0.0) ;
    Return_t  t(0.0) ;
    Return_t step(0.0) , tmpd(0.0) ;
    tmpd = lineProduct(Parms , a_gy , 0.0) ;
    if (tmpd > 0.0)
    {
      a_restart = 1 ;
      return 0.0 ;
    }
    int its(0);
    Return_t x = a_lastx / a_gtyp ;
    Return_t s = lineProduct(Parms , a_gx , x) ;
    if (s < 0 && TargetFunc->IsValid)
    {
      do
      {
        y = x * a_linmin_g1 ;
        t = lineProduct(Parms , a_gy , y) ;
        if ((t >= 0.0) || (!TargetFunc->IsValid))
        {
          if (a_verbose > 1)
            printf("s < 0: s = %6.3g: t = %6.3g; x = %6.3g y = %6.3g\n",s, t , x , y);
          break ;
        }
        if (t==t)
        {
          if (a_verbose > 1)
            printf("s < 0: s = %6.3g: t = %6.3g; x = %6.3g y = %6.3g\n",s, t , x , y);
          x = y ;
          s = t ;
          a_gunused = a_gx ;
          a_gx = a_gy ;
          a_gy = a_gunused ;
        }
        else
        {
          y=x;
          t=s;
          its = a_linmin_maxits;
        }
        its++ ;
      }
      while ((its <= a_linmin_maxits)&&(TargetFunc->IsValid)) ;
    }
    else
      if (s > 0 && TargetFunc->IsValid)
      {
        do
        {
          y = x * a_linmin_g3 ;
          t = lineProduct(Parms , a_gy , y) ;
          if ((t <= 0.0) || (!TargetFunc->IsValid))
          {
            if (a_verbose > 1)
              printf("s > 0: s = %6.3g: t = %6.3g; x = %6.3g y = %6.3g\n",s, t , x , y);
            break ;
          }
          if (t==t)
          {
            if (a_verbose > 1)
              printf("s > 0: s = %6.3g: t = %6.3g; x = %6.3g y = %6.3g\n",s, t , x , y);
            x = y ;
            s = t ;
            a_gunused = a_gx ;
            a_gx = a_gy ;
            a_gy = a_gunused ;
          }
          else
          {
            y=x;
            t=s;
            its = a_linmin_maxits;
          }
          its ++ ;
        }
        while ((its <= a_linmin_maxits)&&(TargetFunc->IsValid)) ;
      }
      else
      {
        if (a_verbose > 1)
          printf("s = 0: s = %6.3g: t = %6.3g; x = %6.3g y = %6.3g\n",s, t , x , y);
        t = 1.0 ;
        y = x;
      }
    Return_t u(0.0);
    Return_t m(0.5*(x+y));
    int xyit(0);
    if ((xyit<xycleanup) && (std::abs(s-t)>GradTol) && (s*t<0.0) &&(TargetFunc->IsValid))
    {
      int XYBisectCounter=xybisect;
      do
      {
        if (XYBisectCounter)
        {
          if (a_verbose > 1)
            printf("XYBisectCounter %d of %d  :", XYBisectCounter ,xybisect) ;
          m = 0.5*(x+y) ;
          XYBisectCounter--;
        }
        else
        {
          Return_t ms(std::abs(s)), mt(std::abs(t));
          m= (ms*y + mt*x)/(ms+mt);
        }
        u=lineProduct(Parms , a_gunused , m) ;
        if (a_verbose > 1)
          printf("xyit %d of %d : u = %6.3g m = %6.3g: \n", xyit ,xycleanup, u , m) ;
        if (u==u)
        {
          if (u*s >= 0.0  && TargetFunc->IsValid)
          {
            s=u;
            x=m;
            a_gx = a_gunused;
            if (a_verbose > 1)
              printf("s = %6.3g: t = %6.3g; x = %6.3g y = %6.3g\n",s, t , x , y);
          }
          else
            if (u*t >= 0.0 && TargetFunc->IsValid)
            {
              t=u;
              y=m;
              a_gy = a_gunused ;
              if (a_verbose > 1)
                printf("s = %6.3g: t = %6.3g; x = %6.3g y = %6.3g\n",s, t , x , y);
            }
        }
        else
          if (std::abs(t)>std::abs(s))
          {
            s=t;
            x=y;
            xyit=xycleanup;
          }
          else
          {
            t=s;
            y=x;
            xyit=xycleanup;
          }
        xyit++;
      }
      while ((TargetFunc->IsValid) && (std::abs(s-t)>GradTol) && (xyit<xycleanup)) ;
    }
    else
      if ((s*t>0.0) || (!TargetFunc->IsValid))
      {
        return -1.0;
      }
    if (TargetFunc->IsValid)
    {
      Return_t ms(std::abs(s)), mt(std::abs(t));
      m= (ms*y + mt*x)/(ms+mt);
      step = 0.0;
      for (int i = 0 ; i < NumParams ; i ++)
      {
        tmpd = m * a_h[i] ;
        Parms[i] += tmpd ;
        step += std::abs(tmpd) ;
        a_xi[i] = (ms/(ms+mt)) * a_gy[i] + (mt/(ms+mt)) * a_gx[i] ;
      }
      if (a_rich || a_restart)
        dfunc(Parms , a_xi);
      Past_a_xi.push_back(a_xi);
      a_lastx = m * a_linmin_g2 *  a_gtyp ;
      for (int j=0; j<NumParams; j++)
        TargetFunc->Params(j)=Parms[j];
      TargetFunc->Report();
    }
    return (step / (1.0*NumParams)) ;
  }


  Return_t lineProduct(std::vector<Return_t> p , std::vector<Return_t> &gy , Return_t y)
  {
    Return_t s = 0.0 ;
    for (int i = 0 ; i < NumParams ; i ++)
      a_pt[i] = p[i] + y * a_h[i] ;
    dfunc(a_pt , gy) ;
    for (int i = 0 ; i < NumParams ; i ++)
      s += gy[i] * a_h[i] ;
    return s ;
  }

  void restart(int start)
  {
    if (start == 0)
      a_lastx = a_lastx_default;
    for (int j = 0 ; j < NumParams ; j ++)
      a_h[j] = -a_xi[j];
    a_restart = 0;
  }

};

#endif

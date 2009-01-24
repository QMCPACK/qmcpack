
#ifndef QMCPLUSPLUS_MACOPT_OPTIMIZATION_H
#define QMCPLUSPLUS_MACOPT_OPTIMIZATION_H

#include "Optimize/OptimizeBase.h"
#include "OhmmsData/ParameterSet.h"

template<class T>
class MacOptimization: public MinimizerBase<T> {
  public:
    
    typedef T Return_t;
    typedef typename MinimizerBase<T>::ObjectFuncType ObjectFuncType;
    using MinimizerBase<T>::msg_stream;
    
    /** Function to be optimized.  */
  ObjectFuncType* TargetFunc;
  
  /** number of minimization steps */
  int NumSteps;
  /** number of parameters to vary */
  int NumParams;
  /** displacement to evaluate the gradients numerically */
  Return_t Displacement;
  /** tolerance the converged value */
  Return_t CostTol;
  /** tolerance for the converged gamma = (g*g-g*cg)/g0*g0 */
  Return_t GammaTol;
  /** tolerance for the converged gradients ||g|| < GradTol */
  Return_t GradTol;
  /** tolerance for the converged gradient component max|g[i]| < GradMaxTol */
  Return_t GradMaxTol;
  Return_t LambdaMax;
  
  ///Parameters from MacOpt++

  bool RestartCG;
  int CurStep;
  
  Return_t gdotg, gdotg0, gdoth, gamma;
  Return_t a_tol ;    /* convergence declared when the gradient vector is smaller
  in magnitude than this, or when the mean absolute 
  step is less than this (see above) */
  Return_t a_grad_tol_tiny ; /* if gradient is less than this, we definitely 
  stop, even if we are not using a gradient 
  tolerance */
  Return_t a_step_tol_tiny ; /* if step is less than this, we stop, even if 
  we are not using a step tolerance */
  int a_end_if_small_step ; /* defines the role of tol -- alternative is
  end_on_small_grad */
  int a_its ;               /* number of its */
  int a_itmax ;             /* max */
  int a_rich ; /* whether to do the extra gradient evaluation at the beginning 
  of each new line min */
  int a_verbose ; 
  Return_t a_stepmax ;        /* largest step permitted (not used in macopt) */

  int a_linmin_maxits ;     /* in maclinmin */
  Return_t a_linmin_g1 ;      /* factors for growing and shrinking the interval */
  Return_t a_linmin_g2 ;
  Return_t a_linmin_g3 ;
  Return_t a_lastx     ;      /* keeps track of typical step length */
  Return_t a_lastx_default ;  /* if maclinmin is reset, lastx is set to this */

  int xycleanup;
  Return_t a_gtyp ; /* stores the rms gradient for linmin */
  vector<Return_t>  a_pt ,  a_gx ,  a_gy , a_gunused ;
  vector<Return_t>  a_xi , a_g , a_h, Parms ; 
  int a_restart ;           /* whether to restart macopt - fresh cg directions */
  //If we want to do a steepest descent minimization
  bool SD;








/** constructor
* @param atarget target function to optimize
*/
  MacOptimization(ObjectFuncType* atarget=0):  TargetFunc(atarget), RestartCG(true),
  NumSteps(100),Displacement(1e-6),
  CostTol(1.e-6),GradTol(1.e-6),GammaTol(1.e-7) // GammaTol was set to 1e-4 originally
  {
    if (atarget) setTarget(atarget);
    a_end_if_small_step = 0 ; /* Change this to 0/1 if you prefer */
  a_grad_tol_tiny = 1e-16 ; /* Probably not worth fiddling with */
  a_step_tol_tiny = 0.0 ;   /* Probably not worth fiddling with */
  a_linmin_maxits = 5 ;    /* Probably not worth fiddling with */
  a_lastx = 1e-2 ;          /* only has a transient effect      */
  a_lastx_default = 1e-5 ;  /* -- defines typical distance in parameter
  space at which the line minimum is expected;
  both these should be set. the default is 
  consulted if something goes badly wrong and 
  a reset is demanded. */

/* don't fiddle with the following, unless you really mean it */
  a_linmin_g1 = 1.5 ; 
  a_linmin_g2 = 1.0 ; 
  a_linmin_g3 = 0.5 ; 
  a_restart = 0 ; 
  xycleanup=3;
  } 
  
  /** destructor */
  ~MacOptimization() {}
  
  /** set the target function
  * @param fn target function to optimize
  *
  * Allocate the internal data.
  */
  void setTarget(ObjectFuncType* fn) {
    TargetFunc=fn;
    NumParams=TargetFunc->NumParams();
    resizeAllArray(NumParams);
    for(int i=0;i<NumParams; i++) 
    {
      Parms[i]= TargetFunc->Params(i);
    }
  }
  
  void resizeAllArray(int newSize)
  {
      a_g.resize(newSize,0.0);
      a_gunused.resize(newSize,0.0);
      a_h.resize(newSize,0.0);
      a_xi.resize(newSize,0.0);
      a_pt.resize(newSize,0.0); /* scratch vector for sole use of macprod */
      a_gx.resize(newSize,0.0); /* scratch gradients             */
      a_gy.resize(newSize,0.0); /* used by maclinmin and macprod */
      Parms.resize(newSize,0.0); 
  }
 
 /** optimize an object function
 * @param fn object function
 */
  bool optimize(ObjectFuncType* fn) {
    setTarget(fn);
    return optimize();
  }
  
  /** optimize TargetFunc */
  bool optimize(){
    bool successM = macoptII();
    //not sure if I need to make sure the parameters in the WF are correct.
    for(int i=0;i<NumParams;i++) {
      TargetFunc->Params(i)=Parms[i];
    }
    TargetFunc->Report();
    return successM;
  }
  
  ///write to a ostream
  bool get(std::ostream& ) const;
  
  ///read from istream
  bool put(std::istream& );
  
  ///read from an xmlNode 
  bool put(xmlNodePtr cur)
  {
  string steepestD("no");
  ParameterSet p;
  p.add(a_itmax,"max_steps","none"); p.add(a_itmax,"maxSteps","none");
  p.add(a_tol,"tolerance","none");
  p.add(GradTol,"tolerance_g","none"); p.add(GradTol,"toleranceG","none");
  p.add(GammaTol,"tolerance_cg","none"); p.add(GammaTol,"toleranceCG","none");
  p.add(Displacement,"epsilon","none");
  p.add(xycleanup,"xypolish","none");
  p.add(steepestD,"SD","none");
  p.add(a_rich,"rich","none");
  p.add(a_verbose,"verbose","none");
  p.add( a_lastx_default,"stepsize","none"); p.add( a_lastx_default,"stepSize","none");
  p.add( a_linmin_maxits,"max_linemin","none"); p.add( a_linmin_maxits,"maxLinemin","none");
  p.put(cur);
  
  SD=(steepestD=="yes");
  a_lastx=a_lastx_default;
  return true;
  }
  
  ///reset member data
  void reset(){}
  
  Return_t func(vector<Return_t> _p){
    for(int i=0;i<NumParams;++i) TargetFunc->Params(i) = _p[i] ;
    return TargetFunc->Cost();
  }
  void dfunc(vector<Return_t> RT, vector<Return_t> &FG)
  {
    //save initial parameters
    vector<Return_t> tmpPs(NumParams,0.0);
    for(int k=0; k<NumParams; k++) tmpPs[k] = TargetFunc->Params(k) ;
    Return_t dh=1.0/(2.0*Displacement);
    for(int i=0; i<NumParams ; i++) {
      for(int j=0; j<NumParams; j++) TargetFunc->Params(j)=RT[j];
      TargetFunc->Params(i) = RT[i]+ Displacement;
      Return_t CostPlus = TargetFunc->Cost(); 
      TargetFunc->Params(i) = RT[i]- Displacement;
      Return_t CostMinus = TargetFunc->Cost(); 
      FG[i]= (CostPlus-CostMinus)*dh;
    }
    //     for(int i=0;i<NumParams;i++) cout<<FG[i]<<"  ";
    //     cout<<NumParams<<endl;
    //leave initial parameters inside the wave function/cost function
    for(int k=0; k<TargetFunc->NumParams(); k++) TargetFunc->Params(k)=tmpPs[k];
  }

    bool macoptII( ) 
    { 
      Return_t gg(0.0) , gam(0.0) , dgg(0.0) ; 
      int end_if_small_grad = 1 - a_end_if_small_step ;
      Return_t step , tmpd ;
      dfunc( Parms, a_xi );
      macopt_restart ( 1 ) ;
      gg = CostTol*2.0;
      for (int a_its = 0 ; ((a_its < a_itmax)&(TargetFunc->IsValid)&(gg>CostTol)) ; a_its++ ) {
	gg = 0.0;
	for (int j = 0 ; j < NumParams ; j ++ )  gg += a_g[j]*a_g[j];
	a_gtyp = sqrt ( gg / (1.0*NumParams) ) ;
	if ( a_verbose > 0 ) 
	  printf ( "mac_it %d of %d : gg = %6.3g tol = %6.3g: ", a_its , a_itmax , gg , a_tol ) ;
	if ( ( end_if_small_grad && ( gg <= a_tol ) ) 
	  || ( gg <= a_grad_tol_tiny ) ) {
	  //      macopt_free ( a ) ;
	if ( a_verbose > 0 ) printf ("\n");
	return true;
	}
       
       step = maclinminII ( Parms ) ; 
       
       if ( a_restart == 0 ) {
	 if ( a_verbose > 1 ) printf (" (step %9.5g)",step);
	 if ( a_verbose > 0 ) printf ("\n");
	 if ( ( a_end_if_small_step  && ( step <= a_tol ) ) 
	   || ( step <= a_step_tol_tiny ) ) {
	   //	macopt_free ( a ) ;
	 return true;
	 }
       }
       
       /* if we are feeling rich, evaluate the gradient at the new
       `minimum'. alternatively, linmin has already estimated this
       gradient by linear combination of the last two evaluations and
       left it in xi */
    if ( a_rich || a_restart || SD) { 
      dfunc( Parms , a_xi  ) ; 
    }
    if ( a_restart ) {
      if(a_verbose > 0)
	fprintf(stderr,"Restarting macopt\n" ) ; 
      macopt_restart ( 0 ) ;
      /* this is not quite right
      should distinguish whether there was an overrun indicating that the 
      value of lastx needs to be bigger / smaller; 
      in which case resetting lastx to default value may be a bad idea, 
					giving an endless loop of resets 
					*/
    } 
    else if (SD)
    {
      fprintf(stderr,"Using Steepest Descent:\n" ) ;
      a_restart = 1 ;
      macopt_restart(1);
      gg = 0.0;
      for (int j = 0 ; j < NumParams ; j ++ ) 
	gg += a_g[j]*a_g[j]; 
    }
    else
    {
      dgg=0.0;
      for (int j = 0 ; j < NumParams ; j ++ ) {
	dgg += ( a_xi[j] + a_g[j] ) * a_xi[j] ;
      }
      gam = dgg / gg ;
      tmpd = 0.0;
      for ( int j = 0 ; j < NumParams ; j ++ ) {
	a_g[j] = -a_xi[j];                /* g stores (-) the most recent gradient */
	a_xi[j] = a_h[j] = a_g[j] + gam * a_h[j] ;
	/* h stores xi, the current line direction */
	/* check that the inner product of gradient and line search is < 0 */
	tmpd -= a_xi[j] * a_g[j] ; 
      }
      gg = 0.0;
      for (int j = 0 ; j < NumParams ; j ++ ) 
	gg += a_g[j]*a_g[j]; 
      if ( tmpd > 0.0  || a_verbose > 2 ) {
	fprintf(stderr,"new line search has inner prod %9.4g\n", tmpd ) ; 
      }
      if ( tmpd > 0.0 ) { 
	if ( a_rich == 0 ) {
	  fprintf (stderr, "Setting rich to 1; " ) ; 
	  a_rich = 1 ; 
	}
	a_restart = 2 ; /* signifies that g[j] = -xi[j] is already done */
	if(a_verbose > 0)
	  fprintf(stderr,"Restarting macopt (2)\n" ) ; 
	macopt_restart ( 0 ) ;
      }
    }
      }
      if((a_verbose> 0 )&( TargetFunc->IsValid))
	fprintf(stderr,"Reached iteration limit in macopt; continuing.\n"); 
      else
	fprintf(stderr,"TargetFunc->IsValid = FALSE\n");
      //  macopt_free ( a ) ;	
      a_lastx *= 0.5;
      return (( TargetFunc->IsValid)&(gg <= a_tol));
    } /* NB this leaves the best value of p in the p vector, but
    the function has not been evaluated there if rich=0     */
 
  Return_t maclinminII(vector<Return_t> &p)
  {
    Return_t y(0.0) ;
    Return_t  t(0.0) , m(0.0) ;
    Return_t step(0.0) , tmpd(0.0) ;  
    
    /* at x=0, the gradient (uphill) satisfies s < 0 */
    if ( a_verbose > 2 ) { 
      fprintf (stderr, "inner product at 0 = %9.4g\n" ,
		tmpd = macprodII ( Parms , a_gy , 0.0 ) ) ; 
		if ( tmpd > 0.0 ) { 
		  a_restart = 1 ; 
		  return 0.0 ; 
		}
    }
    int its(0);
    Return_t x = a_lastx / a_gtyp ;
    Return_t s = macprodII ( Parms , a_gx , x ) ; 
    
    if ( s < 0 & TargetFunc->IsValid)  {  /* we need to go further */
      do {
	y = x * a_linmin_g1 ;
	t = macprodII ( Parms , a_gy , y  ) ; 
	if ( a_verbose > 1 ) 
	  printf ("s < 0: s = %6.3g: t = %6.3g; x = %6.3g y = %6.3g\n",s, t , x , y );
	if (( t >= 0.0 ) | ( !TargetFunc->IsValid)){
	  printf ("CHANGED sign: s = %6.3g: t = %6.3g; x = %6.3g y = %6.3g\n",s, t , x , y );
	  break ;
	}
	if (t==t)
	{
	  x = y ; s = t ; a_gunused = a_gx ; a_gx = a_gy ; a_gy = a_gunused ; 
	}
	else
	{
	  y=x; t=s;
	  its = a_linmin_maxits;
	}
	its++ ;
	/* replaces: for ( i = 1 ; i <= n ; i ++ ) gx[i] = gy[i] ; */
      }
      while ( (its <= a_linmin_maxits)&(TargetFunc->IsValid) ) ;
      
      
    } else if ( s > 0 & TargetFunc->IsValid) { /* need to step back inside interval */
      do {
	y = x * a_linmin_g3 ;
	t = macprodII ( Parms , a_gy , y ) ; 
	if ( a_verbose > 1 ) 
	  printf ("s > 0: s = %6.3g: t = %6.3g; x = %6.3g y = %6.3g\n",s, t , x , y );
	if (( t <= 0.0 ) | ( !TargetFunc->IsValid)){
	  printf ("CHANGED sign: s = %6.3g: t = %6.3g; x = %6.3g y = %6.3g\n",s, t , x , y );
	  break ;
	}
	if (t==t)
	{
	  x = y ; s = t ; a_gunused = a_gx ; a_gx = a_gy ; a_gy = a_gunused ; 
	}
	else
	{
	  y=x; t=s;
	  its = a_linmin_maxits;
	}
	its ++ ;
      } while ( (its <= a_linmin_maxits)&(TargetFunc->IsValid) ) ;
    } else { /* hole in one s = 0.0 */
      t = 1.0 ; y = x;
    }
    
    Return_t u(0.0);
    int xyit(0);
    do {
//       Return_t s2 = std::fabs(s);
//       Return_t t2 = std::fabs(t);
      m = 0.5*( x+y ) ;
      //       s2 /= m ; t2 /= m ;
      //       m =  s2 * y + t2 * x ;
      printf ("Polishing line minimization m= %6.3g", m);
      u=macprodII ( Parms , a_gunused , m ) ;
      if (u==u){
	if ( u*s >= 0.0  & TargetFunc->IsValid )
	{
	  s=u; x=m; a_gx = a_gunused; 
	}
	else if ( u*t >= 0.0  & TargetFunc->IsValid )
	{
	  t=u; y=m; a_gy = a_gunused ;
	}
// 	Return_t gd=0.0;
// 	for (int j = 0 ; j < NumParams ; j ++ ) 
// 	  gd += a_gunused[j]*a_gunused[j];
	
// 	if ( (gd<GradTol) & (TargetFunc->IsValid) )
// 	{
// 	  t=s=u; y=x=m; a_gy =a_gx = a_gunused ;
// 	}
	if ( a_verbose > 1 ) 
	{
	  printf ("s = %6.3g: t = %6.3g; x = %6.3g y = %6.3g\n",s, t , x , y );
	  printf ("u = %6.3g: m = %6.3g; \n",u, m );
	}
      }
      else if (std::fabs(t)>std::fabs(s))
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
    } while ( (TargetFunc->IsValid) & (std::fabs(s-t)>GradTol) & (xyit<xycleanup) ) ;
    
    if ( its > a_linmin_maxits )  {
      fprintf (stderr, "Warning! maclinmin overran" );
      fprintf (stderr, "- inner product at 0 = %9.4g\n" ,
		tmpd = macprodII ( Parms , a_gy , 0.0 ) ) ; 
		if ( tmpd > 0 && a_rich == 0 ) {
		  fprintf (stderr, "setting rich to 1\n" ) ;       a_rich = 1 ; 
		}
		if ( tmpd > 0 ) {
		  fprintf (stderr, "setting a_restart to 1\n" ) ;       a_restart = 1 ;
		} 
    }
    if (TargetFunc->IsValid){
      /*  Linear interpolate between the last two. 
      This assumes that x and y do bracket. */
    s = std::fabs(s);
    t = std::fabs(t);
    m = ( s+t ) ;
    s /= m ; t /= m ;
    
    m =  s * y + t * x ; 
    /* evaluate the step length, not that it necessarily means anything */
//     cout<<"Setting new parameters"<<endl;
    step = 0.0;
    for (int i = 0 ; i < NumParams ; i ++ ) {
      tmpd = m * a_xi[i] ;
      Parms[i] += tmpd ; /* this is the point where the parameter vector steps */
      step += fabs ( tmpd ) ; 
      a_xi[i] = s * a_gy[i] + t * a_gx[i] ;
      /* send back the estimated gradient in xi (NB not like linmin) */
//       cout<<i<<"  "<<Parms[i]<<endl;
    }
    a_lastx = m * a_linmin_g2 *  a_gtyp ;
    for(int j=0; j<NumParams; j++) TargetFunc->Params(j)=Parms[j];
    TargetFunc->Cost();
    TargetFunc->Report();
    }
    else
    {
      cout<<" TargetFunc->IsValid=False" <<endl;
    }
    return ( step / (1.0*NumParams ) ) ; 
  }
  
  
  Return_t macprodII (vector<Return_t> p , vector<Return_t> &gy , Return_t y  )
  { 
    Return_t s = 0.0 ;
    for (int i = 0 ; i < NumParams ; i ++ ) 
      a_pt[i] = p[i] + y * a_xi[i] ;
    dfunc( a_pt , gy ) ;
    for (int i = 0 ; i < NumParams ; i ++ ) 
      s += gy[i] * a_xi[i] ;
    return s ;
  }
  
  void macopt_restart ( int start)
  {  
    if ( start == 0 ) a_lastx = a_lastx_default ; 
    /* it is assumed that dfunc( p , xi  ) ;  has happened */
  for (int j = 0 ; j < NumParams ; j ++ ) {
    if ( a_restart != 2 ) a_g[j] = -a_xi[j] ;
    a_xi[j] = a_h[j] = a_g[j] ;
  }
  a_restart = 0 ; 
  }
  
    };
    
    
    
    // template<class T>
    // MacOptimization<T>::MacOptimization(ObjectFuncType* atarget):
    //   TargetFunc(atarget), RestartCG(true),
	 //   NumSteps(100),Displacement(1e-6),
	 //   CostTol(1.e-6),GradTol(1.e-6),GammaTol(1.e-7) // GammaTol was set to 1e-4 originally
	 // {
	   //   curCost=prevCost=1.0e13;
	   // }
	   
	   // template<class T>
	   // void MacOptimization<T>::setTarget(ObjectFuncType* fn) {
    //   TargetFunc=fn;
    //   NumParams=TargetFunc->NumParams();
    // //   Y.resize(NumParams);
    // //   gY.resize(NumParams,0);
    // //   cgY.resize(NumParams,0);
    // //   for(int i=0;i<NumParams; i++) {
    // //     Y[i]=TargetFunc->Params(i);
    // //   }
    // 
    //   cout<<"setting Target "<<NumParams<<endl;
    //    setParams(NumParams,  2,  CostTol, NumSteps, 1);
    //   
    //     for(int i=0;i<NumParams;++i) {
    //     Parms[i]= TargetFunc->Params(i);
    //   }
    // //   evaluateGradients( a_g);
    // //   a_h=a_xi=a_g;
    // }
    
    // template<class T>
    // bool MacOptimization<T>::optimize() { 
 //   Parms.resize(NumParams);
 //   for(int i=0;i<NumParams;++i) {
 //     Parms[i]=TargetFunc->Params(i);
 // //     cout<<Parms[i]<<endl;
 //   }
 //   macoptII(Parms,NumParams);
 //   for(int i=0;i<NumParams;++i) {
 //     TargetFunc->Params(i)=Parms[i];
 // //     cout<<Parms[i]<<endl;
 //   }
 //   cout<<"Cost function :"<<TargetFunc->Cost()<<endl;
 //   return true;
 // }
 
 // template<class T>
 // void 
 // MacOptimization<T>::evaluateGradients(std::vector<Return_t> &grad) {
 // 
 //   //use targetFunc evaluateGradients if it does it better
 // //   if(TargetFunc->evaluateGradients(grad)) return;
 // 
 //   //do the finite difference method
 //   Return_t dh=1.0/(2.0*Displacement);
 //   for(int i=0; i<NumParams ; i++) {
 //     for(int j=0; j<TargetFunc->NumParams(); j++) TargetFunc->Params(j)=Parms[j];
 //     TargetFunc->Params(i) = Parms[i]+ Displacement;
 //     Return_t CostPlus = TargetFunc->Cost();
 //     
 //     TargetFunc->Params(i) = Parms[i]- Displacement;
 //     Return_t CostMinus = TargetFunc->Cost();
 // //     cout<<CostMinus<<endl<<endl;
 //     grad[i]=(CostMinus-CostPlus)*dh;
 // //     cout<<"Gradient["<<i<<"] = "<<(CostMinus-CostPlus)*dh<<endl;
 //   }
 // 
 // }
 // template<class T>
 // void 
 // MacOptimization<T>::dfunc(vector<Return_t> _p, vector<Return_t> _xi)
 // {
   //   for(int i=0;i<NumParams;i++) TargetFunc->Params(i) = _p[i] ;
   //   evaluateGradients( _xi);
   // }
   
   // template<class T>
   // typename MacOptimization<T>::Return_t 
   // MacOptimization<T>::Func(Return_t dl) {
 //   for(int i=0; i<NumParams; i++) TargetFunc->Params(i)=Y[i]+dl*cgY[i];
 //   return TargetFunc->Cost();
 // }
 
 
 // template<class T>
 // bool MacOptimization<T>::get(std::ostream& os) const {
   //   return true;
   // }
   // 
   // template<class T>
   // bool MacOptimization<T>::put(std::istream& is) {
 //   return true;
 // }
 // 
 // template<class T>
 // void MacOptimization<T>::reset() {
 // }
 
 // template<class T>
 // bool MacOptimization<T>::put(xmlNodePtr cur) {
 //   ParameterSet p;
 //   p.add(NumSteps,"max_steps","none"); p.add(NumSteps,"maxSteps","none");
 //   p.add(CostTol,"tolerance","none");
 //   p.add(GradTol,"tolerance_g","none"); p.add(GradTol,"toleranceG","none");
 //   p.add(GammaTol,"tolerance_cg","none"); p.add(GammaTol,"toleranceCG","none");
 //   p.add(Displacement,"epsilon","none");
 //   p.add( LambdaMax,"stepsize","none"); p.add(this->LambdaMax,"stepSize","none");
 // //   p.add( ITMAX,"max_linemin","none"); p.add(this->ITMAX,"maxLinemin","none");
 //   p.put(cur);
 //   return true;
 // }
 
 // template<class T>
 // void MacOptimization<T>::setParams(int n, 
 // 	       int _verbose,
 // 	       Return_t _tolerance,
 // 	       int _itmax,
 // 	       int _rich) 
 //  {
   //    cout<<"Setting parameters for optimization"<<endl;
   //   NumParams=(n);
   //   a_tol=(_tolerance);
   //   a_itmax=(_itmax);
   //   a_rich=(_rich);
   //   a_verbose=(_verbose);
   // 
   //   a_g.resize(n+1,0.0);
   //   a_h.resize(n+1,0.0);
   //   a_xi.resize(n+1,0.0);
   //   a_pt.resize(n+1,0.0); /* scratch vector for sole use of macprod */
//   a_gx.resize(n+1,0.0); /* scratch gradients             */
//   a_gy.resize(n+1,0.0); /* used by maclinmin and macprod */
//   /* if verbose = 1 then there is one report for each
//   line minimization.
//   if verbose = 2 then there is an additional 
//     report for
//     each step of the line minimization.
//     if verbose = 3 then extra debugging 
//      routines kick in.
//      */
// 
//   a_end_if_small_step = 1 ; /* Change this to 0/1 if you prefer */
//   a_stepmax = 0.5 ; 
// 
//   a_grad_tol_tiny = 1e-16 ; /* Probably not worth fiddling with */
//   a_step_tol_tiny = 0.0 ;   /* Probably not worth fiddling with */
//   a_linmin_maxits = 20 ;    /* Probably not worth fiddling with */
//   a_lastx = 0.01 ;          /* only has a transient effect      */
//   a_lastx_default = 0.01 ;  /* -- defines typical distance in parameter
// 				space at which the line minimum is expected;
// 				both these should be set. the default is 
// 				consulted if something goes badly wrong and 
// 				a reset is demanded. */
// 
// /* don't fiddle with the following, unless you really mean it */
//   a_linmin_g1 = 2.0 ; 
//   a_linmin_g2 = 1.25 ; 
//   a_linmin_g3 = 0.5 ; 
//   a_restart = 0 ; 
// }


// template<class T>
// void MacOptimization<T>::macoptII
//   (vector<Return_t> p,            /* starting vector                                */
//    int    n             /* number of dimensions                           */
//    )                   
// {
  //   int j ;
  //   Return_t gg , gam , dgg ;
  //   vector<Return_t> g(a_g) , h(a_h) , xi(a_xi) ;
  //   int end_if_small_grad = 1 - a_end_if_small_step ;
  //   Return_t step , tmpd ;
  // 
  //   /* A total of 7 Return_t * 1..n are used by this optimizer. 
  //      p           is provided when the optimizer is called 
  //      pt          is used by the line minimizer as the temporary vector. 
  //                     this could be cut out with minor rewriting, using p alone
  //      g, h and xi are used by the cg method as in NR - could one of these
  //                      be cut out?
  //      the line minimizer uses an extra gx and gy to evaluate two gradients. 
  //      */
// 
//    dfunc(p,xi );
// 
//   macopt_restart ( 1 ) ; 
// 
//   for ( a_its = 1 ; a_its <= a_itmax ; a_its ++ ) {
 // 
 //     for ( gg = 0.0 , j = 1 ; j <= n ; j ++ ) 
 //       gg += g[j]*g[j];          /* find the magnitude of the old gradient */
//     a_gtyp = sqrt ( gg / (Return_t)(n) ) ; 
// 
//     if ( a_verbose > 0 ) 
//       printf ( "mac_it %d of %d : gg = %6.3g tol = %6.3g: ", a_its , a_itmax , gg , a_tol ) ;
// 
//     if ( ( end_if_small_grad && ( gg <= a_tol ) ) 
// 	|| ( gg <= a_grad_tol_tiny ) ) {
 //       //      macopt_free ( a ) ;
 //       if ( a_verbose > 0 ) printf ("\n");
 //       return;
 //     }
 // 
 //     step = maclinminII ( p  ) ; 
 // 
 //     if ( a_restart == 0 ) {
 //       if ( a_verbose > 1 ) printf (" (step %9.5g)",step);
 //       if ( a_verbose > 0 ) printf ("\n");
 //       if ( ( a_end_if_small_step  && ( step <= a_tol ) ) 
 // 	  || ( step <= a_step_tol_tiny ) ) {
 // 	//	macopt_free ( a ) ;
 // 	return;
 //       }
 //     }
 // 
 //     /* if we are feeling rich, evaluate the gradient at the new
 //        `minimum'. alternatively, linmin has already estimated this
 //        gradient by linear combination of the last two evaluations and
 //        left it in xi */
//     if ( a_rich || a_restart ) { 
 //       dfunc( p , xi  ) ; 
 //     }
 //     if ( a_restart ) {
 //       if(a_verbose > 0)
 // 	fprintf(stderr,"Restarting macopt\n" ) ; 
 //       macopt_restart ( 0 ) ;
 // /* this is not quite right
 //    should distinguish whether there was an overrun indicating that the 
 //    value of lastx needs to be bigger / smaller; 
 //    in which case resetting lastx to default value may be a bad idea, 
 //    giving an endless loop of resets 
 // */
//     } else {
  //       dgg=0.0;
  //       for ( j = 1 ; j <= n ; j ++ ) {
 // 	dgg += ( xi[j] + g[j] ) * xi[j] ;
 //       }
 //       gam = dgg / gg ;
 //       for ( tmpd = 0.0 , j = 1 ; j <= n ; j ++ ) {
 // 	g[j] = -xi[j];                /* g stores (-) the most recent gradient */
// 	xi[j] = h[j] = g[j] + gam * h[j] ;
// 	/* h stores xi, the current line direction */
// 	/* check that the inner product of gradient and line search is < 0 */
// 	tmpd -= xi[j] * g[j] ; 
//       }
//       if ( tmpd > 0.0  || a_verbose > 2 ) {
 // 	fprintf(stderr,"new line search has inner prod %9.4g\n", tmpd ) ; 
 //       }
 //       if ( tmpd > 0.0 ) { 
 // 	if ( a_rich == 0 ) {
 // 	  fprintf (stderr, "Setting rich to 1; " ) ; 
 // 	  a_rich = 1 ; 
 // 	}
 // 	a_restart = 2 ; /* signifies that g[j] = -xi[j] is already done */
// 	if(a_verbose > 0)
// 	  fprintf(stderr,"Restarting macopt (2)\n" ) ; 
// 	macopt_restart ( 0 ) ;
//       }
//     }
//   }
//   if(a_verbose> 0)
//     fprintf(stderr,"Reached iteration limit in macopt; continuing.\n"); 
//   //  macopt_free ( a ) ;	
//   return;
// } /* NB this leaves the best value of p in the p vector, but
//      the function has not been evaluated there if rich=0     */

/* maclinmin.
Method: 
evaluate gradient at a sequence of points and calculate the inner 
product with the line search direction. Continue until a 
bracketing is achieved ( i.e a change in sign ). */
// template<class T>
// Return_t MacOptimization<T>::maclinminII (vector<Return_t> p )
// {
  //   int n = NumParams ; 
  // 
  //   Return_t x , y ;
  //   Return_t s , t , m ;
  //   int    its = 1 , i ;
  //   Return_t step , tmpd ; 
  //   vector<Return_t> gx(a_gx) ,  gy(a_gy) ;
  // 
  //   /* at x=0, the gradient (uphill) satisfies s < 0 */
//   if ( a_verbose > 2 ) { /* check this is true: (no need to do this really
// 			   as it is already checked at the end of the main
// 			   loop of macopt) */
// /*
// #define TESTS 5
//     x = a_lastx / a_gtyp ;
//     fprintf (stderr, "inner product at:\n" ) ; 
//     for ( i = -TESTS ; i <= TESTS ; i += 2 ) {
 //       step = x * 2.0 * (Return_t) i / (Return_t) TESTS ; 
 //       fprintf (stderr, "%9.5g %9.5g\n" , step ,
 // 	       tmpd = macprodII ( p , gy , step ) ) ; 
 //     }
 // */
//     fprintf (stderr, "inner product at 0 = %9.4g\n" ,
// 	     tmpd = macprodII ( p , gy , 0.0 ) ) ; 
//     if ( tmpd > 0.0 ) { 
 //       a_restart = 1 ; 
 //       return 0.0 ; 
 //     }
 //   }
 // 
 //   x = a_lastx / a_gtyp ;
 //   s = macprodII ( p , gx , x ) ; 
 //   
 //   if ( s < 0 )  {  /* we need to go further */
//     do {
  //       y = x * a_linmin_g1 ;
  //       t = macprodII ( p , gy , y  ) ; 
  //       if ( a_verbose > 1 ) 
  // 	printf ("s = %6.3g: t = %6.3g; x = %6.3g y = %6.3g\n",s, t , x , y );
  //       if ( t >= 0.0 ) break ;
  //       x = y ; s = t ; a_gunused = gx ; gx = gy ; gy = a_gunused ; 
  //       its++ ;
  // /* replaces: for ( i = 1 ; i <= n ; i ++ ) gx[i] = gy[i] ; */
//     }
//     while ( its <= a_linmin_maxits ) ;
//   } else if ( s > 0 ) { /* need to step back inside interval */
//     do {
  //       y = x * a_linmin_g3 ;
  //       t = macprodII ( p , gy , y ) ; 
  //       if ( a_verbose > 1 ) 
  // 	printf ("s = %6.3g: t = %6.3g; x = %6.3g y = %6.3g\n",s, t , x , y );
  //       if ( t <= 0.0 ) break ;
  //       x = y ; s = t ; a_gunused = gx ; gx = gy ; gy = a_gunused ; 
  //       its ++ ;
  //     } while ( its <= a_linmin_maxits ) ;
  //   } else { /* hole in one s = 0.0 */
//     t = 1.0 ; y = x;
//   }
// 
//   if ( its > a_linmin_maxits )  {
 //     fprintf (stderr, "Warning! maclinmin overran" );
 // /* this can happen where the function goes \_ and doesn't buck up
 //  again; it also happens if the initial `gradient' does not satisfy
 //  gradient.`gradient' > 0, so that there is no minimum in the supposed
 //  downhill direction.  I don't know if this actually happens... If it
 //  does then I guess a_rich should be 1.
 // 
 //  If the overrun is because too big a step was taken then
 //  the interpolation should be made between zero and the most 
 //  recent measurement. 
 // 
 //  If the overrun is because too small a step was taken then 
 //  the best place to go is the most distant point. 
 //  I will assume that this doesn't happen for the moment.
 // 
 //  Also need to check up what happens to t and s in the case of overrun.
 //  And gx and gy. 
 // 
 //  Maybe sort this out when writing a macopt that makes use of the gradient 
 //  at zero? 
 // */
//     fprintf (stderr, "- inner product at 0 = %9.4g\n" ,
// 	     tmpd = macprodII ( p , gy , 0.0 ) ) ; 
//     if ( tmpd > 0 && a_rich == 0 ) {
 //       fprintf (stderr, "setting rich to 1\n" ) ;       a_rich = 1 ; 
 //     }
 //     if ( tmpd > 0 ) a_restart = 1 ; 
 //   }
 // 
 //  /*  Linear interpolate between the last two. 
 //      This assumes that x and y do bracket. */
//   if ( s < 0.0 ) s = - s ;
//   if ( t < 0.0 ) t = - t ;
//   m = ( s + t ) ;
//   s /= m ; t /= m ;
//   
//   m =  s * y + t * x ; 
//   /* evaluate the step length, not that it necessarily means anything */
//   for ( step = 0.0 , i = 1 ; i <= n ; i ++ ) {
 //     tmpd = m * a_xi[i] ;
 //     p[i] += tmpd ; /* this is the point where the parameter vector steps */
//     step += fabs ( tmpd ) ; 
//     a_xi[i] = s * gy[i] + t * gx[i] ;
// /* send back the estimated gradient in xi (NB not like linmin) */
//   }
//   a_lastx = m * a_linmin_g2 *  a_gtyp ;
//   
//   return ( step / (Return_t) ( n ) ) ; 
// }

// template<class T>
// Return_t MacOptimization<T>::macprodII 
// ( 
//  vector<Return_t> p , vector<Return_t> gy , Return_t y  
// 
// ) {
 //   vector<Return_t> pt = a_pt ; 
 //   vector<Return_t> xi = a_xi ; 
 //   /* finds pt = p + y xi and gets gy there, 
 // 				       returning gy . xi */
//   int n = NumParams ; 
// 
//   int i;
//   Return_t s = 0.0 ;
// 
//   for ( i = 1 ; i <= n ; i ++ ) 
//     pt[i] = p[i] + y * xi[i] ;
//   
//   dfunc( pt , gy ) ;
// 
//   for ( i = 1 ; i <= n ; i ++ ) 
//     s += gy[i] * xi[i] ;
// 
//   return s ;
// }

// template<class T>
// void MacOptimization<T>::macopt_restart ( int start ) 
// /* if start == 1 then this is the start of a fresh macopt, not a restart */
// {
  //   int j , n=NumParams ; 
  //   vector<Return_t> g, h, xi ;
  //   g = a_g ;  h = a_h ;  xi = a_xi ; 
  // 
  //   if ( start == 0 ) a_lastx = a_lastx_default ; 
  //   /* it is assumed that dfunc( p , xi  ) ;  has happened */
//   for ( j = 1 ; j <= n ; j ++ ) {
 //     if ( a_restart != 2 ) g[j] = -xi[j] ;
 //     xi[j] = h[j] = g[j] ;
 //   }
 //   a_restart = 0 ; 
 // }
 
 
 // template<class T>
 // void MacOptimization<T>::maccheckgrad 
 // /* Examines objective function and d_objective function to see if 
 //    they agree for a step of size epsilon */
//   (vector<Return_t> p,
//    int    n,
//    Return_t epsilon,
//    int    stopat          /* stop at this component. If 0, do the lot. */
// )
// {
  //   int j;
  //   Return_t f1;
  //   vector<Return_t> g, h;
  //   Return_t tmpp ; 
  // 
  //   f1 = func(p);
  //   //  printf("f1 = %f\n", f1);
  // 
  //   vector<Return_t> h = new vector<Return_t>(n+1,0.0);
  //   vector<Return_t> g = new vector<Return_t>(n+1,0.0);
  //  
  //   dfunc(p, g);
  // 
  //   if ( stopat <= 0 || stopat > n ) stopat = n ; 
  // 
  //   printf("Testing gradient evaluation\n");
  //   printf("      analytic     1st_diffs    difference\n");
  //   for ( j = 1 ; j <= stopat ; j ++ ) {
 //     tmpp = p[j] ; 
 //     p[j] += epsilon ;
 //     h[j] = func(p) - f1 ;
 //     //    printf("h = %f\n", h[j]);
 //     p[j] =  tmpp ;
 // 
 //     printf("%2d %9.5g %9.5g %9.5g\n" , j , g[j] , h[j]/epsilon , g[j] - h[j]/epsilon );
 //     fflush(stdout) ; 
 //   }
 //   free(h);
 //   free(g);
 //   //  free_dvector(h,1,n);
 //   //  free_dvector(g,1,n);
 //   printf("      --------     ---------\n");
 // }
 

#endif

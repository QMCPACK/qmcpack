//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Paul R. C. Kent, kentpr@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Paul R. C. Kent, kentpr@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
    
    
//           http://pathintegrals.info                     //
/////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////
//                                                                         //
// Adaptive Gauss-Kronrod integration                                      //
//                                                                         //
// This C++ version was written by Burkhard Militzer  Livermore 02-20-02   //
// Based on the GNU scientific library                                     //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

#ifndef _GKINTEGRATION_
#define _GKINTEGRATION_

#include <iterator>
#include <list>
#include "Standard.h"

class GK15 {
 public:
  static const int n=8;
  static const double xgk[n];
  static const double wg[n/2];
  static const double wgk[n];
};

class GK21 {
 public:
  static const int n=11;
  static const double xgk[n];
  static const double wg[n/2];
  static const double wgk[n];
};

class GK31 {
 public:
  static const int n=16;
  static const double xgk[n];
  static const double wg[n/2];
  static const double wgk[n];
};

class GK41 {
 public:
  static const int n=21;
  static const double xgk[n];
  static const double wg[(n+1)/2];
  static const double wgk[n];
};

class GK51 {
 public:
  static const int n=26;
  static const double xgk[n];
  static const double wg[n/2];
  static const double wgk[n];
};

class GK61 {
 public:
  static const int n=31;
  static const double xgk[n];
  static const double wg[n/2];
  static const double wgk[n];
};

////////////////////////////////////////////////////////////////////////////////////////

template <class F, class GKRule=GK31> 
class GKIntegration {
 private:

  class IntervalResult {
  public:
    IntervalResult(const double a_, const double b_,
		   const double delta_):a(a_),b(b_),result(0.0),err(0.0),delta(delta_){};
    double a,b,result,err, delta;

    double ErrorL() const {
      return (delta ? err/delta : err);
    }

    friend std::ostream& operator<<(std::ostream &os, const IntervalResult & ir) {
      os << "[a= " << ir.a
	 << " b= " << ir.b
	 << " result= " << ir.result
	 << " error/L= " << ir.err/(ir.b-ir.a)
	 << " error= " << ir.err
	 << " ]";
      return os;
    }
  };

  std::list <IntervalResult> ir;
  F & f; // could be not const in case where calling f() actually modifies the object
  bool relativeErrors;

 public:
  GKIntegration(F & f_):f(f_),relativeErrors(false) {}
  
  void SetAbsoluteErrorMode() {
    relativeErrors=false;
  }
  void SetRelativeErrorMode() {
    relativeErrors=true;
  }

 private:
  // funnel all calls through this function and branch to specfic n knot rule
  void GK(IntervalResult & r) {
    GKGeneral(GKRule::n,GKRule::xgk,GKRule::wg,GKRule::wgk,r);
  }

  // handle all n knot rule with the passed in positions and weights
  void GKGeneral(const int n, 
		 const double xgk[], const double wg[], const double wgk[],
		 IntervalResult & r) {

    const double center     = 0.5 * (r.a + r.b);
    const double halfLength = 0.5 * r.delta;
    const double fCenter = f(center);

    double resultGauss   = 0;
    double resultKronrod = fCenter * wgk[n - 1];

    if (n % 2 == 0) {
      resultGauss = fCenter * wg[n / 2 - 1];
    }

    for (int j=0; j<(n-1)/2; j++) {
      const int jtw = j*2+1;	// j=1,2,3 jtw=2,4,6
      const double xx    = halfLength * xgk[jtw];
      const double fval1 = f(center - xx);
      const double fval2 = f(center + xx);
      const double fsum = fval1 + fval2;
      resultGauss   += wg[j]    * fsum;
      resultKronrod += wgk[jtw] * fsum;
    }
    
    for (int j=0; j<n/2; j++) {
      int jtwm1 = j*2;
      const double xx = halfLength * xgk[jtwm1];
      const double fval1 = f(center - xx);
      const double fval2 = f(center + xx);
      resultKronrod += wgk[jtwm1] * (fval1 + fval2);
    };
    
    /* scale by the width of the integration region */
    resultGauss   *= halfLength;
    resultKronrod *= halfLength;
    
    r.result = resultKronrod;
    r.err = std::abs(resultKronrod - resultGauss);
    //r.err = pow(200.0 * std::abs(resultKronrod - resultGauss), 1.5);
    //    BMWrite(r);
  }

  void PrintList() {
    std::cout << "/------------------------------------------\\" << std::endl;
    int i=0;
    for(typename std::list <IntervalResult>::iterator p=ir.begin(); p!=ir.end(); p++) {
      BMWrite2(i,*p);
      i++;
    }
    std::cout << "\\------------------------------------------/" << std::endl;
  }

  //Print interval with maxium error per interval length
  //(not with maximum error - this is on top of the list)
  void PrintMax() {
    typename std::list <IntervalResult>::iterator rMax(ir.begin());

    for(typename std::list <IntervalResult>::iterator r=ir.begin()++; r!=ir.end(); r++) {
      if ((*r).ErrorL()>(*rMax).ErrorL())
	rMax=r;
    }

    BMWrite(*rMax);
  }

  void CheckList() {
    if (ir.size()==0) return;
    
    for(typename std::list <IntervalResult>::iterator p=ir.begin(); p!=ir.end(); p++) {
      typename std::list <IntervalResult>::iterator pn = p;
      pn++;
      if (pn!=ir.end())
	if ( ((*p).err) < ((*pn).err) ) {
	  PrintList();
	  BMWrite2(*pn,*p);
	  ::error("Ordering problem in list");
	}
    }
  }

  void CheckError(const double err) {
    double errorSum=0.0;

    if (ir.size()>0) {
      for(typename std::list <IntervalResult>::iterator p=ir.begin(); p!=ir.end(); p++) {
	errorSum += (*p).err;
      }
    }
    
    if (errorSum==0.0) {
      if (err!=0.0)
	error("CheckError",errorSum,err);
    } else {
      if (err/errorSum-1.0>1e-8 && std::abs(err-errorSum)>1e-14) 
	error("CheckError",errorSum,err,errorSum-err);
    }

    BMWrite4("PassedErrorCheck",errorSum,err,errorSum-err);
  }

  double RecomputeError() {
    double errorSum=0.0;

    if (ir.size()>0) {
      for(typename std::list <IntervalResult>::iterator p=ir.begin(); p!=ir.end(); p++) {
	errorSum += (*p).err;
      }
    }
    
    return errorSum;
  }

  void Insert(const IntervalResult & r) {
    //    std::cout << "Inserting.." << std::endl;
    //    PrintList();

    if (ir.empty()) {
      ir.push_back(r);
      return;
    }

    if (ir.back().err>=r.err) {
      ir.push_back(r);
      return;
    }

    if (r.err>=ir.front().err) {
      ir.push_front(r);
      return;
    }

    // size must be >=2
    typename std::list <IntervalResult>::iterator p = ir.end();
    p--;

    // p cannot become less the begin() because of check above
    while (r.err > (*p).err)
      p--;

    // go one down because insert put the element before p
    p++;
    ir.insert(p,r);
    //    CheckList();
  }

  double Integrate(const double a, const double b, 
		   const double absError, const bool absErrorFlag, 
		   const double relError, const bool relErrorFlag, 
		   const bool andFlag) {
    ir.clear();

    // #define PRINT_IT
#ifdef PRINT_IT
    std::cout << "Beginning integration" << std::endl;
#endif

    double errorUnresolved=0.0;
    const int iterationMax=30;
    // double lengthMin = (b-a)*pow(0.5,iterationMax);
    double lengthMin = ldexp(b-a,-iterationMax);
    
    IntervalResult r0(a,b,b-a);
    GK(r0);
    double result =r0.result;
    double err    =r0.err;
    double errLast=err;

    ir.push_back(r0);

    bool quitFlag;
    do {
      // PrintList();

      // Test if the interval with the biggest error has already been subdivided
      // the maximum number of times. If this is the case throw it out and print add
      // this contribution to the 'unresolved' errors to be printed at the end
      while (ir.size()>0) {
	IntervalResult & rTest (ir.front());
	double lengthTest = rTest.delta;
	if (lengthTest<lengthMin) {
	  warning("KC:Interval was divided too many times",iterationMax,
		  rTest.a,rTest.b,rTest.err,ir.size());
	  //	  warning("KC:current result=",result,"error=",err);
	  if (absErrorFlag) warning("Absolute accuracy = ",absError);
	  if (relErrorFlag) warning("Relative accuracy = ",relError,
				    "->absolute accuracy=",relError*std::abs(result));
	  // this means there is a problem with the integrand->you could exit here
	  //	  exit(1);
	  errorUnresolved += rTest.err;
	  //	  PrintList();
	  ir.pop_front();
	} else break;
      }
      // do you want to exit with a warning after the first unresolved sub-interval occured?
      if (ir.size()==0 || errorUnresolved>0.0) break;
      // or print as many as occur
      //      if (ir.size()==0) break;

      IntervalResult & r (ir.front());

      double center = 0.5*(r.a+r.b);
      IntervalResult r1(r.a,center,0.5*r.delta);
      IntervalResult r2(center,r.b,0.5*r.delta);

      GK(r1);
      GK(r2);

      // must not use r after popping
      result += r1.result+r2.result - r.result;
      err    += r1.err   +r2.err    - r.err;

#ifdef PRINT_IT
      std::cout.setf(std::ios::scientific);
      std::cout << "Refined [ " << r.a << " " << r.b 
	   << " ] err/L=" << (r1.err +r2.err)/(r.b-r.a)
      	   << " error=" << err << std::endl;
#endif

      // must remove old element first because new ones could move to top
      ir.pop_front();

      Insert(r1);
      Insert(r2);

      // In rare events, the error decreases over may (>10) orders of magnitude
      // during the refinement. Rounding errors from the beginning can prevent
      // err from becoming small enough. Recompute err after a substantial decrease.
      if (err<1e-6*errLast) {
	err     = RecomputeError();
	errLast = err;
      }
 
      //      CheckError(err);
      //       PrintList();
      
      const bool relOk = (err < relError*std::abs(result) || result==0.0);
      const bool absOk = (err < absError);

      if (absErrorFlag && relErrorFlag) {
	quitFlag = andFlag ? (relOk && absOk) : (relOk || absOk);
      } else {
	quitFlag = absErrorFlag ? absOk : relOk;
      }

    } while (!quitFlag);

#ifdef PRINT_IT
    PrintMax();
#endif

    if (errorUnresolved>0.0) {
      warning("KC:Unresolved error sum=",errorUnresolved,"for integration interval",a,b);
      warning("KC:--> Result=",result,"total error=",err,"rel. error=",
	      ((result!=0.0) ? err/std::abs(result) : 0.0));
      //      if (absErrorFlag) warning("Absolute accuracy = ",absError);
      //      if (relErrorFlag) warning("Relative accuracy = ",relError,
      //				"->absolute accuracy=",relError*std::abs(result));
    }

    //    CheckList();
#ifdef PRINT_IT
    std::cout << "End integration" << std::endl;
#endif
    double sum = 0.0;
    int numIntervals = 0;
    for(typename std::list<IntervalResult>::iterator p=ir.begin();p!=ir.end();p++) {
      sum += p->result;
      numIntervals++;
    }

//     if (numIntervals > 2000)
//       std::cerr << "Number of intervals = " << numIntervals << std::endl;

    double badSum = std::abs((result-sum) / sum);
    if ((badSum > 1.0e-7) && (std::abs(result-sum) > absError)) {
      std::cerr << "absError tolerance = " << absError << std::endl;
      std::cerr << "Percent error = " << badSum*100.0 << std::endl;
      std::cerr << "Number of intervals = " << numIntervals << std::endl;
      std::cerr << "result = " << result << " sum = " << sum << std::endl;
    }


    return result;
  }

 public:

  double Integrate(const double a, const double b, 
		   const double acc) {
    return Integrate(a,b,acc,!relativeErrors,acc,relativeErrors,false);
  }
  double Integrate(const double a, const double b, 
		   const double accAbs, const double accRel,
		   const bool andFlag) {
    return Integrate(a,b,accAbs,true,accRel,true,andFlag);
  }

};

#endif // _GKINTEGRATION_

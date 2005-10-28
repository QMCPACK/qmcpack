//////////////////////////////////////////////////////////////////
// (c) Copyright 1998-2002,2003- by Jeongnim Kim 
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
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
//   Department of Physics, Ohio State University
//   Ohio Supercomputer Center
//////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_CG_OPTIMIZATION_NRC_H
#define QMCPLUSPLUS_CG_OPTIMIZATION_NRC_H

#include "Optimize/OptimizeBase.h"
#include "Optimize/NRCOptimization.h"
#include "OhmmsData/ParameterSet.h"

template<class T>
class CGOptimization: public MinimizerBase<T>,
                      private NRCOptimization<T> {
  public:


  /** Function to be optimized.  */
  ObjectFuncType* TargetFunc;

  /** number of minimization steps */
  int NumSteps;
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

  /** constructor
   * @param atarget target function to optimize
   */
  CGOptimization(ObjectFuncType* atarget=0);

  /** destructor */
  ~CGOptimization() {}

  /** set the target function
   * @param fn target function to optimize
   *
   * Allocate the internal data.
   */
  void setTarget(ObjectFuncType* fn);

  /** optimize an object function
   * @param fn object function
   */
  bool optimize(ObjectFuncType* fn) {
    setTarget(fn);
    return optimize();
  }

  /** optimize TargetFunc */
  bool optimize();

  ///write to a ostream
  bool get(std::ostream& ) const;

  ///read from istream
  bool put(std::istream& );

  ///read from an xmlNode 
  bool put(xmlNodePtr cur);

  ///reset member data
  void reset();

 protected:

  bool RestartCG;
  int NumParams;
  int CurStep;

  Return_t gdotg, gdotg0, gdoth, gamma;
  Return_t curCost, prevCost;

  std::vector<Return_t> Y;
  std::vector<Return_t> gY, cgY, gY0;

  /** evaluate the value for y+dl*cg 
   *
   * Lineminimization uses this function to find the minimum along the CG direction
   */
  Return_t Func(Return_t dl);

  /** evaluate the gradients numerically 
   * @param grad container for the gradients
   */
  void evaluateGradients(std::vector<Return_t>& grad);
};

template<class T>
inline T 
dotProduct(const std::vector<T>& a, const std::vector<T>& b) {
  T res=0.0;
  for(int i=0; i<a.size(); i++) res += a[i]*b[i];
  return res;
}

template<class T>
CGOptimization<T>::CGOptimization(ObjectFuncType* atarget):
  TargetFunc(atarget), RestartCG(true),
  NumSteps(100),Displacement(1e-6),
  CostTol(1.e-6),GradTol(1.e-6),GammaTol(1.e-4)
{
  prevCost=10e6; //some large number
}

template<class T>
void CGOptimization<T>::setTarget(ObjectFuncType* fn) {
  TargetFunc=fn;
  NumParams=TargetFunc->NumParams();
  Y.resize(NumParams);
  gY.resize(NumParams,0);
  cgY.resize(NumParams,0);
  for(int i=0;i<NumParams; i++) {
    Y[i]=TargetFunc->Params(i);
  }
}

template<class T>
bool CGOptimization<T>::optimize() {

  CurStep=0;
  do {
    if(RestartCG) { //first time
      evaluateGradients(gY);
      gdotg = dotProduct(gY,gY);
      gdotg0 = gdotg;
      gdoth  = gdotg;
      gamma  = 0.0e0;
      cgY= gY;
      gY0=gY;
      RestartCG=false;
    }

    bool success = this->lineoptimization();
    success &= TargetFunc->IsValid;
    if(success) { //successful lineminimization
      curCost= Func(this->Lambda);
      for(int i=0; i<NumParams; i++) Y[i]+=this->Lambda*cgY[i];
    } else {
      if(msg_stream) {
        *msg_stream << "Stop CGOptimization due to the failure of line optimization" << std::endl;
        *msg_stream << "Total number of steps = " << CurStep << std::endl;
      }
      return false;
    }

    evaluateGradients(gY);
    Return_t fx= fabs(*(std::max_element(gY.begin(), gY.end())));
    gdotg0=gdotg;
    gdotg=dotProduct(gY,gY);
    //Do not check the component yet
    //gdotg=Dot(dY,dY,fx);
    //if(fx<GradMaxTol) {
    if(fx<GradTol) {
      if(msg_stream) *msg_stream << " CGOptimization  has reached gradient max|G| = " << fx << endl;
      return false;
    }
    //if(gdotg < GradTol) {
    //  *msg_stream << " CGOptimization::Converged gradients" << endl;
    //  return false; 
    //}

    gdoth = dotProduct(gY,gY0);
    gamma = (gdotg-gdoth)/gdotg0;
    gY0=gY; //save the current gradient 

    if(fabs(gamma) < GammaTol){
      if(msg_stream) 
      *msg_stream << " CGOptimization::Converged conjugate gradients" << endl;
      return false; 
    }

    if(gamma > 1.0e2) {
      if(msg_stream) 
      *msg_stream << " CGOptimization restart: " << gamma << " is too big." << endl;
      RestartCG = true;
    }

    Return_t dx=fabs((curCost-prevCost)/curCost);
    if(dx <= CostTol) {
      if(msg_stream) 
      *msg_stream << " CGOptimization::Converged cost with " << dx << endl;
      return false; //
    }

    prevCost=curCost;
    for(int i=0; i<NumParams; i++) cgY[i]=gY[i]+gamma*cgY[i];

    CurStep++;

    if(TargetFunc->IsValid) {
      TargetFunc->Report();
    } else {
      if(msg_stream)
        *msg_stream << " CGOptimization stopped due to invalid cost values " << endl;
      return false;
    }

  } while(CurStep<NumSteps);

  if(msg_stream) *msg_stream << " Failed to converged after " << CurStep << " steps." << std::endl;
  return false;
}

template<class T>
void 
CGOptimization<T>::evaluateGradients(std::vector<Return_t>& grad) {
  Return_t dh=1.0/(2.0*Displacement);
  for(int i=0; i<TargetFunc->NumParams() ; i++) {
    for(int j=0; j<TargetFunc->NumParams(); j++) TargetFunc->Params(j)=Y[j];
    TargetFunc->Params(i) = Y[i]+ Displacement;
    Return_t CostPlus = TargetFunc->Cost();
    TargetFunc->Params(i) = Y[i]- Displacement;
    Return_t CostMinus = TargetFunc->Cost();
    grad[i]=(CostMinus-CostPlus)*dh;
  }
}

template<class T>
typename CGOptimization<T>::Return_t 
CGOptimization<T>::Func(Return_t dl) {
  for(int i=0; i<NumParams; i++) TargetFunc->Params(i)=Y[i]+dl*cgY[i];
  return TargetFunc->Cost();
}


template<class T>
bool CGOptimization<T>::get(std::ostream& os) const {
  return true;
}

template<class T>
bool CGOptimization<T>::put(std::istream& is) {
  return true;
}

template<class T>
void CGOptimization<T>::reset() {
}

template<class T>
bool CGOptimization<T>::put(xmlNodePtr cur) {
  ParameterSet p;
  p.add(NumSteps,"min_steps","none");
  p.add(CostTol,"tolerance","none");
  p.add(GradTol,"tolerance_g","none");
  p.add(GammaTol,"tolerance_cg","none");
  p.add(this->LambdaMax,"stepsize","none");
  p.add(Displacement,"epsilon","none");
  p.add(this->ITMAX,"brent_steps","none");
  p.put(cur);
  return true;
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
#endif

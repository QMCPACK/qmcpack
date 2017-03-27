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
    
    


#ifndef QMCPLUSPLUS_STEEPESTDESCENT_OPTIMIZATION_H
#define QMCPLUSPLUS_STEEPESTDESCENT_OPTIMIZATION_H

#include "Optimize/OptimizeBase.h"
#include "OhmmsData/ParameterSet.h"

template<class T>
class DampedDynamics: public MinimizerBase<T>
{

public:

  typedef T Return_t;
  typedef typename MinimizerBase<T>::ObjectFuncType ObjectFuncType;
  using MinimizerBase<T>::msg_stream;

  /** Function to be optimized.  */
  ObjectFuncType* TargetFunc;

  /** number of minimization steps */
  int NumSteps;
  /** displacement to evaluate the gradients numerically */
  Return_t Displacement;
  /** fictitious time step for the damped dynamics, default=1 */
  Return_t Dt;
  /** friction for the damped dynamics, default= 0.5 */
  Return_t Friction;
  /** fictitious mass for the damped dynamics, default=10*/
  Return_t Mass;
  /** tolerance the converged value */
  Return_t CostTol;
  /** tolerance for the converged gradients ||g|| < GradTol */
  Return_t GradTol;
  /** tolerance for the converged gradient component max|g[i]| < GradMaxTol */
  Return_t GradMaxTol;

  /** constructor
   * @param atarget target function to optimize
   */
  DampedDynamics(ObjectFuncType* atarget=0);

  /** destructor */
  ~DampedDynamics() {}

  /** set the target function
   * @param fn target function to optimize
   *
   * Allocate the internal data.
   */
  void setTarget(ObjectFuncType* fn);

  /** optimize an object function
   * @param fn object function
   */
  bool optimize(ObjectFuncType* fn)
  {
    setTarget(fn);
    return optimize();
  }

  /** optimize TargetFunc */
  bool optimize();

  ///write to a std::ostream
  bool get(std::ostream& ) const;

  ///read from std::istream
  bool put(std::istream& );

  ///read from an xmlNode
  bool put(xmlNodePtr cur);

  ///reset member data
  void reset();

protected:

  int NumParams;
  int CurStep;

  Return_t curCost, prevCost;

  std::vector<Return_t> Y;
  std::vector<Return_t> gY, Y0;

  /** evaluate the value for y+dl*cg
   *
   * Lineminimization uses this function to find the minimum along the CG direction
   */
  //Return_t Func(Return_t dl);

  /** evaluate the gradients numerically
   * @param grad container for the gradients
   */
  void evaluateGradients(std::vector<Return_t>& grad);
};

template<class T>
DampedDynamics<T>::DampedDynamics(ObjectFuncType* atarget):
  TargetFunc(atarget), NumSteps(100),
  Displacement(1e-6),Dt(1.0),Friction(0.5),Mass(10.),
  CostTol(1.e-6),GradTol(1.e-6)
{
  prevCost=1e6; //some large number
}

template<class T>
void DampedDynamics<T>::setTarget(ObjectFuncType* fn)
{
  TargetFunc=fn;
  NumParams=TargetFunc->NumParams();
  Y.resize(NumParams);
  gY.resize(NumParams,0);
  Y0.resize(NumParams);
  for(int i=0; i<NumParams; i++)
  {
    Y[i]=TargetFunc->Params(i);
    Y0[i]=Y[i];
  }
}

template<class T>
bool DampedDynamics<T>::optimize()
{
  //T dt=1.0;
  //T friction=0.5;
  T t=Friction*Dt*0.5;
  T c0=2.0/(1+t);
  T c1=1-c0;
  T c2 =Dt*Dt/(1.0+t)/Mass;
  CurStep=0;
  do
  {
    evaluateGradients(gY);
    for(int i=0; i<NumParams; i++)
    {
      T old=Y0[i];
      Y0[i]=Y[i];
      Y[i]=c0*Y0[i]+c1*old+c2*gY[i];
    }
    curCost = TargetFunc->Cost();
    Return_t fx= std::abs(*(std::max_element(gY.begin(), gY.end())));
    if(fx<GradTol)
    {
      if(msg_stream)
        *msg_stream << " CGOptimization  has reached gradient max|G| = " << fx << std::endl;
      return false;
    }
    Return_t dx=std::abs((curCost-prevCost)/curCost);
    if(dx <= CostTol)
    {
      if(msg_stream)
        *msg_stream << " CGOptimization::Converged cost with " << dx << std::endl;
      return false;
    }
    CurStep++;
    if(TargetFunc->IsValid)
    {
      TargetFunc->Report();
    }
    else
    {
      if(msg_stream)
        *msg_stream << " DampedDynamics stopped due to invalid cost values " << std::endl;
      return false;
    }
  }
  while(CurStep<NumSteps);
  if(msg_stream)
    *msg_stream << " Failed to converged after " << CurStep << " steps." << std::endl;
  return false;
}

template<class T>
void
DampedDynamics<T>::evaluateGradients(std::vector<Return_t>& grad)
{
  TargetFunc->GradCost(grad, Y, Displacement);
  for(int i=0; i<grad.size(); i++)
    grad[i] *=-1.0;
//   if (Displacement==0)  TargetFunc->GradCost(grad, Y, Displacement);
//   else
//   {
//     Return_t dh=1.0/(2.0*Displacement);
//     for(int i=0; i<TargetFunc->NumParams() ; i++) {
//       for(int j=0; j<TargetFunc->NumParams(); j++) TargetFunc->Params(j)=Y[j];
//       TargetFunc->Params(i) = Y[i]+ Displacement;
//       Return_t CostPlus = TargetFunc->Cost();
//       TargetFunc->Params(i) = Y[i]- Displacement;
//       Return_t CostMinus = TargetFunc->Cost();
//       grad[i]=(CostMinus-CostPlus)*dh;
//     }
//   }
}

//template<class T>
//typename DampedDynamics<T>::Return_t
//DampedDynamics<T>::Func(Return_t dl) {
//  for(int i=0; i<NumParams; i++) TargetFunc->Params(i)=Y[i]+dl*cgY[i];
//  return TargetFunc->Cost();
//}


template<class T>
bool DampedDynamics<T>::get(std::ostream& os) const
{
  return true;
}

template<class T>
bool DampedDynamics<T>::put(std::istream& is)
{
  return true;
}

template<class T>
void DampedDynamics<T>::reset()
{
}

template<class T>
bool DampedDynamics<T>::put(xmlNodePtr cur)
{
  ParameterSet p;
  p.add(NumSteps,"max_steps","int");
  p.add(CostTol,"tolerance","scalar");
  p.add(GradTol,"tolerance_g","scalar");
  p.add(Displacement,"epsilon","scalar");
  p.add(Dt,"stepsize","scalar");
  p.add(Friction,"friction","scalar");
  p.add(Mass,"mass","scalar");
  p.put(cur);
  return true;
}
#endif


#ifndef QMCPLUSPLUS_COUNTING_FUNCTOR_H
#define QMCPLUSPLUS_COUNTING_FUNCTOR_H

#include "OhmmsData/AttributeSet.h"

namespace qmcplusplus
{

//template <class T> class SigmoidFunctor
//{
//  public:
//
//  SigmoidFunctor() {}
//  bool put(xmlNodePtr cur);
//
//
//
//};


// doesn't inherit from OptimizableFunctorBase since this is a function of the entire position vector
template <class T> class GaussianFunctor: public QMCTraits
{
  typedef optimize::VariableSet::real_type real_type;
  typedef optimize::VariableSet opt_variables_type;


  // id string
  std::string id;

  // most recent evaluations
  RealType Fval;
  GradType Fgrad;
  RealType Flap; 

public: 
  // optimizable variables
  opt_variables_type myVars;

  const std::vector<std::pair<std::string, real_type> >* getNameAndValue() { return &myVars.NameAndValue; }

  GaussianFunctor(std::string fid)
  {
    id = fid;
  }


  void restore(int iat) {}

  void checkInVariables(opt_variables_type& active)
  { 
    active.insertFrom(myVars); 
  }


  void checkOutVariables(const opt_variables_type& active)
  {
    myVars.getIndex(active); 
  }


  void resetParameters(const opt_variables_type& active)
  {
  }

  void reportStatus(std::ostream& os) 
  {
  }

  GaussianFunctor<T>* makeClone() const
  {
  }


  bool put(xmlNodePtr cur)
  {
    return true;
  }


  void multiply_eq(const GaussianFunctor* rhs)
  {
  }


  void divide_eq(const GaussianFunctor* rhs)
  {
  }


  void evaluate(PosType r, RealType& fval, GradType& fgrad, RealType& flap)
  {
  }

  void evaluateLog(PosType r, RealType& lval, GradType& lgrad, RealType& llap)
  {
  }

  void evaluateDerivatives(PosType r, std::vector<RealType>& dfval, std::vector<GradType>& dfgrad, std::vector<RealType>& dflap)
  {
  }

  void evaluateLogDerivatives(PosType r, std::vector<RealType>& dlval, std::vector<GradType>& dlgrad, std::vector<RealType>& dllap)
  {
  }

  void evaluateLogTempDerivatives(PosType r, std::vector<RealType>& dlval)
  {
  }

  void evaluate_print(std::ostream& os, ParticleSet& P)
  {
  }


};

}
#endif

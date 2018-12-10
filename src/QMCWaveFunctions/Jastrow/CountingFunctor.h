
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

  // enumerations for axes and T parameters
  enum A_vars { XX, XY, XZ, YY, YZ, ZZ, NUM_A };
  enum B_vars { X, Y, Z, DIM};

  //TensorType A;
  //PosType B;
  //RealType C;

  // opt variables: vector of bools: one for each parameter
  std::vector<bool> opt_A;
  std::vector<bool> opt_B;
  bool opt_C;

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

  GaussianFunctor<T>* makeClone(std::string fid) const
  {
    app_log() << "  GaussianCountingFunctor::makeClone" << std::endl;
    GaussianFunctor<T>* rptr = new GaussianFunctor<T>(fid);
    // copy class variables defined in put
    //std::copy(A.begin(),A.end(), rptr->A.begin());
    //std::copy(B.begin(),B.end(), rptr->B.begin());
//    for(int i = 0; i < A.size(); ++i)
//      rptr->A[i] = A[i];
//    for(int i = 0; i < B.size(); ++i)
//      rptr->B[i] = B[i];
//    rptr->C = C;
//    rptr->opt_A.resize(opt_A.size());
//    rptr->opt_B.resize(opt_B.size());
//    for(int i = 0; i < opt_A.size(); ++i)
//      rptr->opt_A[i] = opt_A[i];
//    for(int i = 0; i < opt_B.size(); ++i)
//      rptr->opt_B[i] = opt_B[i];
//    //std::copy(opt_A.begin(), opt_A.end(), rptr->opt_A.begin());
//    //std::copy(opt_B.begin(), opt_B.end(), rptr->opt_B.begin());
//    rptr->opt_C = opt_C;
//    //rptr->initialize();
//    rptr->myVars = myVars;
    return rptr;
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

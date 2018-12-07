#ifndef QMC_PLUS_PLUS_COUNTING_FUNCTOR_H
#define QMC_PLUS_PLUS_COUNTING_FUNCTOR_H

#include "CountingFunctor.h"
#include "Particle/ParticleSet.h"

namespace qmcplusplus
{

// T is precision

template <class T> class NormalizedGaussianRegion
{
public:
  typedef GaussianFunctor<T> FunctorType;

  // variables

  // constructor
  NormalizedGaussianRegion(ParticleSet& targetPtcl)
  {
  }

  // destructor
  ~NormalizedGaussianRegion()
  {}

  void addFunc(FunctorType* func, std::string id);
};

template <class T> class SigmoidRegion
{
public:
  typedef SigmoidFunctor<T> FunctorType;

  // variables

  // constructor
  SigmoidRegion(ParticleSet& targetPtcl)
  {
  }

  // destructor
  ~SigmoidRegion()
  {}

  void addFunc(FunctorType* func, std::string id);
  //void addFunc( );
};

}
#endif

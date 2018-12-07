

#include "CountingFunctor.h"

// T is precision

template <class T> class NormalizedGaussianRegion
{
public:
  typedef GaussianFunctor FunctorType;

  // variables

  // constructor
  NormalizedGaussianRegion()
  {
  }

  // destructor
  ~NormalizedGaussianRegion()
  {}

  void addFunc(*FT, std::string id);
}

template <class T> class SigmoidRegion
{
public:
  typedef SigmoidFunctor FunctorType<T>;

  // variables

  // constructor
  SigmoidRegion()
  {
  }

  // destructor
  ~SigmoidRegion()
  {}

  void addFunc(*FT, std::string id);
public:
  //void addFunc( );
}

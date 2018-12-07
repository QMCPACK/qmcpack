

// T is precision

template <class T> class NormalizedGaussianRegion
{
public:
  typedef GaussianFunctor FT;

  // variables

  // constructor
  NormalizedGaussianRegion()
  {
  }

  // destructor
  ~NormalizedGaussianRegion()
  {}

  void addFunc(*FT);
}

template <class T> class SigmoidRegion
{
public:
//  typedef GaussianFunctor FT;

  // variables

  // constructor
  SigmoidRegion()
  {
  }

  // destructor
  ~SigmoidRegion()
  {}

  void addFunc(*FT);
public:
  //void addFunc( );
}

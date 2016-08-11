#ifndef FAKERANDOM_H
#define FAKERANDOM_H

// Deterministic generator for testing.

class FakeRandom
{
public:
  typedef unsigned int uint_type;
  double operator()();
};

#endif

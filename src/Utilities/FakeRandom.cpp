
#include "Utilities/FakeRandom.h"

double FakeRandom::operator()()
{
  // Decided by committee to be 'random enough'.
  return 0.5;
}


#ifndef YLM_H
#define YLM_H

#include <cmath>
#include <complex>
#include "../OhmmsPETE/TinyVector.h"

using namespace std;

namespace qmcplusplus {
  double LegendrePll (int l, double x);

  double LegendrePlm (int l, int m, double x);

  complex<double> Ylm(int l, int m, TinyVector<double,3> r);
}

#endif

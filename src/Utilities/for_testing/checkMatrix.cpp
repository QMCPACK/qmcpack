
#include "checkMatrix.hpp"

namespace qmcplusplus
{
  template bool approxEquality<double>(double val_a, double val_b);
  template bool approxEquality<std::complex<double>>(std::complex<double> val_a, std::complex<double> val_b);
}

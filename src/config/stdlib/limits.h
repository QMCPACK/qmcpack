#ifndef QMCPLUSPLUS_LIMITS_H
#define QMCPLUSPLUS_LIMITS_H
#include <limits>
template<typename T>
inline bool iszero(T a)
{
  return (abs(a)<numeric_limits<T>::epsilon());
}
#endif


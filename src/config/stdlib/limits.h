#ifndef QMCPLUSPLUS_LIMITS_H
#define QMCPLUSPLUS_LIMITS_H
#include <limits>
template<typename T>
inline bool iszero(T a)
{
  return (std::abs(a)<std::numeric_limits<T>::epsilon());
}
#endif


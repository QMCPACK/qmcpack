#include <cmath>

namespace qmcplusplus
{
bool isnan(float a) { return a != a; }
bool isnan(double a) { return a != a; }
} // namespace qmcplusplus

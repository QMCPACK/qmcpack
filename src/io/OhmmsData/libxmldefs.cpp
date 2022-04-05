#include "libxmldefs.h"
#include "ModernStringUtils.hpp"

std::string getNodeName(xmlNodePtr cur)
{
  return qmcplusplus::lowerCase(castXMLCharToChar(cur->name));
}

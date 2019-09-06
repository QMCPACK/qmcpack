
//Prototype code for an engine to handle hybrid optimization


#ifndef QMCPLUSPLUS_HYBRID_ENGINE_HEADER
#define QMCPLUSPLUS_HYBRID_ENGINE_HEADER

#include <vector>
#include <libxml/tree.h>
#include "Message/Communicate.h"
#include "OhmmsData/ParameterSet.h"


namespace qmcplusplus
{
class HybridEngine
{
private:
  Communicate* myComm;

  ParameterSet m_param;


  //number of steps in a descent section of hybrid method
  int descent_len;
  //number of steps in a BLM section of hybrid method
  int blm_len;

  //Vector for storing parameter values for calculating differences to be given to hybrid method
  std::vector<double> paramsForDiff;

  ///process xml node
  bool processXML(const xmlNodePtr cur);

public:
  //Constructor for engine
  HybridEngine(Communicate* comm, const xmlNodePtr cur);
};

} // namespace qmcplusplus
#endif

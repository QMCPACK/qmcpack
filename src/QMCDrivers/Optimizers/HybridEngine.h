
//Prototype code for an engine to handle descent optimization


#ifndef QMCPLUSPLUS_DESCENT_ENGINE_HEADER
#define QMCPLUSPLUS_DESCENT_ENGINE_HEADER

#include <vector>
#include "Message/Communicate.h"



namespace qmcplusplus
{
class HybridEngine
{
private:
  
  Communicate* myComm;

 

//Vector for storing parameter values for calculating differences to be given to hybrid method
  std::vector<double> paramsForDiff;

public:
  //Constructor for engine
  HybridEngine(const bool targetExcited, Communicate* comm);


  ///process xml node
  //bool processXML(xmlNodePtr cur);
  

};

} // namespace qmcplusplus
#endif


//Prototype code for an engine to handle hybrid optimization


#ifndef QMCPLUSPLUS_HYBRID_ENGINE_HEADER
#define QMCPLUSPLUS_HYBRID_ENGINE_HEADER

#include <vector>
#include <libxml/tree.h>
#include "Message/Communicate.h"
#include "Optimize/VariableSet.h"

namespace qmcplusplus
{
class HybridEngine
{
private:
  Communicate* myComm;



  //number of steps in a descent section of hybrid method
  int descent_len;
  //number of steps in a BLM section of hybrid method
  int blm_len;

  //Vector for storing parameter values for calculating differences to be given to hybrid method
  std::vector<double> paramsForDiff;

  //Vector for storing the input vectors to the BLM steps of hybrid method
  std::vector<std::vector<double>> hybridBLM_Input;

  ///process xml node
  bool processXML(const xmlNodePtr cur);

public:
  //Constructor for engine
  HybridEngine(Communicate* comm, const xmlNodePtr cur);


void getInitialParams(const optimize::VariableSet& myVars);

void storeVectors(std::vector<double>& currentParams,int descentCount);

const std::vector<std::vector<double>> retrieveHybridBLM_Input() const {return  hybridBLM_Input;}

const int getDescentLen() const {return descent_len;}

const int getBLMLen() const {return blm_len;}


};

} // namespace qmcplusplus
#endif

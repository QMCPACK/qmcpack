
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

  ///process xml node
  bool processXML(const xmlNodePtr cur);

  ///xml saved node
  std::vector<xmlNodePtr> saved_xml_opt_methods_;

  //number of updates in each individual method of hybrid method
  std::vector<int> num_updates_opt_methods_;

  //inidividual methods used in hybrid optimization
  std::vector<std::string> saved_opt_method_types_;

  int identifyMethodIndex(int counter) const;

public:
  //Constructor for engine
  HybridEngine(Communicate* comm, const xmlNodePtr cur);

  xmlNodePtr getSelectedXML(int counter) const;

  const int getDescentLen() const { return descent_len; }

  const int getBLMLen() const { return blm_len; }

  const bool queryStore(int counter,int store_num,std::string methodType) const;

  const std::string queryMethod(int counter);

};

} // namespace qmcplusplus
#endif

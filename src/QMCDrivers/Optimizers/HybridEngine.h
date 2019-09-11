
//Prototype code for an engine to handle hybrid optimization


#ifndef QMCPLUSPLUS_HYBRID_ENGINE_HEADER
#define QMCPLUSPLUS_HYBRID_ENGINE_HEADER

#include <vector>
#include <libxml/tree.h>
#include "Message/Communicate.h"
#include "Optimize/VariableSet.h"
#include "QMCDrivers/Optimizers/OptimizerTypes.h"


namespace qmcplusplus
{
class HybridEngine
{
private:
  Communicate* myComm;


  ///number of optimization steps taken
  int step_num_;

  ///process xml node
  bool processXML(const xmlNodePtr cur);

  ///list of methods known by hybrid
  std::vector<OptimizerType> opt_methods_;

  ///xml saved node
  std::vector<xmlNodePtr> saved_xml_opt_methods_;

  //number of updates in each individual method of hybrid method
  std::vector<int> num_updates_opt_methods_;

  int identifyMethodIndex() const;

public:
  //Constructor for engine
  HybridEngine(Communicate* comm, const xmlNodePtr cur);

  xmlNodePtr getSelectedXML();

  bool queryStore(int store_num, OptimizerType methodType) const;
};

} // namespace qmcplusplus
#endif



#include <vector>
#include <string>

#include "QMCDrivers/Optimizers/HybridEngine.h"
#include "OhmmsData/ParameterSet.h"
#include "Message/CommOperators.h"


namespace qmcplusplus
{
HybridEngine::HybridEngine(Communicate* comm, const xmlNodePtr cur) : myComm(comm)
{
  processXML(cur);
}


bool HybridEngine::processXML(const xmlNodePtr cur)
{
    ParameterSet m_param;
  //Number of steps in a descent section of hybrid method
  m_param.add(descent_len, "descent_length", "int");
  //Number of steps in a BLM of section of hybrid method
  m_param.add(blm_len, "BLM_length", "int");

  m_param.put(cur);


  return true;
}

void HybridEngine::getInitialParams(const optimize::VariableSet& myVars)
{

    for(int i =0; i < myVars.size(); i++)
    {
	paramsForDiff.push_back(myVars[i]);
    }

}


// Helper method for storing vectors of parameter differences over the course of
// a descent optimization for use in BLM steps of the hybrid method
void HybridEngine::storeVectors(std::vector<double>& currentParams,int descentCount)
{


  std::vector<double> rowVec(currentParams.size());
  std::fill(rowVec.begin(), rowVec.end(), 0.0);

  // Take difference between current parameter values and the values from 20
  // iterations before (in the case descent_len = 100) to be stored as input to BLM.
  // The current parameter values are then copied to paramsForDiff to be used
  // another 20 iterations later.
  for (int i = 0; i < currentParams.size(); i++)
  {
    rowVec[i]        = currentParams[i] - paramsForDiff[i];
    paramsForDiff[i] = currentParams[i];
  }

  // If on first step, clear anything that was in vector
  if ((descentCount + 1) % descent_len == descent_len / 5)
  {
    hybridBLM_Input.clear();
    hybridBLM_Input.push_back(rowVec);
  }
  else
  {
    hybridBLM_Input.push_back(rowVec);
  }

  for (int i = 0; i < hybridBLM_Input.size(); i++)
  {
    std::string entry = "";
    for (int j = 0; j < hybridBLM_Input.at(i).size(); j++)
    {
      entry = entry + std::to_string(hybridBLM_Input.at(i).at(j)) + ",";
    }
    app_log() << "Stored Vector: " << entry << std::endl;
  }
}


} // namespace qmcplusplus

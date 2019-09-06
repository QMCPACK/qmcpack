

#include <vector>
#include <string>

#include "QMCDrivers/Optimizers/HybridEngine.h"
#include "Message/CommOperators.h"


namespace qmcplusplus
{
HybridEngine::HybridEngine(Communicate* comm,const xmlNodePtr cur)
    : myComm(comm)   
{

//Number of steps in a descent section of hybrid method
    m_param.add(descent_len,"descent_length","int");
//Number of steps in a BLM of section of hybrid method
m_param.add(blm_len,"BLM_length","int");
processXML(cur);
}


bool HybridEngine::processXML(const xmlNodePtr cur)
{
m_param.put(cur);


return true;

}



} // namespace qmcplusplus

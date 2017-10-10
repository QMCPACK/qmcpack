
#include "OhmmsData/libxmldefs.h"
#include "OhmmsData/AttributeSet.h"
#include "OhmmsData/ParameterSet.h"

#include "AFQMC/config.h"
#include<Message/MPIObjectBase.h>
#include "AFQMC/Wavefunctions/WavefunctionHandler.h"
#include "AFQMC/Wavefunctions/WavefunctionBase.h"
#include "AFQMC/Wavefunctions/PureSingleDeterminant.h"
#include "AFQMC/Wavefunctions/GeneralSingleDeterminant.h"
#include "AFQMC/Wavefunctions/MultiPureSingleDeterminant.h"

namespace qmcplusplus
{

bool WavefunctionHandler::parse(xmlNodePtr cur)
{

  if(cur == NULL)
    return false;

  std::string str("no");
  xmlNodePtr curRoot=cur;
  OhmmsAttributeSet oAttrib;
  oAttrib.add(name,"name");
  oAttrib.add(str,"new_algo");
  oAttrib.put(cur);

  if(str=="yes") new_algo=true;

  cur = curRoot->children;
  while (cur != NULL) {
    std::string cname((const char*)(cur->name));
    if(cname =="ImpSamp") {
      ImpSampWfn = addWfn(cur);
      ImpSampWfn->wfn_role += std::string("ImpSamp");
    } else if(cname =="Estimator") {
      EstimatorWfn = addWfn(cur); 
      EstimatorWfn->wfn_role += std::string("Estimator");
app_log()<<" Successfully created Estimator wfn: " <<EstimatorWfn->name <<" " <<EstimatorWfn->wfn_role <<std::endl;
    //} else if(cname =="phaseless") {
    //  phaselessWfn = addWfn(cur); 
    }    

    cur = cur->next;
  }

  return true;
  
} 

WavefunctionHandler::WfnPtr WavefunctionHandler::addWfn(xmlNodePtr cur) 
{

  if(cur==NULL) return NULL;

  std::string id("");
  std::string type("");
  OhmmsAttributeSet aAttrib;
  aAttrib.add(type,"type");
  aAttrib.add(id,"name");
  aAttrib.add(id,"id");
  aAttrib.put(cur);

  // look for object with same id
  for(int i=0; i<wfns.size(); i++)
    if(id == wfns[i]->name) return wfns[i]; 

  WfnPtr it;
  if(type == "PureSD") 
    it = (WavefunctionBase*) new PureSingleDeterminant(myComm);
  else if(type == "GenSD") 
    it = (WavefunctionBase*) new GeneralSingleDeterminant(myComm);
  else if(type == "MultiPureSD")
    it = (WavefunctionBase*) new MultiPureSingleDeterminant(myComm);
  
  it->parse(cur);
  wfns.push_back(it);
  return it;

}
 

bool WavefunctionHandler::init(std::vector<int>& TGdata, SPComplexSMVector *v, hdf_archive& read, const std::string& tag, MPI_Comm tg_comm, MPI_Comm node_comm, MPI_Comm node_heads_comm)
{
  MPI_COMM_HEAD_OF_NODES = node_heads_comm;
  head_of_nodes = (TGdata[1]==0);

  for(int i=0; i<wfns.size(); i++) {
    wfns[i]->copyInfo(*this);
    if(!wfns[i]->init(TGdata,v,read,tag,tg_comm,node_comm,node_heads_comm))
      return false;
  } 
  return true;
} 


bool WavefunctionHandler::setup(HamPtr ham)
{
  for(int i=0; i<wfns.size(); i++)
    if(!wfns[i]->setup(ham))
      return false;
  return true;
} 


} // qmcplusplus

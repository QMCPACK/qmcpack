//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



/**@file QMCInterface.cpp
 * @brief Implments QMCInterface operators.
 */
#include "QMCApp/QMCInterface.h"
#include "QMCApp/ParticleSetPool.h"
#include "QMCApp/WaveFunctionPool.h"
#include "QMCApp/HamiltonianPool.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "QMCHamiltonians/QMCHamiltonian.h"
#include "Particle/HDFWalkerIO.h"
#include "QMCApp/InitMolecularSystem.h"
#include "Particle/DistanceTable.h"
#include "QMCDrivers/QMCDriver.h"
#include "QMCDrivers/VMC/VMCSingle.h"
#include "QMCDrivers/CorrelatedSampling/CSVMC.h"
//#include "QMCDrivers/VMC/VMCMultiple.h"
#include "QMCDrivers/RQMCMultiple.h"
#include "Message/Communicate.h"
#include <queue>
#include "OhmmsData/AttributeSet.h"
#include <sstream>

namespace qmcplusplus
{

typedef std::map<std::string,HamiltonianFactory*> PoolType;

QMCInterface::QMCInterface(Communicate* c): QMCDriverFactory(c), QMCAppBase(),
  FirstQMC(true)
{
  app_log() << "\n=========================================================\n"
            <<   "                   qmcplusplus 0.2                       \n"
            << "\n  (c) Copyright 2003-  qmcplusplus developers          \n"
            <<   "=========================================================\n";
  app_log().flush();
}

///destructor
QMCInterface::~QMCInterface()
{
}


bool QMCInterface::initialize(int myProc, int numProcs)
{
  if(XmlDocStack.empty())
  {
    ERRORMSG("No valid input file exists! Aborting QMCInterface::initialize")
    return false;
  }
  OHMMS::Controller->setNodeID(myProc);
  OHMMS::Controller->setNumNodes(numProcs);
  std::ostringstream newTitle;
  newTitle << myProject.m_title << "." << OHMMS::Controller->rank();
  myProject.setTitle(newTitle.str());
  //validate the input file
  bool success = validateXML();
  if(!success)
  {
    ERRORMSG("Input document does not contain valid objects")
    return false;
  }
  //initialize all the instances of distance tables and evaluate them
  ptclPool->reset();
  OHMMS::Controller->barrier();
  //write stuff
  app_log() << "=========================================================\n";
  app_log() << " QMC system initialized:\n";
  app_log() << "=========================================================\n";
  ptclPool->get(app_log());
  hamPool->get(app_log());
  return true;
}

bool QMCInterface::SetVMC(double dt, int w, int steps, int nblocks)
{
  //cerr << "    QMCInterface::SetVMC: Creating new VMC driver...";
  if(qmcDriver != NULL)
  {
    //cerr << "qmcpack: deleting previous driver" << std::endl;
    delete qmcDriver;
  }
  qmcDriver = new VMCSingle(*ptclPool->getWalkerSet("e"),*psiPool->getPrimary(),*hamPool->getPrimary(),*psiPool);
  //cerr << " done." << std::endl;
  bool append_run = false;
  qmcDriver->setStatus(myProject.CurrentRoot(),PrevConfigFile, append_run);
  //cerr << "QMCInterface::SetVMC: Using Default Parameters:" << std::endl;
  //cerr << "  timestep = 0.03; walkers = 100; steps = 100;" << std::endl;
  qmcDriver->setValue("timeStep", dt);
  qmcDriver->setValue("walkers", w);
  qmcDriver->setValue("steps", steps);
  qmcDriver->setValue("blocks", nblocks);
  //cerr << "  blocks = " << nblocks << std::endl;
  //cerr << "  Setting xmlNodepointer...";
  runInfoNode = xmlNewNode(NULL,(const xmlChar*)"vmc");
  //cerr << " Done.  Ready to run." << std::endl;
  return true;
}

bool QMCInterface::SetVMCMultiple(double dt, int w, int steps, int nblocks)
{
  //cerr << "    QMCInterface::SetVMCMultiple: Creating new VMCMultiple driver...";
  if(qmcDriver != NULL)
  {
    //cerr << "qmcpack: deleting previous driver" << std::endl;
    delete qmcDriver;
  }
  qmcDriver = new CSVMC(*ptclPool->getWalkerSet("e"),*psiPool->getPrimary(),*hamPool->getPrimary(),*psiPool);
  //qmcDriver = new VMCMultiple(*ptclPool->getWalkerSet("e"),*psiPool->getPrimary(),*hamPool->getPrimary());
  //cerr << " done." << std::endl;
  // get second psi, hamiltonian
  QMCHamiltonian* secondHam = hamPool->getHamiltonian("h1");
  TrialWaveFunction* secondPsi = psiPool->getWaveFunction("psi1");
  // add them
  qmcDriver->add_H_and_Psi(secondHam,secondPsi);
  bool append_run = false;
  qmcDriver->setStatus(myProject.CurrentRoot(),PrevConfigFile, append_run);
  qmcDriver->setValue("timeStep", dt);
  qmcDriver->setValue("walkers", w);
  qmcDriver->setValue("steps", steps);
  qmcDriver->setValue("blocks", nblocks);
  runInfoNode = xmlNewNode(NULL,(const xmlChar*)"vmc-multi");
  return true;
}

bool QMCInterface::SetRQMCMultiple(double dt, int chains, int steps, int nblocks)
{
  bool append_run;
  bool isNewDriver = false;
  if(qmcDriver == NULL)
  {
    //cerr << "Creating new RQMC driver" << std::endl;
    qmcDriver = new RQMCMultiple(*ptclPool->getWalkerSet("e"),*psiPool->getPrimary(),*hamPool->getPrimary(),*psiPool);
    // get second psi, hamiltonian
    QMCHamiltonian* secondHam = hamPool->getHamiltonian("h1");
    TrialWaveFunction* secondPsi = psiPool->getWaveFunction("psi1");
    // add them
    qmcDriver->add_H_and_Psi(secondHam,secondPsi);
    append_run = false;
    isNewDriver = true;
  }
  else
  {
    //cerr << "REUSING RQMC driver" << std::endl;
    append_run = false;
  }
  qmcDriver->setStatus(myProject.CurrentRoot(),PrevConfigFile, append_run);
  qmcDriver->setValue("timeStep", dt);
  qmcDriver->setValue("chains", chains);
  qmcDriver->setValue("steps", steps);
  qmcDriver->setValue("blocks", nblocks);
  runInfoNode = xmlNewNode(NULL,(const xmlChar*)"rqmc-multi");
  //cerr << " done." << std::endl;
  //return true;
  return isNewDriver;
}

/*
   bool QMCInterface::SetRQMCMultiple(double dt, int chains, int steps, int nblocks){
   if(qmcDriver != NULL){
   delete qmcDriver;
   }
   qmcDriver = new RQMCMultiple(*ptclPool->getWalkerSet("e"),*psiPool->getPrimary(),*hamPool->getPrimary());
//cerr << " done." << std::endl;

// get second psi, hamiltonian
QMCHamiltonian* secondHam = hamPool->getHamiltonian("h1");

TrialWaveFunction* secondPsi = psiPool->getWaveFunction("psi1");

// add them
qmcDriver->add_H_and_Psi(secondHam,secondPsi);

bool append_run = false;
qmcDriver->setStatus(myProject.CurrentRoot(),PrevConfigFile, append_run);

qmcDriver->setValue("timeStep", dt);
qmcDriver->setValue("chains", chains);
qmcDriver->setValue("steps", steps);
qmcDriver->setValue("blocks", nblocks);

runInfoNode = xmlNewNode(NULL,(const xmlChar*)"rqmc-multi");
return true;
}
*/

bool QMCInterface::process()
{
  //cerr << "    QMCInterface::process: Processing..." << std::endl;
  qmcDriver->process(runInfoNode);
  //DO NOT NEED THIS
  //qmcDriver->Estimators->setAccumulateMode(true);
  //cerr << " done." << std::endl;
  return true;
}

bool QMCInterface::execute()
{
  // SERIOUSLY SPECIALIZED HACK FOR H2 WITH MOLECULE-CENTERED GAUSSIAN ORBIATLS
  //ParticleSet* refSet = ptclPool->getParticleSet("i");
  //ParticleSet* pseudoSet = ptclPool->getParticleSet("pseudo");
  //ParticleSet* refSet1 = ptclPool->getParticleSet("i1");
  //ParticleSet* pseudoSet1 = ptclPool->getParticleSet("pseudo1");
  //cerr << "	qmcpack: updating pseudo center from " << pseudoSet->R[0];
  //pseudoSet->R[0] = 0.5*(refSet->R[0] + refSet->R[1]);
  //cerr << " to " << pseudoSet->R[0] << " given H1 at " << refSet->R[0] << " and H2 at " << refSet->R[1] << std::endl;
  //pseudoSet1->R[0] = 0.5*(refSet1->R[0] + refSet1->R[1]);
  qmcDriver->run();
  return true;
}

// returns a pointer to a vector containing energy estimator output
/*
   std::vector<double>* QMCInterface::GetData(){
   return qmcDriver->energyData.Send();
   }
   */

/*
   double QMCInterface::GetData( std::string estimator, std::string tag){
   int EstIndex = qmcDriver->Estimators.add(estimator);
   int ValueIndex = qmcDriver->Estimators[EstIndex].add(tag);
   double value = qmcDriver->Estimators[EstIndex].average(ValueIndex);
   return value;
   }
   */
// sets ptcl coordinates in set "i" for ion
void QMCInterface::SetPtclPos(int id, double* newR)
{
  ParticleSet* mySet = ptclPool->getParticleSet("i");
  //cerr << "    QMCInterface::SetPtclPos: updated ion " << id << " from " << mySet->R[id];
  qmcplusplus::TinyVector<double,3> newVec(newR[0],newR[1],newR[2]);
  mySet->R[id] = newVec;
  //cerr << " to " << mySet->R[id] << std::endl;
}

// sets ptcl coordinates in speciefied set
void QMCInterface::SetPtclPos( std::string set, int id, double* newR)
{
  //cerr << "int SetPtclPos, looking to move " << id << " of set " << set << std::endl;
  ParticleSet* mySet = ptclPool->getParticleSet(set);
  //cerr << "    QMCInterface::SetPtclPos: updated particle " << id << " of ParticleSet " << set << " from " << mySet->R[id];
  qmcplusplus::TinyVector<double,3> newVec(newR[0],newR[1],newR[2]);
  mySet->R[id] = newVec;
  //cerr << " to " << mySet->R[id] << std::endl;
}

/** validate the Interface document
 * @return false, if any of the basic objects is not properly created.
 *
 * Current xml schema is changing. Instead validating the input file,
 * we use xpath to create necessary objects. The order follows
 * - project: determine the simulation title and sequence
 * - random: initialize random number generator
 * - particleset: create all the particleset
 * - wavefunction: create wavefunctions
 * - hamiltonian: create hamiltonians
 * Finally, if /simulation/mcwalkerset exists, read the configurations
 * from the external files.
 */
bool QMCInterface::validateXML()
{
  xmlXPathContextPtr m_context = XmlDocStack.top()->getXPathContext();
  OhmmsXPathObject result("//project",m_context);
  if(result.empty())
  {
    app_warning() << "Project is not defined" << std::endl;
    myProject.reset();
  }
  else
  {
    myProject.put(result[0]);
  }
  app_log() << std::endl;
  myProject.get(app_log());
  app_log() << std::endl;
  //initialize the random number generator
  xmlNodePtr rptr = myRandomControl.initialize(m_context);
  //preserve the input order
  xmlNodePtr cur=XmlDocStack.top()->getRoot()->children;
  while(cur != NULL)
  {
    std::string cname((const char*)cur->name);
    if(cname == "particleset")
    {
      ptclPool->put(cur);
    }
    else
      if(cname == "wavefunction")
      {
        psiPool->put(cur);
      }
      else
        if(cname == "hamiltonian")
        {
          hamPool->put(cur);
        }
        else
          if(cname == "include")
            //file is provided
          {
            const xmlChar* a=xmlGetProp(cur,(const xmlChar*)"href");
            if(a)
            {
              pushDocument((const char*)a);
              processPWH(XmlDocStack.top()->getRoot());
              popDocument();
            }
          }
          else
            if(cname == "qmcsystem")
            {
              processPWH(cur);
            }
            else
              if(cname == "init")
              {
                InitMolecularSystem moinit(ptclPool);
                moinit.put(cur);
              }
    cur=cur->next;
  }
  if(ptclPool->empty())
  {
    ERRORMSG("Illegal input. Missing particleset ")
    return false;
  }
  if(psiPool->empty())
  {
    ERRORMSG("Illegal input. Missing wavefunction. ")
    return false;
  }
  if(hamPool->empty())
  {
    ERRORMSG("Illegal input. Missing hamiltonian. ")
    return false;
  }
  setMCWalkers(m_context);
  return true;
}



/** grep basic objects and add to Pools
 * @param cur current node
 *
 * Recursive search  all the xml elements with particleset, wavefunction and hamiltonian
 * tags
 */
void QMCInterface::processPWH(xmlNodePtr cur)
{
  if(cur == NULL)
    return;
  cur=cur->children;
  while(cur != NULL)
  {
    std::string cname((const char*)cur->name);
    if(cname == "simulationcell")
    {
      ptclPool->putLattice(cur);
    }
    else
      if(cname == "particleset")
      {
        ptclPool->put(cur);
      }
      else
        if(cname == "wavefunction")
        {
          psiPool->put(cur);
        }
        else
          if(cname == "hamiltonian")
          {
            hamPool->put(cur);
          }
    cur=cur->next;
  }
}

/** prepare for a QMC run
 * @param cur qmc element
 * @return true, if a valid QMCDriver is set.
 */
bool QMCInterface::runQMC(xmlNodePtr cur)
{
  OHMMS::Controller->barrier();
  bool append_run = setQMCDriver(myProject.m_series,cur);
  if(qmcDriver)
  {
    app_log() << std::endl;
    myProject.get(app_log());
    app_log() << std::endl;
    //advance the project id
    //if it is NOT the first qmc node and qmc/@append!='yes'
    if(!FirstQMC && !append_run)
      myProject.advance();
    qmcDriver->setStatus(myProject.CurrentRoot(),PrevConfigFile, append_run);
    qmcDriver->putWalkers(m_walkerset_in);
    qmcDriver->process(cur);
    //set up barrier
    OHMMS::Controller->barrier();
    qmcDriver->run();
    //keeps track of the configuration file
    PrevConfigFile = myProject.CurrentRoot();
    //do not change the input href
    //change the content of mcwalkerset/@file attribute
    for(int i=0; i<m_walkerset.size(); i++)
    {
      xmlSetProp(m_walkerset[i],
                 (const xmlChar*)"href", (const xmlChar*)myProject.CurrentRoot());
    }
    //curMethod = what;
    return true;
  }
  else
  {
    return false;
  }
}

bool QMCInterface::setMCWalkers(xmlXPathContextPtr context_)
{
  OhmmsXPathObject result("/simulation/mcwalkerset",context_);
  if(result.empty())
  {
    if(m_walkerset.empty())
    {
      result.put("//qmc",context_);
      xmlNodePtr newnode_ptr = xmlNewNode(NULL,(const xmlChar*)"mcwalkerset");
      if(result.empty())
      {
        xmlAddChild(XmlDocStack.top()->getRoot(),newnode_ptr);
      }
      else
      {
        xmlAddPrevSibling(result[0],newnode_ptr);
      }
      m_walkerset.push_back(newnode_ptr);
    }
  }
  else
  {
    for(int iconf=0; iconf<result.size(); iconf++)
    {
      xmlNodePtr mc_ptr = result[iconf];
      m_walkerset.push_back(mc_ptr);
      m_walkerset_in.push_back(mc_ptr);
    }
  }
  return true;
}
}

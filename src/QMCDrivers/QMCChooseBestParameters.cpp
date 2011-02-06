//////////////////////////////////////////////////////////////////
// (c) Copyright 2005- by Jeongnim Kim                                                
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim                                                                     
//   National Center for Supercomputing Applications &                                
//   Materials Computation Center                                                     
//   University of Illinois, Urbana-Champaign                                         
//   Urbana, IL 61801                                                                 
//   e-mail: jnkim@ncsa.uiuc.edu                                                      
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)                                   
//                                                                                    
// Supported by                                                                       
//   National Center for Supercomputing Applications, UIUC                            
//   Materials Computation Center, UIUC                                               
//////////////////////////////////////////////////////////////////                    
// -*- C++ -*-                                                                        
#include "QMCDrivers/QMCChooseBestParameters.h"                                             
#include "OhmmsData/AttributeSet.h"                                                   
#include "OhmmsData/ParameterSet.h"  
#include "Message/CommOperators.h"                                                    
#if defined(ENABLE_OPENMP)                                                            
#include "QMCDrivers/VMC/VMCSingleOMP.h"                                              
#include "QMCDrivers/QMCCostFunctionOMP.h"                                            
#endif                                                                                
#include "QMCDrivers/VMC/VMCSingle.h"                                                 
#include "QMCDrivers/QMCCostFunctionSingle.h"                                         
#include "Numerics/Blasf.h"                                                           
#include <cassert>   
#if defined(QMC_CUDA)
  #include "QMCDrivers/VMC/VMC_CUDA.h"
  #include "QMCDrivers/QMCCostFunctionCUDA.h"
#endif

                                                                 
namespace qmcplusplus                                                                 
{
  
  QMCChooseBestParameters::QMCChooseBestParameters(MCWalkerConfiguration& w,
    TrialWaveFunction& psi, QMCHamiltonian& h, HamiltonianPool& hpool, WaveFunctionPool& ppool): QMCDriver(w,psi,h,ppool), CloneManager(hpool), 
     hamPool(hpool),vmcEngine(0), WF(&psi), WarmupBlocks(10), alpha(1)
    {
      //set the optimization flag
      QMCDriverMode.set(QMC_OPTIMIZE,1);
      //read to use vmc output (just in case)                                                                                   
      RootName = "pot";                                                                                    
      QMCType ="QMCChooseBestParameters";   
    }
    
    /** Clean up the vector */
    QMCChooseBestParameters::~QMCChooseBestParameters()
    {                                      
      delete vmcEngine;                     
    }                                      
                                                    
    
    bool QMCChooseBestParameters::run()
    {
      Timer t1;                              
      app_log() << "<optimization-report>" << endl;                                                                       
        
        if (W.getActiveWalkers()>NumOfVMCWalkers)
        {                                      
          W.destroyWalkers(W.getActiveWalkers()-NumOfVMCWalkers);
          app_log() << "  QMCChooseBestParameters::generateSamples removed walkers." << endl;
          app_log() << "  Number of Walkers per node " << W.getActiveWalkers() << endl;
        }                                                                              
        
        vmcEngine->QMCDriverMode.set(QMC_OPTIMIZE,1);
        vmcEngine->QMCDriverMode.set(QMC_WARMUP,0);  
        
        vmcEngine->setValue("current",0);//reset CurrentStep 
        app_log() << "<vmc stage=\"main\" blocks=\"" << nBlocks << "\">" << endl;
        t1.restart();           
                                                             
        branchEngine->flush(0);                                               
        branchEngine->reset();                                                
        vmcEngine->run();                                                        
        app_log() << "  Execution time = " << t1.elapsed() << endl;              
        app_log() << "</vmc>" << endl;                                           
        
        
        opt_variables_type OptVariablesForPsi;
        OptVariablesForPsi.clear();
        WF->checkInVariables(OptVariablesForPsi);
        
        //write parameter history and energies to the parameter file in the trial wave function through opttarget
        RealType e,w,var;
        vmcEngine->Estimators->getEnergyAndWeight(e,w,var);
        WF->coefficientHistory.addParams(OptVariablesForPsi,e,var);
        
        app_log()<<"Using alpha: "<<alpha<<endl;
        app_log()<<"  energy:"<<1.0-alpha<<"  variance:"<<alpha<<endl;
        //choose best set of parameters
        opt_variables_type bestCoeffs = WF->coefficientHistory.getBestCoefficients(1.0-alpha,alpha,!myComm->rank());
        
        //check back into the WF
        WF->resetParameters(bestCoeffs);
        for (int i=0; i<psiClones.size(); ++i)
           psiClones[i]->resetParameters(bestCoeffs);
           
        //for (int i=0; i<psiClones.size(); ++i) 
        //{
        //  app_log()<<i<<endl;
        //  psiClones[i]->reportStatus(app_log());
        //}
        
        if (!myComm->rank())
        {
          app_log()<<"Best Parameters are:"<<endl;
          bestCoeffs.print(app_log()); 
          //app_log()<<"WF params"<<endl;
          //WF->reportStatus(app_log());
        }
      return true;
    }                                                                       
    
    /** Parses the xml input file for parameter definitions for the wavefunction optimization.
    * @param q current xmlNode                                                               
    * @return true if successful                                                             
    */                                                                                       
    bool                                                                                      
    QMCChooseBestParameters::put(xmlNodePtr q)                                                      
    {                                                                                         
      string useGPU("no");
      string vmcMove("pbyp");
      OhmmsAttributeSet oAttrib;
      oAttrib.add(useGPU,"gpu");
      oAttrib.add(vmcMove,"move");
      oAttrib.put(q);             
      
      xmlNodePtr qsave=q;
      xmlNodePtr cur=qsave->children;
      
      ParameterSet pAttrib;
      pAttrib.add(alpha,"alpha","scalar");
      pAttrib.put(q);
      
      int pid=OHMMS::Controller->rank();
      while (cur != NULL)               
      {                               
        string cname((const char*)(cur->name));
        if (cname == "mcwalkerset")            
        {                                    
          mcwalkerNodePtr.push_back(cur);    
        }                                                                                         
        cur=cur->next;
      }
      //no walkers exist, add 10
      if (W.getActiveWalkers() == 0) addWalkers(omp_get_max_threads());
      
      NumOfVMCWalkers=W.getActiveWalkers();
      
      //create VMC engine
      if(vmcEngine ==0) {
#if defined (QMC_CUDA)
        if (useGPU == "yes")
	  vmcEngine = new VMCcuda(W,Psi,H,psiPool);
        else
#endif
          vmcEngine = new VMCSingleOMP(W,Psi,H,hamPool,psiPool);

        vmcEngine->setUpdateMode(vmcMove[0] == 'p');
        vmcEngine->initCommunicator(myComm);
      }
      vmcEngine->setStatus(RootName,h5FileRoot,AppendRun);
      vmcEngine->process(qsave);
      
      bool success=true;
      return success;
    }
    
}
/***************************************************************************
* $RCSfile$   $Author: jnkim $
* $Revision: 1286 $   $Date: 2006-08-17 12:33:18 -0500 (Thu, 17 Aug 2006) $
* $Id: QMCChooseBestParameters.cpp 1286 2006-08-17 17:33:18Z jnkim $
***************************************************************************/

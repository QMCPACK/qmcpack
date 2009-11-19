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
#include "QMCDrivers/QMCLinearOptimize.h"                                             
#include "Particle/HDFWalkerIO.h"                                                     
#include "Particle/DistanceTable.h"                                                   
#include "OhmmsData/AttributeSet.h"                                                   
#include "Message/CommOperators.h"                                                    
#if defined(ENABLE_OPENMP)                                                            
#include "QMCDrivers/VMC/VMCSingleOMP.h"                                              
#include "QMCDrivers/QMCCostFunctionOMP.h"                                            
#endif                                                                                
#include "QMCDrivers/VMC/VMCSingle.h"                                                 
#include "QMCDrivers/QMCCostFunctionSingle.h"                                         
#include "QMCApp/HamiltonianPool.h"                                                   
#include "Numerics/Blasf.h"                                                           
#include <cassert>                                                                    
namespace qmcplusplus                                                                 
{
  
  QMCLinearOptimize::QMCLinearOptimize(MCWalkerConfiguration& w,
    TrialWaveFunction& psi, QMCHamiltonian& h, HamiltonianPool& hpool): QMCDriver(w,psi,h),
    PartID(0), NumParts(1), WarmupBlocks(10),                                                                               
    SkipSampleGeneration("no"), hamPool(hpool),                                                                             
    optTarget(0), vmcEngine(0), Max_iterations(1),                                                                          
    wfNode(NULL), optNode(NULL), exp0(-8), allowedCostDifference(2.0e-6),
    nstabilizers(3), stabilizerScale(4.0), bigChange(50), eigCG(1)
    {                                                                                                                           
      //set the optimization flag                                                                                               
      QMCDriverMode.set(QMC_OPTIMIZE,1);                                                                                        
      //read to use vmc output (just in case)                                                                                   
      RootName = "pot";                                                                                                         
      QMCType ="QMCLinearOptimize";                                                                                             
      optmethod = "Linear";                                                                                                     
      m_param.add(WarmupBlocks,"warmupBlocks","int");                                                                           
      m_param.add(SkipSampleGeneration,"skipVMC","string");                                                                     
      m_param.add(Max_iterations,"max_its","int"); 
      
      m_param.add(exp0,"exp0","int");                                                                        
      m_param.add(nstabilizers,"nstabilizers","int");                                                                        
      m_param.add(stabilizerScale,"stabilizerscale","double");
      m_param.add(allowedCostDifference,"alloweddifference","double"); 
      m_param.add(bigChange,"bigchange","double"); 
      m_param.add(eigCG,"eigcg","int");
      //Set parameters for line minimization:
      
    }
    
    /** Clean up the vector */
    QMCLinearOptimize::~QMCLinearOptimize()
    {                                      
      delete vmcEngine;                    
      delete optTarget;                    
    }                                      
    
    
    QMCLinearOptimize::RealType QMCLinearOptimize::Func(RealType dl)
    {
      for (int i=0;i<optparm.size();i++) optTarget->Params(i) = optparm[i] + dl*optdir[i];
      return optTarget->Cost(false);
      // return 0;                                                                            
    }                                                                                     
    
    /** Add configuration files for the optimization
    * @param a root of a hdf5 configuration file   
    */                                             
    void QMCLinearOptimize::addConfiguration(const string& a)
    {                                                        
      if (a.size()) ConfigFile.push_back(a);                 
    }                                                        
    
    bool QMCLinearOptimize::run()
    {                            
      optTarget->initCommunicator(myComm);
      //close files automatically generated by QMCDriver
      //branchEngine->finalize();                       
      
      
      //generate samples
      generateSamples();
      
      app_log() << "<opt stage=\"setup\">" << endl;
      app_log() << "  <log>"<<endl;                
      
      //reset the rootname
      optTarget->setRootName(RootName);
      optTarget->setWaveFunctionNode(wfNode);
      
      app_log() << "   Reading configurations from h5FileRoot " << endl;
      //get configuration from the previous run                         
      Timer t1;                                                         
      
      optTarget->getConfigurations(h5FileRoot);
      optTarget->checkConfigurations();        
      
      app_log() << "  Execution time = " << t1.elapsed() << endl;
      app_log() << "  </log>"<<endl;                             
      app_log() << "</opt>" << endl;                             
      
      app_log() << "<opt stage=\"main\" walkers=\""<< optTarget->getNumSamples() << "\">" << endl;
      app_log() << "  <log>" << endl;                                                             
      
      optTarget->setTargetEnergy(branchEngine->getEref());
      
      t1.restart();
      
      ///Here is our optimization routine
      bool Valid(true);                  
      int Total_iterations(0);           
      
      TOL = allowedCostDifference;         

      int N=optTarget->NumParams() + 1;
      int numParams = optTarget->NumParams();
      optdir.resize(numParams,0);

      vector<RealType> currentParameterDirections(N,0);
      vector<RealType> currentParameters(numParams,0);
      for (int i=0;i<numParams; i++) currentParameters[i] = optTarget->Params(i);
      while (Max_iterations>Total_iterations)
      {
        Total_iterations+=1;               
        app_log()<<"Iteration: "<<Total_iterations<<"/"<<Max_iterations<<endl;
        
        
//         store this for use in later tries
        int bestStability(0);
        for(int tries=0;tries<eigCG;tries++){
          RealType lastCost(optTarget->Cost(false));
          RealType newCost(lastCost);
          
          Matrix<RealType> Ham(N,N);
          Matrix<RealType> Ham2(N,N);
          Matrix<RealType> S(N,N);
          Matrix<RealType> S2(N,N);
        
          for (int i=0;i<numParams; i++) optTarget->Params(i) = currentParameters[i];
          optTarget->fillOverlapHamiltonianMatrix(Ham, S);
          optTarget->fillOverlapHamiltonianSquaredMatrix(Ham2, S2);
          
          vector<RealType> bestParameters(currentParameters);
          
          for(int stability=0;stability<((tries==0)?nstabilizers:1);stability++){
            Matrix<RealType> HamT(N,N), ST(N,N);
            Matrix<RealType> HamT2(N,N), ST2(N,N);
            for (int i=0;i<N;i++)
              for (int j=0;j<N;j++)
              {
                HamT(i,j)= (Ham)(j,i);        
                ST(i,j)= (S)(j,i);
                HamT2(i,j)= (Ham2)(j,i);        
                ST2(i,j)= (S2)(j,i);
              }
              RealType Xs(0);
              if (tries==0) Xs = std::pow(10.0,exp0 + stabilizerScale*stability);
              else Xs = std::pow(10.0,exp0 + stabilizerScale*bestStability);
              for (int i=1;i<N;i++) HamT(i,i) += Xs;
              for (int i=1;i<N;i++) HamT2(i,i) += Xs;
              
              
              char jl('N');
              char jr('V');
              vector<RealType> alphar(N),alphai(N),beta(N);
              Matrix<RealType> eigenT(N,N);                
              int info;                                    
              int lwork(-1);                               
              vector<RealType> work(1);                    
              
              RealType tt(0);
              int t(1);
              dggev(&jl, &jr, &N, HamT.data(), &N, ST.data(), &N, &alphar[0], &alphai[0], &beta[0],&tt,&t, eigenT.data(), &N, &work[0], &lwork, &info);
              lwork=work[0];
              work.resize(lwork);
              
              //here is the pure energy part
              if (false)
              {
                ///RealType==double to use this one, ned to write case where RealType==float                                                             
                dggev(&jl, &jr, &N, HamT.data(), &N, ST.data(), &N, &alphar[0], &alphai[0], &beta[0],&tt,&t, eigenT.data(), &N, &work[0], &lwork, &info);
                assert(info==0);
                
                vector<std::pair<RealType,int> > mappedEigenvalues(numParams);
                for (int i=0;i<numParams;i++)
                {
                  mappedEigenvalues[i].first=alphar[i]/beta[i];
                  mappedEigenvalues[i].second=i;
                }
                std::sort(mappedEigenvalues.begin(),mappedEigenvalues.end());
                for (int i=0;i<N;i++) currentParameterDirections[i] = eigenT(mappedEigenvalues[tries].second,i)/eigenT(mappedEigenvalues[tries].second,0);
              }

              
              //here is the mixed variance and energy part
              if (true)
              {
                // CG like algorithm, we try to move in maxtries directions. eigenvalues are orthogonal.
                vector<RealType> alphar2(N),alphai2(N),beta2(N);
                Matrix<RealType> eigenT2(N,N);              
                dggev(&jl, &jr, &N, HamT2.data(), &N, ST2.data(), &N, &alphar2[0], &alphai2[0], &beta2[0],&tt,&t, eigenT2.data(), &N, &work[0], &lwork, &info);
                assert(info==0);
                
                vector<std::pair<RealType,int> > mappedEigenvalues2(numParams);
                for (int i=0;i<numParams;i++)
                {
                  mappedEigenvalues2[i].first=alphar2[i]/beta2[i];
                  mappedEigenvalues2[i].second=i;
                }
                std::sort(mappedEigenvalues2.begin(),mappedEigenvalues2.end());
                for (int i=0;i<N;i++) currentParameterDirections[i] = eigenT2( mappedEigenvalues2[tries].second,i)/eigenT2(mappedEigenvalues2[tries].second,0);
              }
              
              
//               if(mappedEigenvalues[0].first<0)
//               {
//                 app_log()<<mappedEigenvalues[0].first<<endl;
//                 for (int i=0;i<N;i++) HamT(i,i) -= mappedEigenvalues[0].first;
//                 dggev(&jl, &jr, &N, HamT.data(), &N, ST.data(), &N, &alphar[0], &alphai[0], &beta[0],&tt,&t, eigenT.data(), &N, &work[0], &lwork, &info);
//                 assert(info==0);
//                 for (int i=0;i<numParams;i++)
//                 {
//                   mappedEigenvalues[i].first=alphar[i]/beta[i];
//                   mappedEigenvalues[i].second=i;
//                 }
//                 std::sort(mappedEigenvalues.begin(),mappedEigenvalues.end());
//               }

              
              
//               if (false)
//               {
//                 //Umrigar and Sorella suggest using 0.5 for xi.                                        
//                 RealType xi=0.5;                                                                       
//                 RealType D(1.0);                                                                       
//                 for (int i=0;i<numParams;i++)                                                                
//                 {                                                                                    
//                   if (optTarget->getType(i) != 2)                                                    
//                   {                                                                                
//                     for (int j=0;j<numParams;j++)                                                        
//                     {                                                                            
//                       if (optTarget->getType(j) != 2) D += S(j+1,i+1)*currentParameterDirections[i+1]*currentParameterDirections[j+1];           
//                     }                                                                            
//                   }                                                                                
//                 }                                                                                    
//                 D = std::sqrt(std::abs(D));                                                            
//                 
//                 vector<RealType> N_i(numParams,0);
//                 for (int i=0;i<numParams;i++)     
//                 {                         
//                   RealType tsumN(0);      
//                   for (int j=0;j<numParams;j++) 
//                   {                     
//                     if (optTarget->getType(j) != 2)
//                     {                            
//                       tsumN += S(i+1,j+1)*currentParameterDirections[j+1];
//                     }                             
//                   }                                 
//                   N_i[i] += (1-xi)*tsumN  / (xi*D + (1-xi));
//                 }                                           
//                 
//                 RealType rescale(1);
//                 for (int j=0;j<numParams;j++) rescale -= N_i[j]*currentParameterDirections[j+1] ;
//                 rescale = 1.0/rescale;                             
//                 if ((rescale==rescale)&&(rescale!=0))              
//                 {                                                
//                   for (int i=0;i<numParams; i++)                     
//                   {                                            
//                     if (optTarget->getType(i) != 2) currentParameterDirections[i+1]*=rescale;
//                   }                                                  
//                 }                                                      
//                 
//               }

 //If not rescaling and linear parameters, step size and grad are the same.
              LambdaMax = 1.0;
              optparm= currentParameters;
              for (int i=0;i<numParams; i++) optdir[i] = currentParameterDirections[i+1];
//               if (tries==0) quadstep=0.2;
//               else quadstep=0.01;
              lineoptimization2();
              
              if (Lambda==Lambda)        
              {
                for (int i=0;i<numParams; i++) optTarget->Params(i) = optparm[i] + Lambda * optdir[i];
                newCost = optTarget->Cost(false);
                app_log()<<" OldCost: "<<lastCost<<" NewCost: "<<newCost<<endl;
//                 quit if newcost is greater than lastcost. E(Xs) looks quadratic (between steepest descent and parabolic)
                
                if ((newCost > (lastCost - bigChange))&&(newCost < lastCost)&&(newCost==newCost))
                {
                  //Move was acceptable 
                  for (int i=0;i<numParams; i++) bestParameters[i] = optTarget->Params(i);
                  bestStability=stability; lastCost=newCost;
                }
                else if (newCost>lastCost+0.001) stability = nstabilizers;
              }
            }
            // Let us try a single steepest descent step just in case
//         optparm= currentParameters;
//         LambdaMax = 0.02;
//         optTarget->GradCost(optdir, optparm, 0);
//         lineoptimization2();
//         if (Lambda==Lambda)        
//         {
//           for (int i=0;i<numParams; i++) optTarget->Params(i) = optparm[i] + Lambda * optdir[i];
//           newCost = optTarget->Cost(false);
//           app_log()<<" OldCost: "<<lastCost<<" NewCost: "<<newCost<<endl;
//           if ((newCost > lastCost - bigChange)&&(newCost < lastCost)&&(newCost==newCost))
//           {
//             //Move was acceptable 
//             for (int i=0;i<numParams; i++) bestParameters[i] = optTarget->Params(i);
//             lastCost=newCost;
//           }
//         }
//         
        for (int i=0;i<numParams; i++) optTarget->Params(i) = bestParameters[i]; 
        currentParameters=bestParameters;
        }
      }
      
      MyCounter++;
      app_log() << "  Execution time = " << t1.elapsed() << endl;
      app_log() << "  </log>" << endl;                           
      optTarget->reportParameters();                             
      app_log() << "</opt>" << endl;                             
      app_log() << "</optimization-report>" << endl;             
      return (optTarget->getReportCounter() > 0);                
    }                                                            
    
    void QMCLinearOptimize::generateSamples()
    {                                        
      Timer t1;                              
      app_log() << "<optimization-report>" << endl;
      //if(WarmupBlocks)                           
      //{                                          
        //  app_log() << "<vmc stage=\"warm-up\" blocks=\"" << WarmupBlocks << "\">" << endl;
        //  //turn off QMC_OPTIMIZE                                                          
        //  vmcEngine->setValue("blocks",WarmupBlocks);                                      
        //  vmcEngine->QMCDriverMode.set(QMC_WARMUP,1);                                      
        //  vmcEngine->run();                                                                
        //  vmcEngine->setValue("blocks",nBlocks);                                           
        //  app_log() << "  Execution time = " << t1.elapsed() << endl;                      
        //  app_log() << "</vmc>" << endl;                                                   
        //}                                                                                  
        
        if (W.getActiveWalkers()>NumOfVMCWalkers)
        {                                      
          W.destroyWalkers(W.getActiveWalkers()-NumOfVMCWalkers);
          app_log() << "  QMCLinearOptimize::generateSamples removed walkers." << endl;
          app_log() << "  Number of Walkers per node " << W.getActiveWalkers() << endl;
        }                                                                              
        
        vmcEngine->QMCDriverMode.set(QMC_OPTIMIZE,1);
        vmcEngine->QMCDriverMode.set(QMC_WARMUP,0);  
        
        //vmcEngine->setValue("recordWalkers",1);//set record
        vmcEngine->setValue("current",0);//reset CurrentStep 
        app_log() << "<vmc stage=\"main\" blocks=\"" << nBlocks << "\">" << endl;
        t1.restart();                                                            
        //     W.reset();                                                            
        //     branchEngine->flush(0);                                               
        //     branchEngine->reset();                                                
        vmcEngine->run();                                                        
        app_log() << "  Execution time = " << t1.elapsed() << endl;              
        app_log() << "</vmc>" << endl;                                           
        
        //branchEngine->Eref=vmcEngine->getBranchEngine()->Eref;
        branchEngine->setTrialEnergy(vmcEngine->getBranchEngine()->getEref());
        //set the h5File to the current RootName                              
        h5FileRoot=RootName;                                                  
    }                                                                       
    
    /** Parses the xml input file for parameter definitions for the wavefunction optimization.
    * @param q current xmlNode                                                               
    * @return true if successful                                                             
    */                                                                                       
    bool                                                                                      
    QMCLinearOptimize::put(xmlNodePtr q)                                                      
    {                                                                                         
      
      string vmcMove("pbyp");
      OhmmsAttributeSet oAttrib;
      oAttrib.add(vmcMove,"move");
      oAttrib.put(q);             
      
      xmlNodePtr qsave=q;
      xmlNodePtr cur=qsave->children;
      
      
      int pid=OHMMS::Controller->rank();
      while (cur != NULL)               
      {                               
        string cname((const char*)(cur->name));
        if (cname == "mcwalkerset")            
        {                                    
          mcwalkerNodePtr.push_back(cur);    
        }                                    
        else if (cname == "optimizer")         
        {                                    
          xmlChar* att= xmlGetProp(cur,(const xmlChar*)"method");
          if (att)                                               
          {                                                    
            optmethod = (const char*)att;                      
          }                                                    
          optNode=cur;                                           
        }                                                        
        else if (cname == "optimize")                              
        {                                                        
          xmlChar* att= xmlGetProp(cur,(const xmlChar*)"method");
          if (att)                                               
          {                                                    
            optmethod = (const char*)att;                      
          }                                                    
        }                                                        
        cur=cur->next;
      }
      //no walkers exist, add 10
      if (W.getActiveWalkers() == 0) addWalkers(omp_get_max_threads());
      
      NumOfVMCWalkers=W.getActiveWalkers();
      
      //create VMC engine
      if (vmcEngine ==0)
      {
        #if defined(ENABLE_OPENMP)
        if (omp_get_max_threads()>1)
          //             vmcEngine = new DMCOMP(W,Psi,H,hamPool);
        vmcEngine = new VMCSingleOMP(W,Psi,H,hamPool);
        else
          #endif
          vmcEngine = new VMCSingle(W,Psi,H);
        vmcEngine->setUpdateMode(vmcMove[0] == 'p');
        vmcEngine->initCommunicator(myComm);
      }
      
      vmcEngine->setStatus(RootName,h5FileRoot,AppendRun);
      vmcEngine->process(qsave);
      
      bool success=true;
      if (optTarget == 0)
      {
        #if defined(ENABLE_OPENMP)
        if (omp_get_max_threads()>1)
        {
          optTarget = new QMCCostFunctionOMP(W,Psi,H,hamPool);
        }
        else
          #endif
          optTarget = new QMCCostFunctionSingle(W,Psi,H);
        optTarget->setStream(&app_log());
        success=optTarget->put(q);
      }
      return success;
    }
}
/***************************************************************************
* $RCSfile$   $Author: jnkim $
* $Revision: 1286 $   $Date: 2006-08-17 12:33:18 -0500 (Thu, 17 Aug 2006) $
* $Id: QMCLinearOptimize.cpp 1286 2006-08-17 17:33:18Z jnkim $
***************************************************************************/

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
#if defined(QMC_CUDA)
  #include "QMCDrivers/VMC/VMC_CUDA.h"
  #include "QMCDrivers/QMCCostFunctionCUDA.h"
#endif

                                                                 
namespace qmcplusplus                                                                 
{
  
  QMCLinearOptimize::QMCLinearOptimize(MCWalkerConfiguration& w,
    TrialWaveFunction& psi, QMCHamiltonian& h, HamiltonianPool& hpool): QMCDriver(w,psi,h),
    PartID(0), NumParts(1), WarmupBlocks(10), 
    SkipSampleGeneration("no"), hamPool(hpool),
    optTarget(0), vmcEngine(0), Max_iterations(1),
    wfNode(NULL), optNode(NULL), exp0(-8), allowedCostDifference(2.0e-6), 
    nstabilizers(3), stabilizerScale(4.0), bigChange(1), eigCG(1), w_beta(1),
    UseQuarticMin("yes")
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
      m_param.add(w_beta,"beta","double");
      quadstep=3.0;
      m_param.add(quadstep,"stepsize","double");
      m_param.add(UseQuarticMin,"UseQuarticMin","string");
      m_param.add(LambdaMax,"LambdaMax","double");
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
      validFuncVal= optTarget->IsValid;
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
 
// mmorales
        if(!Valid) {
          app_log() <<"Aborting current opt cycle due to small wfn overlap during correlated sampling. If this happens to frequently, try reducing the step size of the line minimization or reduce the number of cycles. " <<endl; 
          break;
        }       
 
        
//         store this for use in later tries
        int bestStability(0);
        vector<vector<RealType> > LastDirections;
        RealType deltaPrms(-1.0);
        for(int tries=0;tries<eigCG;tries++){
          Matrix<RealType> Ham(N,N);
          Matrix<RealType> Ham2(N,N);
          Matrix<RealType> S(N,N);
          vector<RealType> BestDirection(N,0);

          for (int i=0;i<numParams; i++) optTarget->Params(i) = currentParameters[i];

          RealType lastCost(optTarget->Cost(true));
          RealType newCost(lastCost);

// mmorales
          if(!optTarget->IsValid) 
          {
            Valid=false;
            break;
          }
          
// correlated sampling was called in Cost, so the call here is redundant
          optTarget->fillOverlapHamiltonianMatrices(Ham2, Ham, S);

// mmorales
          if(!optTarget->IsValid) 
          {
            Valid=false;
            break;
          }
          
          vector<RealType> bestParameters(currentParameters);
          bool acceptedOneMove(false);
          for(int stability=0;stability<((tries==0)?nstabilizers:1);stability++){

            Matrix<RealType> HamT(N,N), ST(N,N), HamT2(N,N);
            for (int i=0;i<N;i++)
              for (int j=0;j<N;j++)
              {
                HamT(i,j)= (Ham)(j,i);        
                ST(i,j)= (S)(j,i);
                HamT2(i,j)= (Ham2)(j,i);        
              }
              RealType Xs(0);
              if (tries==0) Xs = std::pow(10.0,exp0 + stabilizerScale*stability);
              else Xs = std::pow(10.0,exp0 + stabilizerScale*bestStability);
              for (int i=1;i<N;i++) HamT(i,i) += Xs;

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

              Matrix<RealType> ST2(N,N);
              RealType H2rescale=1.0/std::sqrt(std::abs(HamT2(0,0)));
              for (int i=0;i<N;i++)  for (int j=0;j<N;j++) HamT2(i,j) *= H2rescale;
              for (int i=0;i<N;i++)  for (int j=0;j<N;j++) ST2(i,j) = (1.0-w_beta)*ST(i,j) + w_beta*HamT2(i,j);

              dggev(&jl, &jr, &N, HamT.data(), &N, ST2.data(), &N, &alphar[0], &alphai[0], &beta[0],&tt,&t, eigenT.data(), &N, &work[0], &lwork, &info);
              assert(info==0);

                vector<std::pair<RealType,int> > mappedEigenvalues(N);
                for (int i=0;i<N;i++)
                {
                  RealType evi(alphar[i]/beta[i]);
                  mappedEigenvalues[i].first=(evi-Ham(0,0))*(evi-Ham(0,0));
                  mappedEigenvalues[i].second=i;
//                  app_log()<<evi<<" "<<(evi-Ham(0,0))*(evi-Ham(0,0))<<endl;
                }
                app_log()<<endl;
                std::sort(mappedEigenvalues.begin(),mappedEigenvalues.end());

                for (int i=0;i<N;i++) currentParameterDirections[i] = eigenT( mappedEigenvalues[0].second,i)/eigenT(mappedEigenvalues[0].second,0);
                //eigenCG part
                RealType nrmold(0);

                for(int ldi=0;ldi<LastDirections.size();ldi++)
                {
                  for (int i=1;i<N;i++) nrmold += LastDirections[ldi][i]*LastDirections[ldi][i];
                  RealType ovlpold(0), nrmnew(0);
                  for (int i=1;i<N;i++) nrmnew +=currentParameterDirections[i]*currentParameterDirections[i];
                  for (int i=1;i<N;i++) ovlpold += LastDirections[ldi][i]*currentParameterDirections[i];
                  for (int i=1;i<N;i++) currentParameterDirections[i] -= ovlpold/nrmold * LastDirections[ldi][i];
                  nrmnew=0;
                  for (int i=1;i<N;i++) nrmnew +=currentParameterDirections[i]*currentParameterDirections[i];
                  //rescale to be something that worked before
                  for (int i=1;i<N;i++) currentParameterDirections[i] *= nrmold/nrmnew;
                }
 //If not rescaling and linear parameters, step size and grad are the same.
//              LambdaMax=1.0;
              optparm= currentParameters;
              for (int i=0;i<numParams; i++) optdir[i] = currentParameterDirections[i+1];
              
              RealType dopt(0);
              for (int i=0;i<numParams; i++) dopt += optdir[i] * optdir[i];
              dopt =std::sqrt(dopt)/RealType(numParams);
              
              largeQuarticStep=1e3;
              if (deltaPrms>0) quadstep=deltaPrms/dopt;
              if(UseQuarticMin=="yes") {
                Valid=lineoptimization();
                if (dopt*std::abs(Lambda)>bigChange)
                  Valid=lineoptimization2();
              } else {
                Valid=lineoptimization2();
              }
              if(!Valid) break;

              dopt *= std::abs(Lambda);
                
              if ( (Lambda==Lambda)&&(dopt<bigChange))
              {
                app_log() <<endl <<"lambda: " <<Lambda <<endl;
                for (int i=0;i<numParams; i++) {
                  app_log() <<"param: i, currValue, optdir: " 
                            <<i <<"  " <<optparm[i] <<"  " <<optdir[i] <<endl;
                }


                for (int i=0;i<numParams; i++) optTarget->Params(i) = optparm[i] + Lambda * optdir[i];
                newCost = optTarget->Cost(false);

// mmorales
                if(!optTarget->IsValid)
                {
                  Valid=false;
                  break;
                }
                app_log()<<" OldCost: "<<lastCost<<" NewCost: "<<newCost<<" RMS step size:"<<dopt<<endl;
//                 quit if newcost is greater than lastcost. E(Xs) looks quadratic (between steepest descent and parabolic)
                
                if ((newCost < lastCost)&&(newCost==newCost))
                {
                  //Move was acceptable 
                  for (int i=0;i<numParams; i++) bestParameters[i] = optTarget->Params(i);
                  bestStability=stability; lastCost=newCost;
                  BestDirection=currentParameterDirections;
                  acceptedOneMove=true;
                  
                  deltaPrms=dopt;
                }
//                else if (newCost>lastCost+0.001) stability = nstabilizers;
              }
              else
              {
                app_log()<<"  Failed Step. RMS step Size:"<<dopt<<endl;
              }
            }
          if(acceptedOneMove && Valid)
          {
            for (int i=0;i<numParams; i++) optTarget->Params(i) = bestParameters[i]; 
            currentParameters=bestParameters;
            LastDirections.push_back(BestDirection);
//             app_log()<< " Wave Function Parameters updated."<<endl;
//             optTarget->reportParameters();
          }
          else
          {
            for (int i=0;i<numParams; i++) optTarget->Params(i) = currentParameters[i];   
            tries=eigCG;
          }
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
      string useGPU("no");
      string vmcMove("pbyp");
      OhmmsAttributeSet oAttrib;
      oAttrib.add(useGPU,"gpu");
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
      if(vmcEngine ==0) {
#if defined (QMC_CUDA)
        if (useGPU == "yes")
	  vmcEngine = new VMCcuda(W,Psi,H);
        else
#endif
//#if defined(ENABLE_OPENMP)
//        if(omp_get_max_threads()>1)
//          vmcEngine = new VMCSingleOMP(W,Psi,H,hamPool);
//        else
//#endif
//          vmcEngine = new VMCSingle(W,Psi,H);
        vmcEngine = new VMCSingleOMP(W,Psi,H,hamPool);
        vmcEngine->setUpdateMode(vmcMove[0] == 'p');
        vmcEngine->initCommunicator(myComm);
      }
      vmcEngine->setStatus(RootName,h5FileRoot,AppendRun);
      vmcEngine->process(qsave);
      
      bool success=true;

      if(optTarget == 0) {
#if defined (QMC_CUDA)
        if (useGPU == "yes") 
  	  optTarget = new QMCCostFunctionCUDA(W,Psi,H,hamPool);
        else
#endif
#if defined(ENABLE_OPENMP)
        if(omp_get_max_threads()>1) {
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

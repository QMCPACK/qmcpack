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
#include "QMCDrivers/QMCFixedSampleLinearOptimize.h"
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
#include "Numerics/MatrixOperators.h"
#include <cassert>
#if defined(QMC_CUDA)
#include "QMCDrivers/VMC/VMC_CUDA.h"
#include "QMCDrivers/QMCCostFunctionCUDA.h"
#endif
#include "Numerics/LinearFit.h"
#include <iostream>
#include <fstream>

/*#include "Message/Communicate.h"*/

namespace qmcplusplus
{

QMCFixedSampleLinearOptimize::QMCFixedSampleLinearOptimize(MCWalkerConfiguration& w,
                                     TrialWaveFunction& psi, QMCHamiltonian& h, HamiltonianPool& hpool, WaveFunctionPool& ppool): QMCLinearOptimize(w,psi,h,hpool,ppool),
        Max_iterations(1), exp0(-16), exp1(0),  nstabilizers(10), stabilizerScale(0.5), bigChange(1), eigCG(1), TotalCGSteps(2), w_beta(0.0),
        MinMethod("quartic"), GEVtype("mixed"), StabilizerMethod("best"), GEVSplit("no")
{
    //set the optimization flag
    QMCDriverMode.set(QMC_OPTIMIZE,1);
    //read to use vmc output (just in case)
    RootName = "pot";
    QMCType ="QMCFixedSampleLinearOptimize";
    m_param.add(WarmupBlocks,"warmupBlocks","int");
    m_param.add(Max_iterations,"max_its","int");
    m_param.add(nstabilizers,"nstabilizers","int");
    m_param.add(stabilizerScale,"stabilizerscale","double");
    m_param.add(bigChange,"bigchange","double");
    m_param.add(eigCG,"eigcg","int");
    m_param.add(TotalCGSteps,"cgsteps","int");
    m_param.add(w_beta,"beta","double");
    quadstep=-1.0;
    m_param.add(quadstep,"stepsize","double");
    m_param.add(exp0,"exp0","double");
    m_param.add(exp1,"exp1","double");
    m_param.add(MinMethod,"MinMethod","string");
    m_param.add(GEVtype,"GEVMethod","string");
    m_param.add(GEVSplit,"GEVSplit","string");
    m_param.add(StabilizerMethod,"StabilizerMethod","string");
    m_param.add(LambdaMax,"LambdaMax","double");
    //Set parameters for line minimization:
    this->add_timers(myTimers);
}

/** Clean up the vector */
QMCFixedSampleLinearOptimize::~QMCFixedSampleLinearOptimize()
{
}

QMCFixedSampleLinearOptimize::RealType QMCFixedSampleLinearOptimize::Func(RealType dl)
{
    for (int i=0; i<optparm.size(); i++) optTarget->Params(i) = optparm[i] + dl*optdir[i];
    QMCLinearOptimize::RealType c = optTarget->Cost(false);
    //only allow this to go false if it was true. If false, stay false
    if (validFuncVal) validFuncVal = optTarget->IsValid;
    return c;
}

bool QMCFixedSampleLinearOptimize::run()
{
    start();
    bool Valid(true);
    int Total_iterations(0);
    savedQuadstep=quadstep;

//size of matrix
    numParams = optTarget->NumParams();
    N = numParams + 1;

//  solve CSFs and other parameters separately then rescale elements accordingly
    int first,last;
    getNonLinearRange(first,last);
//  There is only one type of parameter and all are non-linear, don't split it up
    if (last-first==numParams) GEVSplit=="no";

//     initialize our parameters
    vector<RealType> currentParameterDirections(N,0);
    vector<RealType> currentParameters(numParams,0);
    optdir.resize(numParams,0);
    optparm.resize(numParams,0);
    for (int i=0; i<numParams; i++) currentParameters[i] = optTarget->Params(i);

    Matrix<RealType> Ham(N,N);
    Matrix<RealType> Ham2(N,N);
    Matrix<RealType> Var(N,N);
    Matrix<RealType> S(N,N);
    vector<RealType> BestDirection(N,0);
    vector<RealType> bestParameters(currentParameters);
    vector<RealType> GEVSplitParameters(numParams,0);


    while (Total_iterations < Max_iterations)
    {
        Total_iterations+=1;
        app_log()<<"Iteration: "<<Total_iterations<<"/"<<Max_iterations<<endl;

// mmorales
        if (!ValidCostFunction(Valid)) continue;

//this is the small amount added to the diagonal to stabilize the eigenvalue equation. 10^stabilityBase
        RealType stabilityBase(exp0);
        //This is the amount we add to the linear parameters
        RealType linearStabilityBase(exp1);

        vector<vector<RealType> > LastDirections;
        RealType deltaPrms(-1.0);
        for (int tries=0; tries<TotalCGSteps; tries++)
        {
            bool acceptedOneMove(false);
            int tooManyTries(20);

            Matrix<RealType> Left_tmp(N,N);
            Matrix<RealType> Left(N,N);
            Matrix<RealType> Right(N,N);

            vector<std::pair<RealType,RealType> > mappedStabilizers;
            if (nstabilizers<5)
            {
                if (StabilizerMethod=="fit") app_log()<<" Need 5 stabilizers minimum for the fit"<<endl;
                StabilizerMethod="best";
            }

            for (int i=0; i<numParams; i++) optTarget->Params(i) = currentParameters[i];

            myTimers[4]->start();
            RealType lastCost(optTarget->Cost(true));
            myTimers[4]->start();

            // mmorales
            Valid=optTarget->IsValid;
            if (!ValidCostFunction(Valid)) continue;
            
            RealType newCost(lastCost);
            optTarget->fillOverlapHamiltonianMatrices(Ham2, Ham, Var, S);
            RealType H2rescale=1.0;

            if (GEVtype=="H2")
            {
                Left_tmp=Ham;
                H2rescale=1.0/Ham2(0,0);
                Right=(1-w_beta)*S + w_beta*H2rescale*Ham2;
            }
            else
            {
                Right=S;
                Left_tmp=(1.0-w_beta)*Ham + w_beta*Var;
            }
            
            bool apply_inverse(true);
            if(apply_inverse)
            {
              invert_matrix(Right,false);
              MatrixOperators MO;
              MO.product(Right,Left_tmp,Left);
            }
            else
              Left=Left_tmp;

            
            //Find largest off-diagonal element compared to diagonal element.
            //This gives us an idea how well conditioned it is and can be used to stabilize.
            RealType od_largest(0);
            for (int i=0; i<N; i++) for (int j=0; j<N; j++)
              od_largest=std::max( std::max(od_largest,std::abs(Left(i,j))-std::abs(Left(i,i))), std::abs(Left(i,j))-std::abs(Left(j,j)));
            if (od_largest>0) od_largest = std::log(od_largest);
            else od_largest=1.0;

            RealType safe = Left(0,0);
            for (int stability=0; stability<nstabilizers; stability++)
            {
                Matrix<RealType> LeftT(N,N);//, RightT(N,N);
                for (int i=0; i<N; i++)
                    for (int j=0; j<N; j++)
                    {
                        LeftT(i,j)= Left(j,i);
                        //RightT(i,j)= Right(j,i);
                    }


                RealType XS(0);
//                 if ((StabilizerMethod=="fit")&&(stability==nstabilizers-1))
//                 {
//                     //Quartic fit the stabilizers we have tried and try to choose the best we can
//                     int nms=mappedStabilizers.size();
//                     
//                     
//                     int cTerms(3);
//                     vector<RealType>  Y(nms), Coefs(cTerms);
//                     Matrix<RealType> X(nms,cTerms);
//                     for (int i=0; i<nms; i++) X(i,0)=1.0;
//                     for (int i=0; i<nms; i++) X(i,1)=mappedStabilizers[i].second;
//                     for (int i=0; i<nms; i++) X(i,2)=X(i,1)*X(i,1);
// //                     for (int i=0; i<nms; i++) X(i,3)=X(i,2)*X(i,1);
// //                     for (int i=0; i<nms; i++) X(i,4)=X(i,3)*X(i,1);
//                     for (int i=0; i<nms; i++) Y[i]=mappedStabilizers[i].first;
//                     LinearFit(Y,X,Coefs);
//                     //lowest we will allow is a little less than the bare base stabilizer
//                     RealType dltaBest=std::max(stabilityBase-0.1, QuarticMinimum(Coefs));
//                     XS = dltaBest;
//                     stability=nstabilizers;
//                 }
            int nms=mappedStabilizers.size();
            if ((StabilizerMethod=="fit")&&(stability==nstabilizers-1))
            {
              if (nms>=5)
              {//Quartic fit the stabilizers we have tried and try to choose the best we can
                vector<RealType>  Y(nms), Coefs(5);
                Matrix<RealType> X(nms,5);
                for (int i=0; i<nms; i++) X(i,0)=1.0;
                for (int i=0; i<nms; i++) X(i,1)=mappedStabilizers[i].second;
                for (int i=0; i<nms; i++) X(i,2)=X(i,1)*X(i,1);
                for (int i=0; i<nms; i++) X(i,3)=X(i,2)*X(i,1);
                for (int i=0; i<nms; i++) X(i,4)=X(i,3)*X(i,1);
                for (int i=0; i<nms; i++) Y[i]=mappedStabilizers[i].first;
                LinearFit(Y,X,Coefs);
    //lowest we will allow is a little less than the bare base stabilizer
                RealType dltaBest=std::max(stabilityBase-0.1, QuarticMinimum(Coefs));
                XS = dltaBest;
              }
              else
              {//Quadratic fit the stabilizers we have tried and try to choose the best we can
                std::sort(mappedStabilizers.begin(),mappedStabilizers.end());
                vector<RealType>  Y(nms), Coefs(3);
                Matrix<RealType> X(nms,3);
                for (int i=0; i<nms; i++) X(i,0)=1.0;
                for (int i=0; i<nms; i++) X(i,1)=mappedStabilizers[i].second;
                for (int i=0; i<nms; i++) X(i,2)=X(i,1)*X(i,1);
                for (int i=0; i<nms; i++) Y[i]=mappedStabilizers[i].first;
                LinearFit(Y,X,Coefs);
                
                RealType quadraticMinimum(-1.0*Coefs[1]/Coefs[2]);
                RealType dltaBest=std::max(stabilityBase-0.1, quadraticMinimum);
    //               app_log()<<"smallest XS:      "<<X(0,1)<<endl;
    //               app_log()<<"quadraticMinimum: "<<quadraticMinimum<<endl;
                XS = dltaBest;
              }
               stability=nstabilizers;
            }


                RealType lowestEV(0);
//                 if ((GEVSplit=="rescale")||(GEVSplit=="freeze"))
//                 {
//                   //These are experimental and aren't very good.
//                     //dummy bool
//                     bool CSF_lower(true);
//                     lowestEV = getSplitEigenvectors(first,last,LeftT,RightT,currentParameterDirections,GEVSplitParameters,GEVSplit,CSF_lower);
//                 }
//                 else if (GEVSplit=="stability") //This seems to work pretty well.
//                 {
//                     if (XS==0)
//                     {
//                         XS  = std::exp(stabilityBase + stability*od_largest/nstabilizers);
//                         for (int i=first; i<last; i++) LeftT(i+1,i+1) += XS;
//                         
//                         RealType XS_lin = std::exp(linearStabilityBase + stabilityBase + stability*od_largest/nstabilizers);
//                         if (first==0) for (int i=last; i<N; i++) LeftT(i+1,i+1) += XS_lin;
//                         else for (int i=0; i<first; i++) LeftT(i+1,i+1) += XS_lin;
//                     }
//                     else //else XS is from the quartic fit
//                     {
//                       //Not sure how to control for the quartic fit and the two different stabilizers. This seems ok.
//                       //Better algorithm exists?
//                         for (int i=first; i<last; i++) LeftT(i+1,i+1) += std::exp(XS);
//                       
//                         RealType XS_lin = std::exp(linearStabilityBase+XS);
//                         if (first==0) for (int i=last; i<N; i++) LeftT(i+1,i+1) += XS_lin;
//                         for (int i=0; i<first; i++) LeftT(i+1,i+1) += XS_lin;
//                     }

//                     if (stability==0)
//                     {
//                          Only need to do this the first time we step into the routine
//                         bool CSF_lower(true);
//                         lowestEV=getSplitEigenvectors(first,last,LeftT,RightT,currentParameterDirections,GEVSplitParameters,GEVSplit,CSF_lower);
//                         if (tooLow(safe,lowestEV))
//                         {
//                             if (CSF_lower)
//                             {
//                                 linearStabilityBase+=stabilizerScale;
//                                 app_log()<<"Probably will not converge: CSF Eigenvalue="<<lowestEV<<" LeftT(0,0)="<<safe<<endl;
//                             }
//                             else
//                             {
//                                 linearStabilityBase-=stabilizerScale;
//                                 stabilityBase+=stabilizerScale;
//                                 app_log()<<"Probably will not converge: Jas Eigenvalue="<<lowestEV<<" LeftT(0,0)="<<safe<<endl;
//                             }
//                             maintain same number of "good" stability tries
//                             stability-=1;
//                             continue;
//                         }
//                     }
//                     
//                     myTimers[2]->start();
//                     lowestEV =getLowestEigenvector(LeftT,RightT,currentParameterDirections);
//                     myTimers[2]->stop();
//                 }
//                 else
//                 {
                    if (XS==0)
                    {
                      XS     = std::exp(stabilityBase +  stability*od_largest/nstabilizers);
                      for (int i=1; i<N; i++) LeftT(i,i) += XS;
                    }
                    else
                    {
                      //else XS is from the quartic fit
                      for (int i=1; i<N; i++) LeftT(i,i) += std::exp(XS);
                    }
                    
                    myTimers[2]->start();
//                     lowestEV =getLowestEigenvector(LeftT,RightT,currentParameterDirections);
                    lowestEV =getLowestEigenvector(LeftT,currentParameterDirections);
                    myTimers[2]->stop();
//                 }

                if (tooLow(safe,lowestEV))
                {
                    tooManyTries--;
                    if (tooManyTries>0)
                    {
                      if (stability==0)
                      {
                        app_log()<<"Probably will not converge: Eigenvalue="<<lowestEV<<" LeftT(0,0)="<<safe<<endl;
                        //try a larger stability base and repeat
                        stabilityBase+=stabilizerScale;
                        //maintain same number of "good" stability tries
                        stability-=1;
                      }
                      else
                      {
                        app_log()<<"Probably will not converge: Eigenvalue="<<lowestEV<<" LeftT(0,0)="<<safe<<endl;
                        //try a larger stability base and repeat
                        stabilityBase-=0.66*stabilizerScale;
                        //maintain same number of "good" stability tries
                        stability-=1;               
                      }
                    }
                    else
                    {
                      app_log()<<"Too many tries: Moving on to next step"<<endl;
                      stability=nstabilizers;
                    }
                    
                    continue;
                }

                if (MinMethod=="rescale")
                {
//                   method from umrigar
                    myTimers[3]->start();
                    Lambda = H2rescale*getNonLinearRescale(currentParameterDirections,S);
                    myTimers[3]->stop();

                    RealType bigVec(0);
                    for (int i=0; i<numParams; i++) bigVec = std::max(bigVec,std::abs(currentParameterDirections[i+1]));
                    if (Lambda*bigVec>bigChange)
                    {
                        app_log()<<"  Failed Step. Largest parameter change: "<<Lambda*bigVec<<endl;
                        tooManyTries--;
                        if (tooManyTries>0)
                        {
                            stabilityBase+=stabilizerScale;
                            stability-=1;
                            app_log()<<" Re-run with larger stabilityBase"<<endl;
                            continue;
                        }
                    }
                    else
                        for (int i=0; i<numParams; i++) optTarget->Params(i) = currentParameters[i] + Lambda*currentParameterDirections[i+1];
                }
                else
                {
                    //eigenCG part
                    for (int ldi=0; ldi < std::min(eigCG,int(LastDirections.size())); ldi++)
                    {
                        RealType nrmold(0), ovlpold(0);
                        for (int i=1; i<N; i++) nrmold += LastDirections[ldi][i]*LastDirections[ldi][i];
                        for (int i=1; i<N; i++) ovlpold += LastDirections[ldi][i]*currentParameterDirections[i];
                        ovlpold*=1.0/nrmold;
                        for (int i=1; i<N; i++) currentParameterDirections[i] -= ovlpold * LastDirections[ldi][i];
                    }

                    
                    //if we chose to "freeze" the CSF solutions at their minimum 
                    //  then we must add them in to the fixed part of the parameter changes
                    for (int i=0; i<numParams; i++) optparm[i] = currentParameters[i] + GEVSplitParameters[i];
                    for (int i=0; i<numParams; i++) optdir[i] = currentParameterDirections[i+1];
                    RealType bigVec(0);
                    for (int i=0; i<numParams; i++) bigVec = std::max(bigVec,std::abs(optdir[i]));

                    TOL = 1e-8/bigVec;

                    largeQuarticStep=bigChange/bigVec;
                    if (savedQuadstep>0)
                      quadstep=savedQuadstep/bigVec;
                    else if (deltaPrms>0)
                      quadstep=deltaPrms/bigVec;
                    else 
                      quadstep = 0.5*H2rescale*getNonLinearRescale(currentParameterDirections,S);
                    
//                  initial guess for line min bracketing
                    LambdaMax = 0.1/bigVec;
                    
                    myTimers[3]->start();
                    if (MinMethod=="quartic")
                      Valid=lineoptimization();
                    else if (MinMethod=="quartic_u")
                    {
                      int npts(9); int offset(2);
                      quadstep *= 2.0/npts;
                      Valid=lineoptimization3(npts,offset);
                    }
                    else Valid=lineoptimization2();
                    myTimers[3]->stop();

                    RealType biggestParameterChange = bigVec*std::abs(Lambda);
                    if ( (!Valid) || (biggestParameterChange>bigChange ))
                    {
                        app_log()<<"  Failed Step. Largest parameter change:"<<biggestParameterChange<<endl;
//                     optTarget->printEstimates();
                        tooManyTries--;
                        if (tooManyTries>0)
                        {
                            stabilityBase+=stabilizerScale;
                            stability-=1;
                            app_log()<<" Re-run with larger stabilityBase"<<endl;
                            continue;
                        }
                    }
                    else for (int i=0; i<numParams; i++) optTarget->Params(i) = optparm[i] + Lambda * optdir[i];
                    //Save this value in here for later
                    Lambda = biggestParameterChange;
                }
                //get cost at new minimum
                newCost = optTarget->Cost(false);

                // mmorales
                Valid=optTarget->IsValid;
                if (!ValidCostFunction(Valid)) continue;

                if (StabilizerMethod=="fit")
                {
                    std::pair<RealType,RealType> ms;
                    ms.first=newCost;
//                     the log fit seems to work best
                    ms.second=std::log(XS);
                    mappedStabilizers.push_back(ms);
                }


                app_log()<<" OldCost: "<<lastCost<<" NewCost: "<<newCost<<endl;
                optTarget->printEstimates();
//                 quit if newcost is greater than lastcost. E(Xs) looks quadratic (between steepest descent and parabolic)

                if ((newCost < lastCost)&&(newCost==newCost))
                {
                    //Move was acceptable
                    for (int i=0; i<numParams; i++) bestParameters[i] = optTarget->Params(i);
                    lastCost=newCost;
                    BestDirection=currentParameterDirections;
                    acceptedOneMove=true;

                    deltaPrms= Lambda;
                }
                else if (newCost>lastCost+1.0e-4)
                {
                    int neededForGoodQuarticFit=3;
                    if ((StabilizerMethod=="fit")&&(stability+1 < neededForGoodQuarticFit))
                    {
                        app_log()<<"Small change, but need "<< neededForGoodQuarticFit+1 <<" values for a good quartic stability fit."<<endl;
                    }
                    else if ((StabilizerMethod=="fit")&&(stability+1 >= neededForGoodQuarticFit))
                    {
                        stability = max(nstabilizers-2,stability);
                        if (stability==nstabilizers-2) app_log()<<"Small change, moving on to quartic fit."<<endl;
                        else app_log()<<"Moving on to next eigCG or iteration."<<endl;
                    }
                    else
                    {
                        stability = nstabilizers;
                        app_log()<<"Small change, moving on to next eigCG or iteration."<<endl;
                    }
                }
            }

            if (acceptedOneMove)
            {
                for (int i=0; i<numParams; i++) optTarget->Params(i) = bestParameters[i];
                currentParameters=bestParameters;
                LastDirections.push_back(BestDirection);
//             app_log()<< " Wave Function Parameters updated."<<endl;
//             optTarget->reportParameters();
            }
            else
            {
                for (int i=0; i<numParams; i++) optTarget->Params(i) = currentParameters[i];
                tries=TotalCGSteps;
            }

        }
    }
    finish();
    return (optTarget->getReportCounter() > 0);
}

/** Parses the xml input file for parameter definitions for the wavefunction optimization.
* @param q current xmlNode
* @return true if successful
*/
bool
QMCFixedSampleLinearOptimize::put(xmlNodePtr q)
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
//         else if (cname == "optimizer")
//         {
//             xmlChar* att= xmlGetProp(cur,(const xmlChar*)"method");
//             if (att)
//             {
//                 optmethod = (const char*)att;
//             }
//             optNode=cur;
//         }
//         else if (cname == "optimize")
//         {
//             xmlChar* att= xmlGetProp(cur,(const xmlChar*)"method");
//             if (att)
//             {
//                 optmethod = (const char*)att;
//             }
//         }
        cur=cur->next;
    }
    //no walkers exist, add 10
    if (W.getActiveWalkers() == 0) addWalkers(omp_get_max_threads());

    NumOfVMCWalkers=W.getActiveWalkers();

    //create VMC engine
    if (vmcEngine ==0)
    {
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
        vmcEngine = new VMCSingleOMP(W,Psi,H,hamPool,psiPool);
        vmcEngine->setUpdateMode(vmcMove[0] == 'p');
        vmcEngine->initCommunicator(myComm);
    }
    vmcEngine->setStatus(RootName,h5FileRoot,AppendRun);
    vmcEngine->process(qsave);

    bool success=true;

    if (optTarget == 0)
    {
#if defined (QMC_CUDA)
        if (useGPU == "yes")
            optTarget = new QMCCostFunctionCUDA(W,Psi,H,hamPool);
        else
#endif
// #if defined(ENABLE_OPENMP)
//             if (omp_get_max_threads()>1)
//             {
//                 optTarget = new QMCCostFunctionOMP(W,Psi,H,hamPool);
//             }
//             else
// #endif
//                 optTarget = new QMCCostFunctionSingle(W,Psi,H);
        optTarget = new QMCCostFunctionOMP(W,Psi,H,hamPool);
        optTarget->setStream(&app_log());
        success=optTarget->put(q);
    }
    return success;
}

}
/***************************************************************************
* $RCSfile$   $Author: jnkim $
* $Revision: 1286 $   $Date: 2006-08-17 12:33:18 -0500 (Thu, 17 Aug 2006) $
* $Id: QMCFixedSampleLinearOptimize.cpp 1286 2006-08-17 17:33:18Z jnkim $
***************************************************************************/

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
#include "QMCDrivers/QMCCSLinearOptimize.h"
#include "Particle/HDFWalkerIO.h"
#include "Particle/DistanceTable.h"
#include "OhmmsData/AttributeSet.h"
#include "Message/CommOperators.h"
#include "QMCApp/HamiltonianPool.h"
#include "Numerics/Blasf.h"
#include <cassert>
#include "Numerics/LinearFit.h"
#include <iostream>
#include <fstream>

/*#include "Message/Communicate.h"*/

namespace qmcplusplus
{

QMCCSLinearOptimize::QMCCSLinearOptimize(MCWalkerConfiguration& w,
        TrialWaveFunction& psi, QMCHamiltonian& h, HamiltonianPool& hpool): QMCDriver(w,psi,h),
        PartID(0), NumParts(1), WarmupBlocks(10), SkipSampleGeneration("no"), hamPool(hpool),
        optTarget(0), vmcEngine(0), Max_iterations(1), wfNode(NULL), optNode(NULL), allowedCostDifference(1.0e-6),
        exp0(-16), exp1(0),  nstabilizers(10), stabilizerScale(0.5), bigChange(1), eigCG(1), TotalCGSteps(2), w_beta(0.0),
        MinMethod("quartic"), GEVtype("mixed"), StabilizerMethod("best"), GEVSplit("no")
{
    //set the optimization flag
    QMCDriverMode.set(QMC_OPTIMIZE,1);
    //read to use vmc output (just in case)
    RootName = "pot";
    QMCType ="QMCCSLinearOptimize";
    optmethod = "Linear";
    m_param.add(WarmupBlocks,"warmupBlocks","int");
    m_param.add(SkipSampleGeneration,"skipVMC","string");
    m_param.add(Max_iterations,"max_its","int");
    m_param.add(nstabilizers,"nstabilizers","int");
    m_param.add(stabilizerScale,"stabilizerscale","double");
    m_param.add(allowedCostDifference,"alloweddifference","double");
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
QMCCSLinearOptimize::~QMCCSLinearOptimize()
{
    delete vmcEngine;
    delete optTarget;
}

void QMCCSLinearOptimize::add_timers(vector<NewTimer*>& timers)
{
    timers.push_back(new NewTimer("QMCCSLinearOptimize::GenerateSamples"));
    timers.push_back(new NewTimer("QMCCSLinearOptimize::Initialize"));
    timers.push_back(new NewTimer("QMCCSLinearOptimize::Eigenvalue"));
    timers.push_back(new NewTimer("QMCCSLinearOptimize::Line_Minimization"));
    timers.push_back(new NewTimer("QMCCSLinearOptimize::GradCost"));
    for (int i=0; i<timers.size(); ++i) TimerManager.addTimer(timers[i]);
}

QMCCSLinearOptimize::RealType QMCCSLinearOptimize::Func(RealType dl)
{
    for (int i=0; i<optparm.size(); i++) optTarget->Params(i) = optparm[i] + dl*optdir[i];
    QMCCSLinearOptimize::RealType c = optTarget->Cost(false);
    //only allow this to go false if it was true. If false, stay false
    if (validFuncVal) validFuncVal = optTarget->IsValid;
    return c;
}

/** Add configuration files for the optimization
* @param a root of a hdf5 configuration file
*/
void QMCCSLinearOptimize::addConfiguration(const string& a)
{
    if (a.size()) ConfigFile.push_back(a);
}

void QMCCSLinearOptimize::start()
{
    optTarget->initCommunicator(myComm);

    //generate Matrix
    myTimers[0]->start();
    generateSamples();
    myTimers[0]->stop();

    app_log() << "<opt stage=\"setup\">" << endl;
    app_log() << "  <log>"<<endl;

    //reset the rootname
    optTarget->setRootName(RootName);
    optTarget->setWaveFunctionNode(wfNode);

    app_log() << "  </log>"<<endl;
    app_log() << "</opt>" << endl;

    app_log() << "  <log>" << endl;
    app_log()<<"  GEV method "<<GEVtype<<endl;
    app_log()<<"  Split EV   "<<GEVSplit<<endl;
    app_log()<<"  Line Minimization method "<<MinMethod<<endl;

}

bool QMCCSLinearOptimize::ValidCostFunction(bool valid)
{
    if (!valid) app_log()<<" Cost Function is Invalid. If this frequently, try reducing the step size of the line minimization or reduce the number of cycles. " <<endl;
    return valid;
}

bool QMCCSLinearOptimize::run()
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
    vector<vector<RealType> > savedCSparameters;
    optdir.resize(numParams,0);
    optparm.resize(numParams,0);
    for (int i=0; i<numParams; i++) currentParameters[i] = optTarget->Params(i);
    savedCSparameters.push_back(currentParameters);
    
    Matrix<RealType> Ham(N,N);
    Matrix<RealType> Ham2(N,N);
    Matrix<RealType> Var(N,N);
    Matrix<RealType> S(N,N);
    vmcEngine->fillMatrices(Ham2,Ham,Var,S);
    RealType lastCost(vmcEngine->Cost());
    RealType  newCost(lastCost);

    vector<RealType> bestParameters(currentParameters);
    vector<RealType> GEVSplitParameters(numParams,0);


//this is the small amount added to the diagonal to stabilize the eigenvalue equation. 10^stabilityBase

    RealType stabilityBase(exp0);
    //This is the amount we add to the linear parameters
    RealType linearStabilityBase(exp1);
    int tooManyTries(200);

    Matrix<RealType> Left(N,N);
    Matrix<RealType> Right(N,N);

    vector<std::pair<RealType,RealType> > mappedStabilizers;
    if (nstabilizers<3)
    {
        if (StabilizerMethod=="fit") app_log()<<" Need 3 stabilizers minimum for the fit"<<endl;
        StabilizerMethod="best";
    }

    if (GEVtype=="H2")
    {
        Left=Ham;
        RealType H2rescale=1.0/Ham2(0,0);
        Right=(1-w_beta)*S + w_beta*H2rescale*Ham2;
    }
    else
    {
        Right=S;
        Left=(1.0-w_beta)*Ham + w_beta*Var;
    }
    
//     Invert(Right.data(), N, N);

    //Find largest off-diagonal element compared to diagonal element.
    //This gives us an idea how well conditioned it is and can be used to stabilize.
    RealType od_largest(0);
    for (int i=0; i<N; i++) for (int j=0; j<N; j++)
            od_largest=std::max( std::max(od_largest,std::abs(Left(i,j))-std::abs(Left(i,i))), std::abs(Left(i,j))-std::abs(Left(j,j)));
    if (od_largest>0) od_largest = std::log(od_largest);
    else od_largest=1.0;

    RealType safe = Left(0,0);
    if (StabilizerMethod=="cs") nstabilizers = (omp_get_max_threads()-1)*(nstabilizers/(omp_get_max_threads()-1));
    for (int stability=0; stability<nstabilizers; stability++)
    {
      app_log()<<"Iteration: "<<stability+1<<"/"<<nstabilizers<<endl;
        Matrix<RealType> LeftT(N,N), RightT(N,N);
        for (int i=0; i<N; i++)
            for (int j=0; j<N; j++)
            {
                LeftT(i,j)= Left(j,i);
                RightT(i,j)= Right(j,i);
            }


        RealType XS(0);
        if ((StabilizerMethod=="fit")&&(stability==nstabilizers-1))
        {
            int nms=mappedStabilizers.size();
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
              stability=nstabilizers;
            }
            else
            {//Quadratic fit the stabilizers we have tried and try to choose the best we can
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
              stability=nstabilizers;
            }
        }

        RealType lowestEV(0);
        if ((GEVSplit=="rescale")||(GEVSplit=="freeze"))
        {
//These are experimental and aren't very good.
//dummy bool
            bool CSF_lower(true);
            lowestEV = getSplitEigenvectors(first,last,LeftT,RightT,currentParameterDirections,GEVSplitParameters,GEVSplit,CSF_lower);
        }
        else if (GEVSplit=="stability") //This seems to work pretty well.
        {
            if (XS==0)
            {
                XS  = std::exp(stabilityBase + stability*od_largest/nstabilizers);
                for (int i=first; i<last; i++) LeftT(i+1,i+1) += XS;

                RealType XS_lin = std::exp(linearStabilityBase + stabilityBase + stability*od_largest/nstabilizers);
                if (first==0) for (int i=last; i<N; i++) LeftT(i+1,i+1) += XS_lin;
                else for (int i=0; i<first; i++) LeftT(i+1,i+1) += XS_lin;
            }
            else //else XS is from the quartic fit
            {
//Not sure how to control for the quartic fit and the two different stabilizers. This seems ok.
//Better algorithm exists?
                for (int i=first; i<last; i++) LeftT(i+1,i+1) += std::exp(XS);

                RealType XS_lin = std::exp(linearStabilityBase+XS);
                if (first==0) for (int i=last; i<N; i++) LeftT(i+1,i+1) += XS_lin;
                for (int i=0; i<first; i++) LeftT(i+1,i+1) += XS_lin;
            }

            if (stability==0)
            {
//  Only need to do this the first time we step into the routine
                bool CSF_lower(true);
                lowestEV=getSplitEigenvectors(first,last,LeftT,RightT,currentParameterDirections,GEVSplitParameters,GEVSplit,CSF_lower);
                if (tooLow(safe,lowestEV))
                {
                    if (CSF_lower)
                    {
                        linearStabilityBase+=stabilizerScale;
                        app_log()<<"Probably will not converge:\n  CSF Ev="<<lowestEV<<" LeftT(0,0)="<<safe<<" exp0: "<<stabilityBase<<" exp1: "<<linearStabilityBase<<endl;
                    }
                    else
                    {
                        linearStabilityBase-=stabilizerScale;
                        stabilityBase+=stabilizerScale;
                        app_log()<<"Probably will not converge:\n  Jas Ev="<<lowestEV<<" LeftT(0,0)="<<safe<<" exp0: "<<stabilityBase<<" exp1: "<<linearStabilityBase<<endl;
                    }
//maintain same number of "good" stability tries
                    stability-=1;
                    continue;
                }
            }

            myTimers[2]->start();
            lowestEV =getLowestEigenvector(LeftT,RightT,currentParameterDirections);
            myTimers[2]->stop();
        }
        else
        {
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
            lowestEV =getLowestEigenvector(LeftT,RightT,currentParameterDirections);
            myTimers[2]->stop();
        }

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
        
        RealType Lambda_Last(Lambda);
        myTimers[3]->start();
        if (GEVtype=="H2")
        {
//           the rescaling isn't right for this
          RealType bigVec(0);
          for (int i=0; i<numParams; i++) bigVec = std::max(bigVec,std::abs(currentParameterDirections[i+1]));
          Lambda=0.5*bigChange/bigVec;
        }
          else
        Lambda = getNonLinearRescale(currentParameterDirections,S);
        myTimers[3]->stop();
//         app_log()<<"Computed Lambda is: "<<Lambda<<endl;

        if (MinMethod=="rescale")
        {
//                   method from umrigar
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
//             app_log()<<" Umrigar Lambda: "<<Lambda<<endl;
            optTarget->resetPsi(false);
            vmcEngine->clearComponentMatrices();
            vmcEngine->run(false);
            newCost = vmcEngine->Cost();
//             app_log()<<lastCost<<" "<<newCost<<endl;
        }
        else
        {
//here is correlated sampling routine
             int nthreads = omp_get_max_threads();
             std::vector<RealType> lambdas(nthreads);
             if (savedCSparameters.size()>1)
             {
               if (Lambda_Last>0)
                 for(int i=0;i<nthreads;i++) lambdas[i] = i/(nthreads-1.0)*Lambda_Last;
               else if (Lambda_Last<0)
                 for(int i=0;i<nthreads;i++) lambdas[i] = Lambda_Last - i/(nthreads-1.0)*Lambda_Last;
               else
                 for(int i=0;i<nthreads;i++) lambdas[i] = i/(nthreads-1.0)*Lambda;
             }
             else for(int i=0;i<nthreads;i++) lambdas[i] = i/(nthreads-1.0)*Lambda;
             
             newCost = vmcEngine->runCS(currentParameters,currentParameterDirections,lambdas);
             Lambda=lambdas[0];
             if (Lambda==0)
             {
               app_log()<<"  Failed Step. Lambda=0: "<<endl;
               tooManyTries--;
               if (tooManyTries>0)
               {
                   stabilityBase+=stabilizerScale;
                   stability-=1;
                   app_log()<<" Re-run with larger stabilityBase"<<endl;
                   continue;
               }
             }
             vector<RealType> cs(numParams);
             for (int i=0; i<numParams; i++) cs[i] = optTarget->Params(i)=currentParameters[i] + Lambda*currentParameterDirections[i+1];
             optTarget->resetPsi(false);
             savedCSparameters.push_back(cs);
        }


        if (StabilizerMethod=="fit")
        {
            std::pair<RealType,RealType> ms;
            ms.first=newCost;
//                     the log fit seems to work best
            ms.second=std::log(XS);
            mappedStabilizers.push_back(ms);
        }

        if (StabilizerMethod=="cs")
        {
          if(savedCSparameters.size()==omp_get_max_threads())
          {
            app_log()<<"   Choosing best"<<endl;
            RealType error;
            int bestP = vmcEngine->runCS(savedCSparameters,error);
            if (bestP<0)
            {
              app_log()<<"   Error in CS cost function. Unchanged parameters."<<endl;
              bestP=0;
            }
            vector<RealType> cs(numParams);
            for (int i=0; i<numParams; i++) cs[i] = bestParameters[i] = savedCSparameters[bestP][i];
            savedCSparameters.clear();
            savedCSparameters.push_back(cs);
            Lambda=Lambda_Last;
          }
        }
        else if (newCost < lastCost)
        {
//Move was acceptable
            for (int i=0; i<numParams; i++) bestParameters[i] = optTarget->Params(i);
            lastCost=newCost;
        }
        else if (newCost>lastCost)
        {
          int neededForGoodFit=2;//really one more so if 5, then 6 values kept. 3 is minimum.
          if (stability ==0 )
          {
              stability-=1;
              stabilityBase+=stabilizerScale;
          }
          else if ((StabilizerMethod=="fit")&&(stability < neededForGoodFit))
          {
              app_log()<<"Small change, but need "<< neededForGoodFit+1 <<" values for a good quartic stability fit."<<endl;
          }
          else if ((StabilizerMethod=="fit")&&(stability >= neededForGoodFit))
          {
              stability = max(nstabilizers-2,stability);
              if (stability==nstabilizers-2) app_log()<<"Small change, moving on to fit."<<endl;
              else app_log()<<"Small change, moving on to next iteration."<<endl;
          }
          else
          {
              stability = nstabilizers;
              app_log()<<"Small change, moving on to next iteration."<<endl;
          }
        }
    }
    
    for (int i=0; i<numParams; i++) optTarget->Params(i) = bestParameters[i];
    vmcEngine->clearComponentMatrices();


    
    finish();
    if (tooManyTries==0) return false;
    return true;
}

void QMCCSLinearOptimize::finish()
{
    MyCounter++;
    TimerManager.print(myComm);
    TimerManager.reset();
    Lambda=0;
    app_log() << "  </log>" << endl;
    optTarget->reportParameters();
    app_log() << "</opt>" << endl;
    app_log() << "</optimization-report>" << endl;
}

void QMCCSLinearOptimize::generateSamples()
{

    app_log() << "<optimization-report>" << endl;

    if (W.getActiveWalkers()>NumOfVMCWalkers)
    {
        W.destroyWalkers(W.getActiveWalkers()-NumOfVMCWalkers);
        app_log() << "  QMCCSLinearOptimize::generateSamples removed walkers." << endl;
        app_log() << "  Number of Walkers per node " << W.getActiveWalkers() << endl;
    }

    vmcEngine->QMCDriverMode.set(QMC_OPTIMIZE,1);
    vmcEngine->QMCDriverMode.set(QMC_WARMUP,0);

    //vmcEngine->setValue("recordWalkers",1);//set record
    vmcEngine->setValue("current",0);//reset CurrentStep
    t1.restart();
    //     W.reset();
    branchEngine->flush(0);
    branchEngine->reset();
    vmcEngine->run();
    app_log() << "  Execution time = " << t1.elapsed() << endl;
    app_log() << "</vmc>" << endl;

    //write parameter history and energies to the parameter file in the trial wave function through opttarget
    RealType e,w,var;
    vmcEngine->Estimators->getEnergyAndWeight(e,w,var);
    optTarget->recordParametersToPsi(e,var);

    //branchEngine->Eref=vmcEngine->getBranchEngine()->Eref;
//         branchEngine->setTrialEnergy(vmcEngine->getBranchEngine()->getEref());
    //set the h5File to the current RootName
    h5FileRoot=RootName;
}

QMCCSLinearOptimize::RealType QMCCSLinearOptimize::getLowestEigenvector(Matrix<RealType>& A, Matrix<RealType>& B, vector<RealType>& ev)
{
    int Nl(ev.size());
    //Tested the single eigenvalue speed and It was no faster.
    //segfault issues with single eigenvalue problem for some machines
    bool singleEV(false);
    if (singleEV)
    {
        Matrix<double> TAU(Nl,Nl);
        int INFO;
        int LWORK(-1);
        vector<RealType> WORK(1);
        //optimal work size
        dgeqrf( &Nl, &Nl, B.data(), &Nl, TAU.data(), &WORK[0], &LWORK, &INFO);
        LWORK=int(WORK[0]);
        WORK.resize(LWORK);
        //QR factorization of S, or H2 matrix. to be applied to H before solve.
        dgeqrf( &Nl, &Nl, B.data(), &Nl, TAU.data(), &WORK[0], &LWORK, &INFO);

        char SIDE('L');
        char TRANS('T');
        LWORK=-1;
        //optimal work size
        dormqr(&SIDE, &TRANS, &Nl, &Nl, &Nl, B.data(), &Nl, TAU.data(), A.data(), &Nl, &WORK[0], &LWORK, &INFO);
        LWORK=int(WORK[0]);
        WORK.resize(LWORK);
        //Apply Q^T to H
        dormqr(&SIDE, &TRANS, &Nl, &Nl, &Nl, B.data(), &Nl, TAU.data(), A.data(), &Nl, &WORK[0], &LWORK, &INFO);

        //now we have a pair (A,B)=(Q^T*H,Q^T*S) where B is upper triangular and A is general matrix.
        //reduce the matrix pair to generalized upper Hesenberg form
        char COMPQ('N'), COMPZ('I');
        int ILO(1);
        int LDQ(Nl);
        Matrix<double> Z(Nl,Nl), Q(Nl,LDQ); //starts as unit matrix
        for (int zi=0; zi<Nl; zi++) Z(zi,zi)=1;
        dgghrd(&COMPQ, &COMPZ, &Nl, &ILO, &Nl, A.data(), &Nl, B.data(), &Nl, Q.data(), &LDQ, Z.data(), &Nl, &INFO);

        //Take the pair and reduce to shur form and get eigenvalues
        vector<RealType> alphar(Nl),alphai(Nl),beta(Nl);
        char JOB('S');
        COMPQ='N';
        COMPZ='V';
        LWORK=-1;
        //get optimal work size
        dhgeqz(&JOB, &COMPQ, &COMPZ, &Nl, &ILO, &Nl, A.data(), &Nl, B.data(), &Nl, &alphar[0], &alphai[0], &beta[0], Q.data(), &LDQ, Z.data(), &Nl, &WORK[0], &LWORK, &INFO);
        LWORK=int(WORK[0]);
        WORK.resize(LWORK);
        dhgeqz(&JOB, &COMPQ, &COMPZ, &Nl, &ILO, &Nl, A.data(), &Nl, B.data(), &Nl, &alphar[0], &alphai[0], &beta[0], Q.data(), &LDQ, Z.data(), &Nl, &WORK[0], &LWORK, &INFO);
        //find the best eigenvalue
        vector<std::pair<RealType,int> > mappedEigenvalues(Nl);
        for (int i=0; i<Nl; i++)
        {
            RealType evi(alphar[i]/beta[i]);
            if (abs(evi)<1e10)
            {
                mappedEigenvalues[i].first=evi;
                mappedEigenvalues[i].second=i;
            }
            else
            {
                mappedEigenvalues[i].first=1e100;
                mappedEigenvalues[i].second=i;
            }
        }
        std::sort(mappedEigenvalues.begin(),mappedEigenvalues.end());
        int BestEV(mappedEigenvalues[0].second);

//                   now we rearrange the  the matrices
        if (BestEV!=0)
        {
            bool WANTQ(false);
            bool WANTZ(true);
            int ILST(1);
            int IFST(BestEV+1);
            LWORK=-1;

            dtgexc(&WANTQ, &WANTZ, &Nl, A.data(), &Nl, B.data(), &Nl, Q.data(), &Nl, Z.data(), &Nl, &IFST, &ILST, &WORK[0], &LWORK, &INFO);
            LWORK=int(WORK[0]);
            WORK.resize(LWORK);
            dtgexc(&WANTQ, &WANTZ, &Nl, A.data(), &Nl, B.data(), &Nl, Q.data(), &Nl, Z.data(), &Nl, &IFST, &ILST, &WORK[0], &LWORK, &INFO);
        }
        //now we compute the eigenvector
        SIDE='R';
        char HOWMNY('S');
        int M(0);
        Matrix<double> Z_I(Nl,Nl);
        bool SELECT[Nl];
        for (int zi=0; zi<Nl; zi++) SELECT[zi]=false;
        SELECT[0]=true;

        WORK.resize(6*Nl);
        dtgevc(&SIDE, &HOWMNY, &SELECT[0], &Nl, A.data(), &Nl, B.data(), &Nl, Q.data(), &LDQ, Z_I.data(), &Nl, &Nl, &M, &WORK[0], &INFO);

        std::vector<RealType> evec(Nl,0);
        for (int i=0; i<Nl; i++) for (int j=0; j<Nl; j++) evec[i] += Z(j,i)*Z_I(0,j);
        for (int i=0; i<Nl; i++) ev[i] = evec[i]/evec[0];
//     for (int i=0; i<Nl; i++) app_log()<<ev[i]<<" ";
//     app_log()<<endl;
        return mappedEigenvalues[0].first;
    }
    else
    {
// OLD ROUTINE. CALCULATES ALL EIGENVECTORS
//   Getting the optimal worksize
        char jl('N');
        char jr('V');
        vector<RealType> alphar(Nl),alphai(Nl),beta(Nl);
        Matrix<RealType> eigenT(Nl,Nl);
        int info;
        int lwork(-1);
        vector<RealType> work(1);

        RealType tt(0);
        int t(1);
        dggev(&jl, &jr, &Nl, A.data(), &Nl, B.data(), &Nl, &alphar[0], &alphai[0], &beta[0],&tt,&t, eigenT.data(), &Nl, &work[0], &lwork, &info);
        lwork=int(work[0]);
        work.resize(lwork);

        //~ //Get an estimate of E_lin
        //~ Matrix<RealType> H_tmp(HamT);
        //~ Matrix<RealType> S_tmp(ST);
        //~ dggev(&jl, &jr, &Nl, H_tmp.data(), &Nl, S_tmp.data(), &Nl, &alphar[0], &alphai[0], &beta[0],&tt,&t, eigenT.data(), &Nl, &work[0], &lwork, &info);
        //~ RealType E_lin(alphar[0]/beta[0]);
        //~ int e_min_indx(0);
        //~ for (int i=1; i<Nl; i++)
        //~ if (E_lin>(alphar[i]/beta[i]))
        //~ {
        //~ E_lin=alphar[i]/beta[i];
        //~ e_min_indx=i;
        //~ }
        dggev(&jl, &jr, &Nl, A.data(), &Nl, B.data(), &Nl, &alphar[0], &alphai[0], &beta[0],&tt,&t, eigenT.data(), &Nl, &work[0], &lwork, &info);
        if (info!=0)
        {
            APP_ABORT("Invalid Matrix Diagonalization Function!");
        }

        vector<std::pair<RealType,int> > mappedEigenvalues(Nl);
        for (int i=0; i<Nl; i++)
        {
            RealType evi(alphar[i]/beta[i]);
            if (abs(evi)<1e10)
            {
                mappedEigenvalues[i].first=evi;
                mappedEigenvalues[i].second=i;
            }
            else
            {
                mappedEigenvalues[i].first=1e100;
                mappedEigenvalues[i].second=i;
            }
        }
        std::sort(mappedEigenvalues.begin(),mappedEigenvalues.end());

        for (int i=0; i<Nl; i++) ev[i] = eigenT(mappedEigenvalues[0].second,i)/eigenT(mappedEigenvalues[0].second,0);
        return mappedEigenvalues[0].first;
    }
}


bool QMCCSLinearOptimize::nonLinearRescale(std::vector<RealType>& dP, Matrix<RealType> S)
{
    RealType rescale = getNonLinearRescale(dP,S);
    for (int i=1; i<dP.size(); i++) dP[i]*=rescale;
    return true;
}


void QMCCSLinearOptimize::getNonLinearRange(int& first, int& last)
{
    std::vector<int> types;
    optTarget->getParameterTypes(types);
    first=0;
    last=types.size();

    //assume all non-linear coeffs are together.
    if (types[0]==optimize::LINEAR_P)
    {
        int i(0);
        while ((types[i]==optimize::LINEAR_P)&&(i<types.size())) i++;
        first=i;
    }
    else
    {
        int i(1);
        while ((types[i]!=optimize::LINEAR_P)&&(i<types.size())) i++;
        last=i;
    }
//     returns the number of non-linear parameters.
//     app_log()<<first<<" "<<last<<endl;
}

QMCCSLinearOptimize::RealType QMCCSLinearOptimize::getNonLinearRescale(std::vector<RealType>& dP, Matrix<RealType> S)
{
    int first(0),last(0);
    getNonLinearRange(first,last);
    if (first==last) return 1.0;

    RealType D(1.0);
    for (int i=first; i<last; i++) for (int j=first; j<last; j++) D += S(j+1,i+1)*dP[i+1]*dP[j+1];

    D = std::sqrt(std::abs(D));
    RealType xi(0.5);
//for the jchem126 paper
//     vector<RealType> N_i(last-first,0);
//     vector<RealType> M_i(last-first,0);
//     for (int i=0; i<last-first; i++)
//     {
//         M_i[i] = xi*D +(1-xi);
//         RealType tsumN(0);
//         for (int j=first; j<last; j++)
//         {
//             tsumN += S(i+first+1,j+1)*dP[j+1];
//         }
//         N_i[i] += (1-xi)*tsumN;
//         N_i[i] *= -1.0/M_i[i];
//     }

//from the prl98
    D = 1.0/((1.0-xi) + xi*D);
    vector<RealType> N_i(dP.size()-1,0);
    for (int j=0; j<dP.size()-1; j++) for (int i=first; i<last; i++) N_i[j]+=S(j+1,i+1)*dP[i+1];
    for (int j=0; j<dP.size()-1; j++) N_i[j] *= D;
    
    RealType rescale(1);
    for (int j=0; j<dP.size()-1; j++) rescale += N_i[j]*dP[j+1];
    rescale = 1.0/rescale;
    return rescale;
}

QMCCSLinearOptimize::RealType QMCCSLinearOptimize::getSplitEigenvectors(int first, int last, Matrix<RealType>& FullLeft, Matrix<RealType>& FullRight, vector<RealType>& FullEV, vector<RealType>& LocalEV, string CSF_Option, bool& CSF_scaled)
{
    vector<RealType> GEVSplitDirection(N,0);
    RealType returnValue;
    int N_nonlin=last-first;
    int N_lin   =N-N_nonlin-1;

//  matrices are one larger than parameter sets
    int M_nonlin=N_nonlin+1;
    int M_lin   =N_lin+1;
//  index mapping for the matrices
    int J_begin(first+1), J_end(last+1);
    int CSF_begin(1), CSF_end(first+1);
    if (first==0)
    {
        CSF_begin=last+1;
        CSF_end=N;
    }
//the Mini matrix composed of just the Nonlinear terms
    Matrix<RealType> LeftTJ(M_nonlin,M_nonlin), RightTJ(M_nonlin,M_nonlin);

//                     assume all jastrow parameters are together either first or last
    LeftTJ(0,0) =  FullLeft(0,0);
    RightTJ(0,0)= FullRight(0,0);
    for (int i=J_begin; i<J_end; i++)
    {
        LeftTJ(i-J_begin+1,0) =  FullLeft(i,0);
        RightTJ(i-J_begin+1,0)= FullRight(i,0);
        LeftTJ(0,i-J_begin+1) =  FullLeft(0,i);
        RightTJ(0,i-J_begin+1)= FullRight(0,i);
        for (int j=J_begin; j<J_end; j++)
        {
            LeftTJ(i-J_begin+1,j-J_begin+1) =  FullLeft(i,j);
            RightTJ(i-J_begin+1,j-J_begin+1)= FullRight(i,j);
        }
    }

    vector<RealType> J_parms(M_nonlin);
    myTimers[2]->start();
    RealType lowest_J_EV =getLowestEigenvector(LeftTJ,RightTJ,J_parms);
    myTimers[2]->stop();

//the Mini matrix composed of just the Linear terms
    Matrix<RealType> LeftTCSF(M_lin,M_lin), RightTCSF(M_lin,M_lin);

    LeftTCSF(0,0) =  FullLeft(0,0);
    RightTCSF(0,0)= FullRight(0,0);
    for (int i=CSF_begin; i<CSF_end; i++)
    {
        LeftTCSF(i-CSF_begin+1,0) =  FullLeft(i,0);
        RightTCSF(i-CSF_begin+1,0)= FullRight(i,0);
        LeftTCSF(0,i-CSF_begin+1) =  FullLeft(0,i);
        RightTCSF(0,i-CSF_begin+1)= FullRight(0,i);
        for (int j=CSF_begin; j<CSF_end; j++)
        {
            LeftTCSF(i-CSF_begin+1,j-CSF_begin+1) =  FullLeft(i,j);
            RightTCSF(i-CSF_begin+1,j-CSF_begin+1)= FullRight(i,j);
        }
    }
    vector<RealType> CSF_parms(M_lin);
    myTimers[2]->start();
    RealType lowest_CSF_EV =getLowestEigenvector(LeftTCSF,RightTCSF,CSF_parms);
    myTimers[2]->stop();

// //                   Now we have both eigenvalues and eigenvectors
//                   app_log()<<" Jastrow eigenvalue: "<<lowest_J_EV<<endl;
//                   app_log()<<"     CSF eigenvalue: "<<lowest_CSF_EV<<endl;

//                We can rescale the matrix and re-solve the whole thing or take the CSF parameters
//                  as solved in the matrix and opt the Jastrow instead

    if (CSF_Option=="freeze")
    {
        returnValue=min(lowest_J_EV,lowest_CSF_EV);
//                   Line minimize for the nonlinear components
        for (int i=J_begin; i<J_end; i++) GEVSplitDirection[i]=J_parms[i-J_begin+1];

//                   freeze the CSF components at this minimum
        for (int i=CSF_begin; i<CSF_end; i++) LocalEV[i-1]=CSF_parms[i-CSF_begin+1];

        FullEV[0]=1.0;
        for (int i=J_begin; i<J_end; i++) FullEV[i]=GEVSplitDirection[i];
    }
    else if (CSF_Option=="rescale")
    {
        RealType matrixRescaler=std::sqrt(std::abs(lowest_CSF_EV/lowest_J_EV));
        for (int i=0; i<N; i++)
            for (int j=0; j<N; j++)
            {
                if ((i>=J_begin)&&(i<J_end))
                {
                    FullLeft(i,j) *=matrixRescaler;
                    FullRight(i,j)*=matrixRescaler;
                }

                if ((j>=J_begin)&&(j<J_end))
                {
                    FullLeft(i,j) *=matrixRescaler;
                    FullRight(i,j)*=matrixRescaler;
                }
            }

        myTimers[2]->start();
        returnValue =getLowestEigenvector(FullLeft,FullRight,FullEV);
        myTimers[2]->stop();
    }
    else if (CSF_Option =="stability")
    {
//       just return the value of the CSF part
        if (lowest_J_EV>lowest_CSF_EV)
            CSF_scaled=true;
        else
            CSF_scaled=false;
        returnValue=min(lowest_J_EV,lowest_CSF_EV);
    }
    return returnValue;
}

/** Parses the xml input file for parameter definitions for the wavefunction optimization.
* @param q current xmlNode
* @return true if successful
*/
bool
QMCCSLinearOptimize::put(xmlNodePtr q)
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
    //no walkers exist, add 10
    if (W.getActiveWalkers() == 0) addWalkers(omp_get_max_threads());

    NumOfVMCWalkers=W.getActiveWalkers();

#if not defined(ENABLE_OPENMP)
    APP_ABORT(" CS LINEAR OPT: DESIGNED FOR OPENMP-MPI HYBRID. ");
#endif

    //create VMC engine
    if (vmcEngine ==0)
    {
        vmcEngine = new VMCLinearOptOMP(W,Psi,H,hamPool);
        vmcEngine->setUpdateMode(vmcMove[0] == 'p');
        vmcEngine->initCommunicator(myComm);
    }
    vmcEngine->setStatus(RootName,h5FileRoot,AppendRun);
    vmcEngine->process(qsave);

    bool success=true;
    if (optTarget == 0)
    {
        optTarget = new QMCCSLinearOptimizeWFmanagerOMP(W,Psi,H,hamPool);
        optTarget->setStream(&app_log());
        success=optTarget->put(q);
    }
    return success;
}

void QMCCSLinearOptimize::resetComponents(xmlNodePtr cur)
{
    m_param.put(cur);
    optTarget->put(cur);
    vmcEngine->resetComponents(cur);
}
}
/***************************************************************************
* $RCSfile$   $Author: jnkim $
* $Revision: 1286 $   $Date: 2006-08-17 12:33:18 -0500 (Thu, 17 Aug 2006) $
* $Id: QMCCSLinearOptimize.cpp 1286 2006-08-17 17:33:18Z jnkim $
***************************************************************************/

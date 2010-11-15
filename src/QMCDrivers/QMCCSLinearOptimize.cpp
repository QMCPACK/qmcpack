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
#include "Optimize/NRCOptimization.h"
#include <iostream>
#include <fstream>

/*#include "Message/Communicate.h"*/

namespace qmcplusplus
{

QMCCSLinearOptimize::QMCCSLinearOptimize(MCWalkerConfiguration& w,
        TrialWaveFunction& psi, QMCHamiltonian& h, HamiltonianPool& hpool, WaveFunctionPool& ppool): QMCDriver(w,psi,h),
        hamPool(hpool), psipool(ppool), optTarget(0), vmcEngine(0), Max_iterations(1), wfNode(NULL), optNode(NULL), 
        exp0(-16), nstabilizers(10), stabilizerScale(0.5), bigChange(1), w_beta(0.0),
        MinMethod("quartic"), GEVtype("mixed")
{
    //set the optimization flag
    QMCDriverMode.set(QMC_OPTIMIZE,1);
    //read to use vmc output (just in case)
    RootName = "pot";
    QMCType ="QMCCSLinearOptimize";
    optmethod = "Linear";
    m_param.add(Max_iterations,"max_its","int");
    m_param.add(nstabilizers,"nstabilizers","int");
    m_param.add(stabilizerScale,"stabilizerscale","double");
    m_param.add(bigChange,"bigchange","double");
    m_param.add(w_beta,"beta","double");
    stepsize=-1.0;
    m_param.add(stepsize,"stepsize","double");
    m_param.add(exp0,"exp0","double");
    m_param.add(MinMethod,"MinMethod","string");
    m_param.add(GEVtype,"GEVMethod","string");
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
    if (omp_get_max_threads()<3)
    {
      app_log()<<" Requires openMP and more than 3 threads!"<<endl;
      APP_ABORT("QMCCSLinearOptimize::start");
    }

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
    app_log()<<"  Line Minimization method "<<MinMethod<<endl;

}

bool QMCCSLinearOptimize::run()
{
    start();
    bool Valid(true);
    int Total_iterations(0);

//size of matrix
    numParams = optTarget->NumParams();
    N = numParams + 1;

app_log()<<"Here 1" <<endl; cout.flush();

//  solve CSFs and other parameters separately then rescale elements accordingly
    int first,last;
    getNonLinearRange(first,last);

app_log()<<"Here 2" <<endl; cout.flush();

//     initialize our parameters
    vector<RealType> currentParameterDirections(N,0);
    vector<RealType> currentParameters(numParams,0);
    for (int i=0; i<numParams; i++) currentParameters[i] = optTarget->Params(i);
    optdir.resize(numParams,0);
    optparm.resize(numParams,0);

    
    Matrix<RealType> Ham(N,N);
    Matrix<RealType> Ham2(N,N);
    Matrix<RealType> Var(N,N);
    Matrix<RealType> S(N,N);
    vmcEngine->fillMatrices(Ham2,Ham,Var,S);

app_log()<<"Here 3" <<endl; cout.flush();

    vector<RealType> bestParameters(currentParameters);


//this is the small amount added to the diagonal to stabilize the eigenvalue equation. 10^stabilityBase
    RealType stabilityBase(exp0);
    int tooManyTries(200);

    Matrix<RealType> Left(N,N);
    Matrix<RealType> Right(N,N);
    
    vector<std::pair<RealType,RealType> > mappedStabilizers;
    vector<vector<RealType> > savedCSparameters;
          
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
    RealType XS(0);

app_log()<<"Here 4" <<endl; cout.flush();
    
    if (nstabilizers<=omp_get_max_threads()+1) nstabilizers=omp_get_max_threads()+1;
    else
    {
      nstabilizers -= omp_get_max_threads()+1;
      int Ns(nstabilizers/(omp_get_max_threads()-2));
      nstabilizers = omp_get_max_threads() + 1 + Ns*(omp_get_max_threads()-2);
    }
    
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


        
        if(savedCSparameters.size()==omp_get_max_threads())
        {
          app_log()<<"   Choosing best"<<endl;
          RealType error(0);
          int bestP = vmcEngine->runCS(savedCSparameters,error);
          if (bestP<0)
          {
            app_log()<<"   Error in CS cost function. Unchanged parameters."<<endl;
            bestP=0;
            for (int i=0; i<numParams; i++) optTarget->Params(i) = currentParameters[i];
            finish();
            return false;
          }
          
          vector<RealType> cs(numParams);
          for (int i=0; i<numParams; i++) cs[i] = bestParameters[i] = savedCSparameters[bestP][i];
          savedCSparameters.clear();
          savedCSparameters.push_back(cs);
          std::pair<RealType,RealType> ms;
          ms.first=mappedStabilizers[bestP].first;
          ms.second=mappedStabilizers[bestP].second;
          
          if (stability==nstabilizers-1) continue;
//           if (MinMethod=="rescale")
//           {
            int nms=mappedStabilizers.size();
            std::vector<RealType> csts(nms);
            vmcEngine->getDeltaCosts(csts);
            
            if (nms>=5)
            {//Quartic fit the stabilizers we have tried and try to choose the best we can
              vector<RealType>  Y(nms), Coefs(5);
              Matrix<RealType> X(nms,5);
              for (int i=0; i<nms; i++) X(i,0)=1.0;
              for (int i=0; i<nms; i++) X(i,1)=mappedStabilizers[i].second;
              for (int i=0; i<nms; i++) X(i,2)=X(i,1)*X(i,1);
              for (int i=0; i<nms; i++) X(i,3)=X(i,2)*X(i,1);
              for (int i=0; i<nms; i++) X(i,4)=X(i,3)*X(i,1);
              for (int i=0; i<nms; i++) Y[i]=csts[i];
              LinearFit(Y,X,Coefs);
  //lowest we will allow is a little less than the bare base stabilizer
              RealType dltaBest=std::max(stabilityBase-0.1, QuarticMinimum(Coefs));
              XS = dltaBest;
            }
            else
            {//Quadratic fit the stabilizers we have tried and try to choose the best we can
              vector<RealType>  Y(nms), Coefs(3);
              Matrix<RealType> X(nms,3);
              for (int i=0; i<nms; i++) X(i,0)=1.0;
              for (int i=0; i<nms; i++) X(i,1)=mappedStabilizers[i].second;
              for (int i=0; i<nms; i++) X(i,2)=X(i,1)*X(i,1);
              for (int i=0; i<nms; i++) Y[i]=csts[i];
              LinearFit(Y,X,Coefs);
              
              RealType quadraticMinimum(-1.0*Coefs[1]/Coefs[2]);
              RealType dltaBest=std::max(stabilityBase-0.1, quadraticMinimum);
  //               app_log()<<"smallest XS:      "<<X(0,1)<<endl;
  //               app_log()<<"quadraticMinimum: "<<quadraticMinimum<<endl;
              XS = dltaBest;
            }
//           }
//           else
//             XS=0;
          mappedStabilizers.clear();
          mappedStabilizers.push_back(ms);
        }
        else
          XS=0;
            
        if (XS==0)
        {
            XS     = std::exp(stabilityBase +  stability*od_largest/nstabilizers);
            for (int i=1; i<N; i++) LeftT(i,i) += XS;
        }
        else
        {
          //else XS is from the fit
          for (int i=1; i<N; i++) LeftT(i,i) += std::exp(XS);
        }

        myTimers[2]->start();
        RealType lowestEV =getLowestEigenvector(LeftT,RightT,currentParameterDirections);
        myTimers[2]->stop();

//         RealType Lambda_Last(Lambda);
        myTimers[3]->start();
//         if (GEVtype=="H2")
//         {
// //           the rescaling isn't right for this?
//           RealType bigVec(0);
//           for (int i=0; i<numParams; i++) bigVec = std::max(bigVec,std::abs(currentParameterDirections[i+1]));
//           Lambda=0.5*bigChange/bigVec;
//         }
//           else 
        Lambda = getNonLinearRescale(currentParameterDirections,S);
        
          myTimers[3]->stop();
//         app_log()<<"Computed Lambda is: "<<Lambda<<endl;
        RealType bigVec(0);
        for (int i=0; i<numParams; i++) bigVec = std::max(bigVec,std::abs(currentParameterDirections[i+1]));
        if (Lambda*bigVec>bigChange)
        {
            app_log()<<"  Failed Step. Largest parameter change: "<<Lambda*bigVec<<endl;
            tooManyTries--;
            if (tooManyTries>0)
            {
                if(stability==0) stabilityBase+=stabilizerScale;
                else stabilityBase-=stabilizerScale;
                stability-=1;
                app_log()<<" Re-run with larger stabilityBase"<<endl;
                continue;
            }
        }

        if (MinMethod=="rescale")
        {
             vector<RealType> cs(numParams);
             for (int i=0; i<numParams; i++) cs[i] = currentParameters[i] + Lambda*currentParameterDirections[i+1];
             app_log()<<" Lambda: "<<Lambda<<endl;
             app_log()<<" Largest parameter change: "<<Lambda*bigVec<<endl;
//              optTarget->resetPsi(false);
             savedCSparameters.push_back(cs);
        }
        else
        {
             int nthreads = omp_get_max_threads();
             std::vector<RealType> lambdas(nthreads);
             if (stepsize>0) Lambda = stepsize/bigVec;
             for(int i=0;i<nthreads;i++) lambdas[i] = i/(nthreads-1.0)*Lambda;
             
             vmcEngine->runCS(currentParameters,currentParameterDirections,lambdas);
             Lambda=lambdas[0];
             if (Lambda==0)
             {
               app_log()<<"  Failed Step. Lambda=0."<<endl;
               tooManyTries--;
               if (tooManyTries>0)
               {
                if(stability==0) stabilityBase+=stabilizerScale;
                else stabilityBase-=stabilizerScale;
                   stability-=1;
//                    app_log()<<" Re-run with larger stabilityBase"<<endl;
                   continue;
               }
             }
             vector<RealType> cs(numParams);
             for (int i=0; i<numParams; i++) cs[i] = currentParameters[i] + Lambda*currentParameterDirections[i+1];
//              optTarget->resetPsi(false);
             savedCSparameters.push_back(cs);
        }
        std::pair<RealType,RealType> ms;
        ms.first=stability;
// the log fit seems to work best
        ms.second=std::log(XS);
        mappedStabilizers.push_back(ms);

    }
    
    for (int i=0; i<numParams; i++) optTarget->Params(i) = bestParameters[i];
//     optTarget->resetPsi(false);
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
        vmcEngine = new VMCLinearOptOMP(W,Psi,H,hamPool,psipool);
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

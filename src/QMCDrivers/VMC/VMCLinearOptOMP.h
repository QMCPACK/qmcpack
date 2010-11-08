//////////////////////////////////////////////////////////////////
// (c) Copyright 2003- by Jeongnim Kim
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
#ifndef QMCPLUSPLUS_VMCLINEAROPT_OMP_H
#define QMCPLUSPLUS_VMCLINEAROPT_OMP_H
#include "QMCDrivers/QMCDriver.h"
#include "QMCDrivers/CloneManager.h"
#include "Message/CommOperators.h"
namespace qmcplusplus
{

/** @ingroup QMCDrivers  ParticleByParticle
 * @brief Implements a VMC using particle-by-particle move. Threaded execution.
 */
class VMCLinearOptOMP: public QMCDriver, public CloneManager
{
public:
    /// Constructor.
    VMCLinearOptOMP(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h,
                    HamiltonianPool& hpool);
    bool run() {
        run(true);
    }
    bool run(bool needMatrix);
    RealType runCS(vector<RealType> curParams, vector<RealType> curDir, vector<RealType>& lambdas);
    int runCS(vector<vector<RealType> > bestParams, RealType& errorbars);
    bool put(xmlNodePtr cur);

    inline RealType Cost()
    {
//       pure energy minimization for line min
        app_log()<<"Energy: "<<E_avg<<endl;
        app_log()<<"Variance: "<<V_avg<<endl;
        return E_avg;
    }

    void fillMatrices(Matrix<RealType>& H2, Matrix<RealType>& Hamiltonian, Matrix<RealType>& Variance, Matrix<RealType>& Overlap);

    void clearComponentMatrices()
    {
        HDiHDj=0;
        DiHDj=0;
        DiHDjE=0;
        DiDj=0;
        DiDjE=0;
        DiDjE2=0;
        for (int i=0; i<NumOptimizables; i++)
        {
            HDi[i]=0;
            HDiE[i]=0;
            Di[i]=0;
            DiE[i]=0;
            DiE2[i]=0;
        }
        DerivRecords=0;
        HDerivRecords=0;
        sE=0;
        sE2=0;
        sE4=0;
        sW=0;
    }


private:
    ///number of warmup steps
    int myWarmupSteps;
    ///number of RN warmup steps
    int myRNWarmupSteps;
    ///option to enable/disable drift equation for VMC
    string UseDrift;
    ///target errorbars to use to determine when to stop for filling matrix and line minimization
    RealType alpha_errorbars, beta_errorbars;

    ///check the run-time environments
    void resetRun();
    ///copy constructor
    VMCLinearOptOMP(const VMCLinearOptOMP& a): QMCDriver(a),CloneManager(a) { }
    /// Copy operator (disabled).
    VMCLinearOptOMP& operator=(const VMCLinearOptOMP&)
    {
        return *this;
    }
    ///Ways to set rn constant
    RealType logoffset,logepsilon;
      
    int NumOptimizables;
    RealType w_beta;
    RealType E_avg, V_avg;
    string GEVtype;

    ///These are the values we collect to build the Matrices GLOBAL
    Matrix<RealType> HDiHDj, DiHDj, DiHDjE, DiDj, DiDjE, DiDjE2;
    std::vector<RealType> HDi, HDiE, Di, DiE, DiE2;
    RealType sE,sE2,sE4,sW;
    ///Temp matrices
    Matrix<RealType> DerivRecords, HDerivRecords;


    void initCS();
    
    void resizeForOpt(int n)
    {
        HDiHDj.resize(n,n);
        DiHDj.resize(n,n);
        DiHDjE.resize(n,n);
        DiDj.resize(n,n);
        DiDjE.resize(n,n);
        DiDjE2.resize(n,n);
        HDi.resize(n);
        HDiE.resize(n);
        Di.resize(n);
        DiE.resize(n);
        DiE2.resize(n);
        DerivRecords.resize(NumThreads,n);
        HDerivRecords.resize(NumThreads,n);
    }

    RealType fillComponentMatrices(bool needMatrix)
    {
        int n(NumOptimizables);
        ///These are the values we collect to build the Matrices LOCAL
        Matrix<RealType> lHDiHDj(n,n), lDiHDj(n,n), lDiHDjE(n,n), lDiDj(n,n), lDiDjE(n,n), lDiDjE2(n,n);
        std::vector<RealType> lHDi(n), lHDiE(n), lDi(n), lDiE(n), lDiE2(n);
        RealType lsE,lsE2,lsE4,lsW;

        for (int ip=0; ip<NumThreads; ip++)
        {
            MCWalkerConfiguration::iterator wit(W.begin()+wPerNode[ip]), wit_end(W.begin()+wPerNode[ip+1]);
            RealType E_L = (*wit)->getPropertyBase()[LOCALENERGY];
            RealType E_L2= E_L*E_L;
            RealType wW  = (*wit)->Weight;
            lsE +=E_L*wW;
            lsE2+=E_L2*wW;
            lsE4+=E_L2*E_L2*wW;
            lsW +=wW;
        }

        if (needMatrix)
        {
            for (int ip=0; ip<NumThreads; ip++)
            {
                MCWalkerConfiguration::iterator wit(W.begin()+wPerNode[ip]), wit_end(W.begin()+wPerNode[ip+1]);
                RealType E_L = (*wit)->getPropertyBase()[LOCALENERGY];
                RealType E_L2= E_L*E_L;
                RealType wW  = (*wit)->Weight;
                for (int i=0; i<NumOptimizables; i++)
                {
                    RealType di  = DerivRecords(ip,i);
                    RealType hdi = HDerivRecords(ip,i);
                    //             vectors
                    lHDiE[i]+= wW*E_L* hdi;
                    lHDi[i] += wW*     hdi;
                    lDiE2[i]+= wW*E_L2*di;
                    lDiE[i] += wW*E_L* di;
                    lDi[i]  += wW*     di;

                    for (int j=0; j<NumOptimizables; j++)
                    {
                        RealType dj  = DerivRecords(ip,j);
                        RealType hdj = HDerivRecords(ip,j);

                        lHDiHDj(i,j) += wW*    hdi*hdj;
                        lDiHDjE(i,j) += wW* E_L*di*hdj;
                        lDiHDj(i,j)  += wW*     di*hdj;
                        lDiDjE2(i,j) += wW*E_L2*di*dj;
                        lDiDjE(i,j)  += wW* E_L*di*dj;
                        lDiDj(i,j)   += wW*     di*dj;
                    }
                }
            }
        }

        //Lazy. Pack these for better performance.
        myComm->allreduce(lsE);
        myComm->allreduce(lsE2);
        myComm->allreduce(lsE4);
        myComm->allreduce(lsW);
        if (needMatrix)
        {
            myComm->allreduce(lHDiE);
            myComm->allreduce(lHDi);
            myComm->allreduce(lDiE2);
            myComm->allreduce(lDiE);
            myComm->allreduce(lDi);
            myComm->allreduce(lHDiHDj);
            myComm->allreduce(lDiHDjE);
            myComm->allreduce(lDiHDj);
            myComm->allreduce(lDiDjE2);
            myComm->allreduce(lDiDjE);
            myComm->allreduce(lDiDj);
        }
        //add locals to globals
        sE +=lsE;
        sE2+=lsE2;
        sE4+=lsE4;
        sW +=lsW;
        if (needMatrix)
        {
            for (int j=0; j<NumOptimizables; j++)
            {
                HDiE[j]+=lHDiE[j];
                HDi[j] +=lHDi[j] ;
                DiE2[j]+=lDiE2[j];
                DiE[j] +=lDiE[j] ;
                Di[j]  +=lDi[j]  ;
            }

            HDiHDj += lHDiHDj;
            DiHDjE += lDiHDjE;
            DiHDj  += lDiHDj ;
            DiDjE2 += lDiDjE2;
            DiDjE  += lDiDjE ;
            DiDj   += lDiDj  ;
        }

        RealType nrm = 1.0/sW;
        E_avg = nrm*sE;
        V_avg = nrm*sE2-E_avg*E_avg;
        
        RealType err_E(std::sqrt(V_avg*nrm));
        RealType err_E2(nrm*sE4-nrm*nrm*sE2*sE2);
        err_E2 *= nrm;
        err_E2 = std::sqrt(err_E2);

        return w_beta*err_E2+(1.0-w_beta)*err_E;
    }
    
    Matrix<long double> CorrelatedH, Norm2s;
    vector<long double> Norms;
    vector<long double> Energies;
    vector<long double> NE_i;
    int CSBlock;
    int minE,nE;
    RealType logpsi2_0_0;
    RealType estimateCS();
    
    bool moved_right;
    bool moved_left;
    bool bracketing(vector<RealType>& lambdas, RealType errorbars);
    
};

}
#endif
/***************************************************************************
 * $RCSfile: VMCLinearOptOMP.h,v $   $Author: jnkim $
 * $Revision: 1.5 $   $Date: 2006/07/17 14:29:40 $
 * $Id: VMCLinearOptOMP.h,v 1.5 2006/07/17 14:29:40 jnkim Exp $
 ***************************************************************************/

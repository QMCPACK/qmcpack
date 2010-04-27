//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef QMCPLUSPLUS_BACKFLOW_TRANSFORMATION_H
#define QMCPLUSPLUS_BACKFLOW_TRANSFORMATION_H
#include "Utilities/ProgressReportEngine.h"
#include "OhmmsData/AttributeSet.h"
#include "QMCWaveFunctions/OrbitalSetTraits.h"
#include "QMCWaveFunctions/Fermion/BackflowFunctionBase.h"
#include "QMCWaveFunctions/Fermion/Backflow_ee.h"
#include "QMCWaveFunctions/Fermion/Backflow_eI.h"
#include "QMCWaveFunctions/Fermion/GaussianFunctor.h"
#include "QMCWaveFunctions/Jastrow/BsplineFunctor.h"
#include "Particle/ParticleSet.h"
#include "Configuration.h"
#include <map>
#include <cmath>
#include "OhmmsPETE/OhmmsArray.h"

namespace qmcplusplus
{

  class BackflowTransformation: public OrbitalSetTraits<QMCTraits::ValueType> 
  {

    public:

    typedef map<string,ParticleSet*>   PtclPoolType;
    typedef Array<HessType,3>       HessArray_t;
    //typedef Array<GradType,3>       GradArray_t;
    //typedef Array<PosType,3>        PosArray_t;

    ///number of quantum particles
    int NumTargets;

    // quasiparticle coordinates
    ParticleSet QP; 

    // number of variational parameters
    int numParams;

    // map index of variables from local arrays to outside world
    map<int,int> optIndexMap;

    // pos of first optimizable variable in global array
    int numVarBefore; 

    ParticleSet& targetPtcl; 
    PtclPoolType& ptclPool;

    // matrix of laplacians
    // /vec{B(i)} = sum_{k} /grad_{k}^2 /vec{x_i} 
    GradVector_t Bmat;

    GradMatrix_t Bmat_full;

    // matrix of first derivatives 
    // A(i,j)[a,b] = (Grad_i)_a (x_j)_b  
    //               i,j:particle index
    //               a,b=(x,y,z)  
// notice that A(i,j) is a symmetric matrix, improve later
    HessMatrix_t Amat;


    // \nabla_a A_{i,j}^{\alpha,\beta}
    // derivative of A matrix with respect to var. prms.
    HessArray_t Xmat;

    // \sum_i \nabla_a B_{i,j}^{\alpha}
    GradMatrix_t Ymat;

    // \nabla_a x_i^{\alpha}   
    GradMatrix_t Cmat;

    // Identity
    HessType HESS_ID;
    HessType DummyHess;

    vector<BackflowFunctionBase*> bfFuns;     

    map<string,int> sources;     
    vector<string> names;     

    BackflowTransformation(ParticleSet& els, PtclPoolType& pool):
      targetPtcl(els),QP(els),ptclPool(pool) {
      NumTargets=els.getTotalNum();
      Bmat.resize(NumTargets);
      Bmat_full.resize(NumTargets,NumTargets);
      Amat.resize(NumTargets,NumTargets);
      HESS_ID.diagonal(1.0);
      DummyHess=0.0;
      numVarBefore=0;
    }

    BackflowTransformation(BackflowTransformation &tr): 
      NumTargets(tr.NumTargets),QP(tr.QP), 
      targetPtcl(tr.targetPtcl),ptclPool(tr.ptclPool) { 
      Bmat.resize(NumTargets);
      Bmat_full.resize(NumTargets,NumTargets);
      Amat.resize(NumTargets,NumTargets);
      HESS_ID.diagonal(1.0);
      DummyHess=0.0;
      numVarBefore=tr.numVarBefore;
      optIndexMap=tr.optIndexMap; 
      bfFuns.resize((tr.bfFuns).size());
      vector<BackflowFunctionBase*>::iterator it((tr.bfFuns).begin());
      for(int i=0; i<(tr.bfFuns).size() ; i++,it++)
        bfFuns[i] = (*it)->makeClone();
    }
    

    BackflowTransformation* makeClone()
    {
       BackflowTransformation *clone = new BackflowTransformation(*this);
       return clone; 
    }

    ~BackflowTransformation() {}; 

    void checkInVariables(opt_variables_type& active)
    {
      for(int i=0; i<bfFuns.size(); i++) bfFuns[i]->checkInVariables(active);
      numVarBefore = bfFuns[0]->indexOffset();
    }

    void checkOutVariables(const opt_variables_type& active)
    {
      for(int i=0; i<bfFuns.size(); i++) bfFuns[i]->checkOutVariables(active);
    }

    bool put(xmlNodePtr cur)
    {
      bool first=true;
      bool success=true;
      xmlNodePtr curRoot=cur;
      string cname;      

      cur = curRoot->children;
      while (cur != NULL)
      {
        getNodeName(cname,cur);
        if (cname == "transformation")
        {
          OhmmsAttributeSet spoAttrib;
          string source("none");
          string name("bf0");
          string type("none");
          RealType cusp=0.0;
          string funct("Gaussian");
          string unique("no");
          spoAttrib.add (name, "name");
          spoAttrib.add (type, "type");
          spoAttrib.add (cusp, "cusp");
          spoAttrib.add (source, "source");
          spoAttrib.add (funct, "function");
          spoAttrib.add (unique, "unique");
          spoAttrib.put(cur);

          sources[source] = names.size();
          names.push_back(name);
          if(type == "e-e") {
            app_log() <<"Adding electron-electron backflow. \n";

            if(funct == "Gaussian") {
   APP_ABORT("Disabled GaussianFunctor for now, \n");
              app_log() <<"Using GaussianFunctor type. \n";
              GaussianFunctor *GaussFT = new GaussianFunctor();
              GaussFT->put(cur);
              BackflowFunctionBase *tbf = (BackflowFunctionBase *) new Backflow_ee<GaussianFunctor>(targetPtcl,targetPtcl,GaussFT);
              bfFuns.push_back(tbf);
            } else if(funct == "Bspline")  {
              app_log() <<"Using BsplineFunctor type. \n";
              BsplineFunctor<double> *bsp = new BsplineFunctor<double>(cusp);
              bsp->put(cur);
              BackflowFunctionBase *tbf = (BackflowFunctionBase *) new Backflow_ee<BsplineFunctor<double> >(targetPtcl,targetPtcl,bsp);
              tbf->numParams = bsp->NumParams;
              tbf->derivs.resize(tbf->numParams);
              bfFuns.push_back(tbf);
            } else {
              APP_ABORT("Unknown function type in e-e BF Transformation.\n");
            }
          } else if(type == "e-e-I") {
            APP_ABORT("e-e-I backflow is not implemented yet. \n");
            app_log() <<"Adding electron-electron-Ion backflow. \n";
          } else if(type == "e-I") {
            app_log() <<"Adding electron-Ion backflow for source:"
                      <<source <<" \n";

            ParticleSet* ions=0;
            PtclPoolType::iterator pit(ptclPool.find(source));
            if(pit == ptclPool.end())
            {
              APP_ABORT("Missing backflow/@source.");
            } else {
              ions=(*pit).second;
            }

            if(funct == "Gaussian") {
   APP_ABORT("Disabled GaussianFunctor for now, \n");
              app_log() <<"Using GaussianFunctor type. \n";
              if(unique == "yes") {
                APP_ABORT("Unique radial functors for e-I backflow not implemented. \n"); 
              } else {
                app_log() <<"Using a single radial functor. \n";
                GaussianFunctor *GaussFT = new GaussianFunctor();
                GaussFT->put(cur);
                BackflowFunctionBase *tbf = (BackflowFunctionBase *) new Backflow_eI<GaussianFunctor>(*ions,targetPtcl,GaussFT);
                bfFuns.push_back(tbf);
              }
            } else if(funct == "Bspline")  {
              app_log() <<"Using BsplineFunctor type. \n";
              if(unique == "yes") {
                //APP_ABORT("Unique radial functors for e-I backflow not implemented. \n");
                app_log() <<"Using unique radial functors for ions of this type. \n";
                BackflowFunctionBase *tbf = (BackflowFunctionBase *) new Backflow_eI<BsplineFunctor<double> >(*ions,targetPtcl);
                Backflow_eI<BsplineFunctor<double> > *dum = (Backflow_eI<BsplineFunctor<double> >*) tbf;
                tbf->uniqueFunctions=true;
                int nI = ions->getTotalNum();
                tbf->numParams=0;
                for(int i=0; i<nI; i++)  {
                  BsplineFunctor<double> *bsp = new BsplineFunctor<double>(cusp);
// right now doesn't work because I need to loop over different xml nodes
                  bsp->put(cur); 
                  dum->RadFun.push_back(bsp);
                  tbf->numParams += bsp->NumParams;
                }
                tbf->derivs.resize(tbf->numParams);
                bfFuns.push_back(tbf);
              } else {
                app_log() <<"Using a single radial functor for all ions of this type. \n";
                BsplineFunctor<double> *bsp = new BsplineFunctor<double>(cusp);
                bsp->put(cur);
                BackflowFunctionBase *tbf = (BackflowFunctionBase *) new Backflow_eI<BsplineFunctor<double> >(*ions,targetPtcl,bsp);
                tbf->numParams = bsp->NumParams;
                tbf->derivs.resize(tbf->numParams); 
                bfFuns.push_back(tbf);
              }
            } else {
              APP_ABORT("Unknown function type in e-I BF Transformation.\n");
            }

          } else {
            APP_ABORT("Unknown backflow type. \n");
          }
          cur = cur->next;
        }
        cur = cur->next;
      }
      return success;
    }

    /** reset the distance table with a new target P
     */
    void resetTargetParticleSet(ParticleSet& P)
    {
      for(int i=0; i<bfFuns.size(); i++)
        bfFuns[i]->resetTargetParticleSet(P);
    }


    void resetParameters(const opt_variables_type& active)
    {
      //reset each unique basis functions
      for(int i=0; i<bfFuns.size(); i++)
        bfFuns[i]->resetParameters(active);
    }    

    
    /** calculate quasi-particle coordinates only
     */
    inline void 
    transformOnly(const ParticleSet& P)  
    {
      for(int i=0; i<NumTargets; i++) QP.R[i] = P.R[i];
      for(int i=0; i<bfFuns.size(); i++) bfFuns[i]->evaluate(P,QP);
      QP.update(0);  // update distance tables
    }

    
    /** calculate quasi-particle coordinates, Bmat and Amat 
     */
    inline void
    evaluate(const ParticleSet& P) 
    {
      Bmat=0.0;
      Amat=0.0;
      Bmat_full=0.0;
      for(int i=0; i<NumTargets; i++) {
        QP.R[i] = P.R[i];
        Amat(i,i).diagonal(1.0);
      } 
      for(int i=0; i<bfFuns.size(); i++) bfFuns[i]->evaluate(P,QP,Bmat_full,Amat);
      QP.update(0);  // update distance tables
    } 

    inline void 
    evaluateDerivatives(const ParticleSet& P)
    {

      if(Cmat.size() == 0) { // initialize in the first call
        // assumes that all BF parameters are packed together in  
        // active variable set. is this always correct???
        numParams=0;
        for(int i=0; i<bfFuns.size(); i++) { 
          int tmp = bfFuns[i]->setParamIndex(numParams);
          numParams+=tmp; 
        } 

        for(int i=0; i<numParams; i++) 
          optIndexMap[i] = i+numVarBefore;

        Cmat.resize(numParams,NumTargets);
        Xmat.resize(numParams,NumTargets,NumTargets);
        Ymat.resize(numParams,NumTargets);
      }
    
      // Uncomment to test calculation of Cmat,Xmat,Ymat
      //testDeriv(P); 

      Bmat=0.0;
      Amat=0.0;
      Bmat_full=0.0;
      Cmat=0.0;
      Ymat=0.0;
      Xmat=DummyHess;
      for(int i=0; i<NumTargets; i++) {
        QP.R[i] = P.R[i];
        Amat(i,i).diagonal(1.0);
      }
      for(int i=0; i<bfFuns.size(); i++) bfFuns[i]->evaluateWithDerivatives(P,QP,Bmat_full,Amat,Cmat,Ymat,Xmat);
      QP.update(0);
    }

    void testDeriv(const ParticleSet& P)
    {

       if(Cmat.size() == 0) { // initialize in the first call
         Cmat.resize(numParams,NumTargets);
         Xmat.resize(numParams,NumTargets,NumTargets);
         Ymat.resize(numParams,NumTargets);
       }

       Bmat=0.0;
       Amat=0.0;
       Bmat_full=0.0;
       Cmat=0.0;
       Ymat=0.0;
       Xmat=DummyHess;
       for(int i=0; i<NumTargets; i++) {
         QP.R[i] = P.R[i];
         Amat(i,i).diagonal(1.0);
       }
       for(int i=0; i<bfFuns.size(); i++) bfFuns[i]->evaluateWithDerivatives(P,QP,Bmat_full,Amat,Cmat,Ymat,Xmat);

       ParticleSet::ParticlePos_t qp_0;
       ParticleSet::ParticlePos_t qp_1;
       ParticleSet::ParticlePos_t qp_2;
       GradMatrix_t Bmat_full_1;
       HessMatrix_t Amat_1;
       GradMatrix_t Bmat_full_2;
       HessMatrix_t Amat_2;
       double dh = 0.00001;

       qp_0.resize(NumTargets);
       qp_1.resize(NumTargets);
       qp_2.resize(NumTargets);
       Bmat_full_1.resize(NumTargets,NumTargets); 
       Bmat_full_2.resize(NumTargets,NumTargets); 
       Amat_1.resize(NumTargets,NumTargets); 
       Amat_2.resize(NumTargets,NumTargets); 

       for(int i=0; i<NumTargets; i++) {
          qp_0[i] = QP.R[i];
       }

       app_log() <<" Testing derivatives of backflow transformation. \n";
       app_log() <<" Numtargets: " <<NumTargets <<endl;
       opt_variables_type wfVars,wfvar_prime;
       checkInVariables(wfVars);
       checkOutVariables(wfVars);
       int Nvars= wfVars.size();
       wfvar_prime= wfVars;
       wfVars.print(cout);

       for(int i=0; i<Nvars; i++) {

         for (int j=0; j<Nvars; j++) wfvar_prime[j]=wfVars[j];
         wfvar_prime[i] = wfVars[i]+ dh;
         resetParameters(wfvar_prime);

         Bmat_full_1=0.0;
         Amat_1=0.0;
         for(int k=0; k<NumTargets; k++) {
           QP.R[k] = P.R[k];
           Amat_1(k,k).diagonal(1.0);
         }
         for(int k=0; k<bfFuns.size(); k++) bfFuns[k]->evaluate(P,QP,Bmat_full_1,Amat_1);
         for(int k=0; k<NumTargets; k++) qp_1[k] = QP.R[k];


         for (int j=0; j<Nvars; j++) wfvar_prime[j]=wfVars[j];
         wfvar_prime[i] = wfVars[i] - dh;
         resetParameters(wfvar_prime);

         Bmat_full_2=0.0;
         Amat_2=0.0;
         for(int k=0; k<NumTargets; k++) {
           QP.R[k] = P.R[k];
           Amat_2(k,k).diagonal(1.0);
         }
         for(int k=0; k<bfFuns.size(); k++) bfFuns[k]->evaluate(P,QP,Bmat_full_2,Amat_2);
         for(int k=0; k<NumTargets; k++) qp_2[k] = QP.R[k];

         app_log() <<"Cmat: \n" 
                   <<"i, AvDiff, max: \n";
         ValueType df,av=0.0,cnt=0.0,maxD=-100.0;
         for(int k=0; k<NumTargets; k++) { 
          for(int q=0; q<3; q++) {
           cnt++;
           df=(( (qp_1[k])[q] - (qp_2[k])[q] )/(2.0*dh)-Cmat(i,k)[q]);
           av+=df;
           if( std::fabs(df) > maxD ) maxD=std::fabs(df); 
           //app_log() <<k <<"  " <<q <<"   "
           //          <<( (qp_1[k])[q] - (qp_2[k])[0] )/(2.0*dh)   <<"  "
           //          <<Cmat(i,k)[q] <<"  " <<(( (qp_1[k])[q] - (qp_2[k])[q] )/(2.0*dh)-Cmat(i,k)[q]) <<endl;
          } 
         } 
         app_log() <<i <<"  " <<av/cnt <<"  " <<maxD <<endl;
         av=cnt=maxD=0.0;

         app_log() <<"Ymat: \n";
         for(int k=0; k<NumTargets; k++) {
           for(int q=0; q<3; q++) {
             RealType dB=0.0;
             for(int j=0; j<NumTargets; j++) dB+= (Bmat_full_1(j,k)[q] - Bmat_full_2(j,k)[q]);
             cnt++;
             df=(dB/(2.0*dh)-Ymat(i,k)[q]);
             av+=df;
             if( std::fabs(df) > maxD ) maxD=std::fabs(df); 
             //app_log() <<k <<"  " <<q <<"   "
             //        <<dB/(2.0*dh)   <<"  "
             //        <<Ymat(i,k)[q] <<"  " <<(dB/(2.0*dh)-Ymat(i,k)[q]) <<endl;
           }         
         }         
         app_log() <<i <<"  " <<av/cnt <<"  " <<maxD <<endl;
         av=cnt=maxD=0.0;

         app_log() <<"Xmat: \n";
         for(int k1=0; k1<NumTargets; k1++) 
          for(int k2=0; k2<NumTargets; k2++) {
           for(int q1=0; q1<3; q1++) {
            for(int q2=0; q2<3; q2++) {
             RealType dB=(Amat_1(k1,k2))(q1,q2) - (Amat_2(k1,k2))(q1,q2);
             cnt++;
             df=(dB/(2.0*dh)-(Xmat(i,k1,k2))(q1,q2));
             av+=df;
             if( std::fabs(df) > maxD ) maxD=std::fabs(df); 
             //app_log() <<k1 <<"  " <<k2 <<"  " <<q1 <<"  " <<q2 <<"   "
             //        <<dB/(2.0*dh)   <<"  "
             //        <<(Xmat(i,k1,k2))(q1,q2) <<"  " <<(dB/(2.0*dh)-(Xmat(i,k1,k2))(q1,q2)) <<endl;
             }
            }
           }
         app_log() <<i <<"  " <<av/cnt <<"  " <<maxD <<endl;
         av=cnt=maxD=0.0;
 
       }


    }

  };

}

#endif

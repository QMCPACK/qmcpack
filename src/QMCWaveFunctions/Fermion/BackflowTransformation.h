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

namespace qmcplusplus
{

  class BackflowTransformation: public OrbitalSetTraits<QMCTraits::ValueType> 
  {

    public:

    typedef map<string,ParticleSet*> PtclPoolType;

    ///number of quantum particles
    int NumTargets;

    // quasiparticle coordinates
    ParticleSet QP; 

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

    // Identity
    HessType HESS_ID;

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
    }

    BackflowTransformation(BackflowTransformation &tr): 
      NumTargets(tr.NumTargets),QP(tr.QP), 
      targetPtcl(tr.targetPtcl),ptclPool(tr.ptclPool) { 
      Bmat.resize(NumTargets);
      Bmat_full.resize(NumTargets,NumTargets);
      Amat.resize(NumTargets,NumTargets);
      HESS_ID.diagonal(1.0);
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
    }

    void checkOutVariables(const opt_variables_type& active)
    {
      for(int i=0; i<bfFuns.size(); i++) bfFuns[i]->checkOutVariables(active);
    }

    bool put(xmlNodePtr cur)
    {

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

app_log() <<"#sources: " <<ptclPool.size() <<endl;
for(PtclPoolType::iterator pit=ptclPool.begin(); pit!=ptclPool.end(); pit++)
  app_log() << (*pit).first <<"   "   <<(*pit).second <<endl; 
app_log() <<"end - sources: " <<endl;

            ParticleSet* ions=0;
            PtclPoolType::iterator pit(ptclPool.find(source));
            if(pit == ptclPool.end())
            {
              APP_ABORT("Missing backflow/@source.");
            } else {
              ions=(*pit).second;
            }

            if(funct == "Gaussian") {
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
                APP_ABORT("Unique radial functors for e-I backflow not implemented. \n");
              } else {
                app_log() <<"Using a single radial functor for all ions of this type. \n";
                BsplineFunctor<double> *bsp = new BsplineFunctor<double>(cusp);
                bsp->put(cur);
                BackflowFunctionBase *tbf = (BackflowFunctionBase *) new Backflow_eI<BsplineFunctor<double> >(targetPtcl,targetPtcl,bsp);
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

  };

}

#endif

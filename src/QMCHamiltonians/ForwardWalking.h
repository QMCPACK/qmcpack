//////////////////////////////////////////////////////////////////
// (c) Copyright 2003  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
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
#ifndef QMCPLUSPLUS_FORWARDWALKING_H
#define QMCPLUSPLUS_FORWARDWALKING_H
#include "Particle/ParticleSet.h"
#include "Particle/WalkerSetRef.h"
#include "QMCHamiltonians/QMCHamiltonianBase.h"
#include "ParticleBase/ParticleAttribOps.h"
#include "OhmmsData/ParameterSet.h"
#include "OhmmsData/AttributeSet.h"

#include "QMCHamiltonians/QMCHamiltonian.h"

namespace qmcplusplus {


  struct ForwardWalking: public QMCHamiltonianBase {
    vector<int> Hindices;
    vector<int> Pindices;
    vector<vector<int> > walkerLengths;
    vector<double> Values;
    vector<string> Names;
    int blockT,nObservables,nValues,FirstHamiltonian;
    Walker<Return_t, ParticleSet::ParticleGradient_t>* tWalker;


    /** constructor
     *
     * Pressure operators need to be re-evaluated during optimization.
     */
    ForwardWalking() {
      UpdateMode.set(OPTIMIZABLE,1);
    }
    
    ForwardWalking(ParticleSet& P) {
      UpdateMode.set(OPTIMIZABLE,1);
    }
    
    ///destructor
    ~ForwardWalking() { }

    void resetTargetParticleSet(ParticleSet& P) {
    }
    
    void setHistories(Walker<Return_t, ParticleSet::ParticleGradient_t>& ThisWalker){
      (tWalker)= &(ThisWalker);
    }

    inline Return_t 
    evaluate(ParticleSet& P) {

      P.phLength++;
      if (P.phLength%blockT == 0){
        P.phLength=0;
        for(int i=0;i<nObservables;i++){
          tWalker->addPropertyHistoryPoint(Pindices[i],P.PropertyList[ Hindices[i] ]);
        }
        
        
        vector<double>::iterator Vit=Values.begin();
        for(int i=0;i<nObservables;i++){
          for(int j=0;j<walkerLengths[i].size();j++,Vit++){
            (*Vit) = tWalker->PropertyHistory[Pindices[i]][walkerLengths[i][j]-1];
            
          }
        }
      }
      return 0.0;
    }

    inline Return_t 
    evaluate(ParticleSet& P, vector<NonLocalData>& Txy) {
      return evaluate(P);
    }

    bool put(xmlNodePtr cur) {return true;}
      
    bool put(xmlNodePtr cur, QMCHamiltonian& h, ParticleSet& P ) {
      FirstHamiltonian = h.startIndex();
      nObservables=0;
      nValues=0;
      blockT=1;
      OhmmsAttributeSet attrib;
      attrib.add(blockT,"blockSize");
      attrib.put(cur);
      app_log()<<"  Forward walking block size is "<< blockT<<"*Tau"<<endl;
      P.phLength=0;

      xmlNodePtr tcur = cur->children;
      while(tcur != NULL) {
        string cname((const char*)tcur->name);
//         app_log()<<cname<<endl;
        if(cname == "Observable") 
        {
          string tagName("none");
          int Hindex(-1);
          OhmmsAttributeSet Tattrib;
          Tattrib.add(tagName,"name");
          Tattrib.put(tcur);
          
          Hindex=h.getObservable(tagName);
          
          if (Hindex<0){
            app_log()<<" Hamiltonian Element "<<tagName<<" does not exist!! "<<Hindex<<endl;
            assert(Hindex>=0);
          }
          Names.push_back(tagName);
          Hindices.push_back(FirstHamiltonian+Hindex);
          app_log()<<" Hamiltonian Element "<<tagName<<" was found at "<<FirstHamiltonian+Hindex<<endl;
          
          vector<int> Parameters;
          putContent(Parameters,tcur);
          int numT=Parameters.size();
          nObservables+=1;
          nValues+=numT;
          
          app_log()<<"   "<<numT<<" values will be calculated at block numbers:"<<endl;
          app_log()<<"      ";
          for(int nm=0;nm<Parameters.size();nm++) app_log()<<Parameters[nm]<<"  ";
          app_log()<<endl;
          walkerLengths.push_back(Parameters);
          
          int maxWsize=Parameters.back();
          Pindices.push_back(P.addPropertyHistory(maxWsize));
          
        }
        tcur = tcur->next;
      }
      app_log()<<"Total number of observables calculated:"<<nObservables<<endl;
      app_log()<<"Total number of values calculated:"<<nValues<<endl;
      for (int i=0;i<nValues;i++) Values.push_back(0.0);
      return true;
    }

    bool get(std::ostream& os) const {
      os << "ForwardWalking";
      return true;
    }

    QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi)
    {
      ForwardWalking* myClone = new ForwardWalking();
      ///Do I need to do this?
      myClone->Hindices=Hindices;
      myClone->Pindices=Pindices;
      myClone->walkerLengths=walkerLengths;
      myClone->blockT=blockT;
      myClone->Names=Names;
      myClone->nValues=nValues;
      myClone->nObservables=nObservables;
      myClone->myIndex=myIndex;
      myClone->FirstHamiltonian=FirstHamiltonian;
      myClone->Values=Values;
      return myClone;
    }

    void addObservables(PropertySetType& plist)
    {
      bool first(true);
      for(int i=0;i<nObservables;i++){
        for(int j=0;j<walkerLengths[i].size();j++){
          std::stringstream sstr;
          sstr << "FW_" << Names[i] << "_" << walkerLengths[i][j];
          app_log()<<"Observables named "<<sstr.str()<<endl;
          if (first) {
            myIndex = plist.add(sstr.str());
            first=false;
          } else {
            plist.add(sstr.str());
          }
        }
      }
    }


    void setObservables(PropertySetType& plist)
    {
      for (int i=0;i<nValues;i++) plist[myIndex+i]=Values[i];
    }
    
    void setParticleSet(PropertySetType& plist, int offset)
    {
      for (int i=0;i<nValues;i++) plist[myIndex+i+offset]=Values[i];
    }

  };
}
#endif

/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1581 $   $Date: 2007-01-04 10:02:14 -0600 (Thu, 04 Jan 2007) $
 * $Id: Pressure.h 1581 2007-01-04 16:02:14Z jnkim $ 
 ***************************************************************************/


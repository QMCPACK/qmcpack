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
#ifndef QMCPLUSPLUS_TRIALDMCCORRECTION_H
#define QMCPLUSPLUS_TRIALDMCCORRECTION_H
#include "Particle/ParticleSet.h"
#include "Particle/WalkerSetRef.h"
#include "QMCHamiltonians/QMCHamiltonianBase.h"
#include "ParticleBase/ParticleAttribOps.h"
#include "OhmmsData/ParameterSet.h"
#include "OhmmsData/AttributeSet.h"

#include "QMCHamiltonians/QMCHamiltonian.h"
#include <cstdlib>

namespace qmcplusplus {


  struct TrialDMCCorrection: public QMCHamiltonianBase {
    vector<int> Hindices;
    vector<int> Pindices;
    vector<vector<int> > walkerLengths;
    vector<double> Values,EValues;
    vector<string> Names;
    int blockT,nObservables,nValues,FirstHamiltonian;
    
    double count;


    /** constructor
     *
     * Pressure operators need to be re-evaluated during optimization.
     */
    TrialDMCCorrection() {
      UpdateMode.set(OPTIMIZABLE,1);
    }
    
    TrialDMCCorrection(ParticleSet& P) {
      UpdateMode.set(OPTIMIZABLE,1);
    }
    
    ///destructor
    ~TrialDMCCorrection() { }

    void resetTargetParticleSet(ParticleSet& P) {
    }
    


    inline Return_t 
    rejectedMove(ParticleSet& P) {
        vector<double>::iterator Vit=Values.begin();
	vector<double>::iterator Vit2=EValues.begin();

        for(int i=0;i<nObservables;i++){
//           app_log()<<"Obs#"<<i;
          (*Vit)=0.0;
	  int k=0;
	  for(int j=0;j<walkerLengths[i].back();j++){
            (*Vit) += tWalker->PropertyHistory[Pindices[i]][j];
	    if(j==walkerLengths[i][k]-1){
	      double Tsum=(*Vit);
	      Vit++;
	      (*Vit)=Tsum;
	      k++;
	      (*Vit2)=Tsum* P.PropertyList[LOCALENERGY];
	      Vit2++;
	    }
//             app_log()<<"  "<<tWalker->PropertyHistory[Pindices[i]][walkerLengths[i][j]-1];
          }
// 	  app_log()<<endl;
        }
//       }
// 	double* wFWval = tWalker->getPropertyBase();
	std::copy(Values.begin(),Values.end(),tWalker->getPropertyBase()+FirstHamiltonian+myIndex );

      return 0.0;
    }
    
    inline Return_t 
    evaluate(ParticleSet& P) {

//       tWalker->phLength++;
//       count=tWalker->phLength;
//       if (tWalker->phLength%blockT == 0){
//         tWalker->phLength=0;
        for(int i=0;i<nObservables;i++){
          tWalker->addPropertyHistoryPoint(Pindices[i],  P.PropertyList[Hindices[i]]);
        }
	
//         for(int i=0;i<nObservables;i++){
// 	app_log()<<" Nobs:"<<i<<endl;
// 	for(int j=0;j<tWalker->PropertyHistory[i].size();j++){
//             app_log()<<"  "<<tWalker->PropertyHistory[i][j];
// 	  }
// 	  app_log()<<endl;
//         }
        
        vector<double>::iterator Vit=Values.begin();
	vector<double>::iterator Vit2=EValues.begin();

        for(int i=0;i<nObservables;i++){
//           app_log()<<"Obs#"<<i;
          (*Vit)=0.0;
	  int k=0;
	  for(int j=0;j<walkerLengths[i].back();j++){
            (*Vit) += tWalker->PropertyHistory[Pindices[i]][j];
	    if(j==walkerLengths[i][k]-1){
	      double Tsum=(*Vit);
	      Vit++;
	      (*Vit)=Tsum;
	      k++;
	      (*Vit2)=Tsum* P.PropertyList[LOCALENERGY];
	      Vit2++;
	    }
//             app_log()<<"  "<<tWalker->PropertyHistory[Pindices[i]][walkerLengths[i][j]-1];
          }
// 	  app_log()<<endl;
        }
//       }
// 	double* wFWval = tWalker->getPropertyBase();
	std::copy(Values.begin(),Values.end(),tWalker->getPropertyBase()+FirstHamiltonian+myIndex );
// 	wFWval += FirstHamiltonian;

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
//       OhmmsAttributeSet attrib;
//       attrib.add(blockT,"blockSize");
//       attrib.put(cur);
//       app_log()<<"  Forward walking block size is "<< blockT<<"*Tau"<<endl;
//       P.phLength=0;
      bool FIRST=true;

      xmlNodePtr tcur = cur->children;
      while(tcur != NULL) {
        string cname((const char*)tcur->name);
//         app_log()<<cname<<endl;
        if(cname == "Observable") 
        {
          string tagName("none");
          int Hindex(-100);
	  int blockSeries(0);
	  int blockFreq(0);
          OhmmsAttributeSet Tattrib;
          Tattrib.add(tagName,"name");
	  Tattrib.add(blockSeries,"max");
	  Tattrib.add(blockFreq,"frequency");
          Tattrib.put(tcur);
	  
	  int numProps = P.PropertyList.Names.size();
// 	  Hindex = P.PropertyList.add(tagName);
	  Hindex = h.getObservable(tagName)+NUMPROPERTIES;
	  if(tagName=="LocalPotential") {
	    Hindex=LOCALPOTENTIAL ;
	    tagName="LocPot";
	  }
	  if (Hindex==(NUMPROPERTIES-1)){
	    app_log()<<"Not a valid H element("<<Hindex<<") Valid names are:";
	    for (int jk=0;jk<h.sizeOfObservables();jk++) app_log()<<"  "<<h.getObservableName(jk);
	    app_log()<<endl;
	    exit(-1);
	  }
/*
          if ((Hindex==-100)){
            app_log()<<" Hamiltonian Element "<<tagName<<" does not exist!! "<<Hindex<<endl;
            assert(Hindex>=0);
          }*/
          Names.push_back(tagName);
          Hindices.push_back( Hindex);
          app_log()<<" Hamiltonian Element "<<tagName<<" was found at "<< Hindex<<endl;
          
          vector<int> Parameters;
          if(blockSeries==0) putContent(Parameters,tcur);
	  else{
	    for( int pl=blockFreq;pl<=blockSeries;pl+=blockFreq) Parameters.push_back(pl);
	  }
          
	  int numT=Parameters.size();
	  nObservables+=1;
          nValues+=numT;
          
          app_log()<<"   "<<numT<<" values will be calculated at block numbers:"<<endl;
          app_log()<<"      ";
          for(int nm=0;nm<Parameters.size();nm++) app_log()<<Parameters[nm]<<"  ";
          app_log()<<endl;
          walkerLengths.push_back(Parameters);
          
          int maxWsize=Parameters.back();
	  int pindx = P.addPropertyHistory(maxWsize);
          Pindices.push_back(pindx);
	  app_log()<<"pindex "<<pindx<<endl;
          
        }
        tcur = tcur->next;
      }
      app_log()<<"Total number of observables calculated:"<<nObservables<<endl;
      app_log()<<"Total number of values calculated:"<<nValues<<endl;
      for (int i=0;i<nValues;i++) Values.push_back(0.0);
      EValues=Values;
      return true;
    }

    bool get(std::ostream& os) const {
      os << "TrialVCorrection";
      return true;
    }

    QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi)
    {
      TrialDMCCorrection* myClone = new TrialDMCCorrection();
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
      myClone->EValues=EValues;
      
      return myClone;
    }

    void addObservables(PropertySetType& plist)
    {
      bool first(true);
      for(int i=0;i<nObservables;i++){
        for(int j=0;j<walkerLengths[i].size();j++){
          std::stringstream sstr,sstr2;
          sstr << "T_" << Names[i] << "_" << walkerLengths[i][j];
          app_log()<<"Observables named "<<sstr.str()<<endl;
	  if (first) {
            myIndex = plist.add(sstr.str());
            first=false;
          } else {
            plist.add(sstr.str());
          }
        }
      }
      for(int i=0;i<nObservables;i++){
        for(int j=0;j<walkerLengths[i].size();j++){
          std::stringstream sstr2;
          sstr2 << "ET_" << Names[i] << "_" << walkerLengths[i][j];
          app_log()<<"Observables named "<<sstr2.str()<<endl;
          plist.add(sstr2.str());
        }
      }
    }

    void setObservables(PropertySetType& plist)
    {
      for (int i=0;i<nValues;i++) plist[myIndex+ i]=Values[i];
      for (int i=0;i<nValues;i++) plist[myIndex+i+nValues]=EValues[i];
 
    }
    
    void setParticlePropertyList(PropertySetType& plist, int offset)
    {
 
      for (int i=0;i<nValues;i++) plist[myIndex+i+offset]=Values[i];
      for (int i=0;i<nValues;i++) plist[myIndex+i+offset+nValues]=EValues[i];
 
    }
  };
}
#endif

/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1581 $   $Date: 2007-01-04 10:02:14 -0600 (Thu, 04 Jan 2007) $
 * $Id: Pressure.h 1581 2007-01-04 16:02:14Z jnkim $ 
 ***************************************************************************/


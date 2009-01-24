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
#include <cstdlib>
#include <set>
#include <string>


namespace qmcplusplus {


  struct ForwardWalking: public QMCHamiltonianBase {
    vector<int> Hindices;
    vector<int> Pindices;
    vector<vector<int> > walkerLengths;
    vector<double> Values;
    vector<string> Names;
    int blockT,nObservables,nValues,FirstHamiltonian;

    double count;


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
    


    inline Return_t 
    rejectedMove(ParticleSet& P) {
        vector<double>::iterator Vit=Values.begin();
        for(int i=0;i<nObservables;i++){
	  for(int j=0;j<walkerLengths[i].size();j++,Vit++){
            (*Vit) = tWalker->PropertyHistory[Pindices[i]][walkerLengths[i][j]-1];
          }
        }
	std::copy(Values.begin(),Values.end(),tWalker->getPropertyBase()+FirstHamiltonian+myIndex);
      return 0.0;
    }
    
    inline Return_t 
    evaluate(ParticleSet& P) {
	for(int i=0;i<nObservables;i++){
          tWalker->addPropertyHistoryPoint(Pindices[i],  P.PropertyList[Hindices[i]]);
        }

        vector<double>::iterator Vit=Values.begin();
        for(int i=0;i<nObservables;i++){
	  for(int j=0;j<walkerLengths[i].size();j++,Vit++){
            (*Vit) = tWalker->PropertyHistory[Pindices[i]][walkerLengths[i][j]-1];
          }
        }

	std::copy(Values.begin(),Values.end(),tWalker->getPropertyBase()+FirstHamiltonian+myIndex);
        return 0.0;
    }

    inline Return_t 
    evaluate(ParticleSet& P, vector<NonLocalData>& Txy) {
      return evaluate(P);
    }

    bool put(xmlNodePtr cur) {return true;}
      
    bool put(xmlNodePtr cur, QMCHamiltonian& h, ParticleSet& P) {
      FirstHamiltonian = h.startIndex();
      nObservables=0;
      nValues=0;
      blockT=1;
//       OhmmsAttributeSet attrib;
//       attrib.add(blockT,"blockSize");
//       attrib.put(cur);
      bool FIRST=true;

      xmlNodePtr tcur = cur->children;
      while(tcur != NULL) {
        string cname((const char*)tcur->name);
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
	  
	  if (tagName.find("*")==string::npos)
	  {
	  //Single Observable case
	  int numProps = P.PropertyList.Names.size();
	  Hindex = h.getObservable(tagName)+NUMPROPERTIES;
	  if(tagName=="LocalPotential") {
	    Hindex=LOCALPOTENTIAL ;
	    tagName="LocPot";
	  }
	  else if(tagName=="LocalEnergy") {
	    Hindex=LOCALENERGY ;
	    tagName="LocEn";
	  }	  
	  else if (Hindex==(NUMPROPERTIES-1)){
	    app_log()<<"Not a valid H element("<<Hindex<<") Valid names are:";
	    for (int jk=0;jk<h.sizeOfObservables();jk++) app_log()<<"  "<<h.getObservableName(jk);
	    app_log()<<endl;
	    exit(-1);
	  }

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
	  
	  }
	  else
	  {
	    bool FOUNDH(false);
	    // 	    Multiple observables for this tag
	    int found=tagName.rfind("*");
	    tagName.replace (found,1,"");
	    
	    int numProps = P.PropertyList.Names.size();
	    for(int j=0;j<h.sizeOfObservables();j++)
	    {
	      string Hname = h.getObservableName(j);
	      if (Hname.find(tagName) != string::npos)
	      {
		vector<int> Parameters;
		if(blockSeries==0) putContent(Parameters,tcur);
		else{
		  for( int pl=blockFreq;pl<=blockSeries;pl+=blockFreq) Parameters.push_back(pl);
		}
		FOUNDH=true;
		app_log()<<" Hamiltonian Element "<<Hname<<" was found at "<< j<<endl;
		Names.push_back(Hname);
		Hindex = j+NUMPROPERTIES;
		Hindices.push_back( Hindex);
		
		walkerLengths.push_back(Parameters);
		int maxWsize=Parameters.size();
		int pindx = P.addPropertyHistory(maxWsize);
		Pindices.push_back(pindx);
		
		nValues+=Parameters.size();
		app_log()<<"   "<<numT<<" values will be calculated at block numbers:"<<endl;
		app_log()<<"      ";
		for(int nm=0;nm<Parameters.size();nm++) app_log()<<Parameters[nm]<<"  ";
		app_log()<<endl;
	      }
	    }
	    if (FOUNDH)
	    {
	      nObservables+=1;
	    }
	    else
	    {
	    app_log()<<"Not a valid H element("<<Hindex<<") Valid names are:";
	    for (int jk=0;jk<h.sizeOfObservables();jk++) app_log()<<"  "<<h.getObservableName(jk);
	    app_log()<<endl;
	    exit(-1);
	    }


	  }
        }
        tcur = tcur->next;
      }
      app_log()<<"Total number of observables calculated:"<<nObservables<<endl;
      app_log()<<"Total number of values calculated:"<<nValues<<endl;
      Values.resize(nValues);
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
    
    void setParticlePropertyList(PropertySetType& plist, int offset)
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


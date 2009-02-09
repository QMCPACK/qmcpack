//////////////////////////////////////////////////////////////////
// (c) Copyright 2008-  by Jeremy McMinis and Jeongnim Kim
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
//   Department of Physics, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#include "QMCHamiltonians/trialDMCcorrection.h"
#include "QMCHamiltonians/QMCHamiltonian.h"
#include "ParticleBase/ParticleAttribOps.h"
#include "OhmmsData/ParameterSet.h"
#include "OhmmsData/AttributeSet.h"
#include <cstdlib>
#include <set>
#include <string>
namespace qmcplusplus {
      
  bool TrialDMCCorrection::putSpecial(xmlNodePtr cur, QMCHamiltonian& h, ParticleSet& P ) 
  {
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
        } else if (Hindex==(NUMPROPERTIES-1)){
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

  QMCHamiltonianBase* TrialDMCCorrection::makeClone(ParticleSet& qp, TrialWaveFunction& psi)
  {
    return new TrialDMCCorrection(*this);
  }

  void TrialDMCCorrection::addObservables(PropertySetType& plist)
  {
  }

  void TrialDMCCorrection::addObservables(PropertySetType& plist, BufferType& collectables)
  {
    myIndex=plist.size();
    int nc=0;
    for(int i=0;i<nObservables;++i)
      for(int j=0;j<walkerLengths[i].size();++j,++nc)
      {
        std::stringstream sstr;
        sstr << "T_" << Names[i] << "_" << walkerLengths[i][j];
        int id= plist.add(sstr.str());
        myIndex=std::min(myIndex,id);
        //app_log()<<"Observables named "<<sstr.str()<< " at " << id <<endl;
      }

    for(int i=0;i<nObservables;++i)
      for(int j=0;j<walkerLengths[i].size();++j,++nc)
      {
        std::stringstream sstr;
        sstr << "ET_" << Names[i] << "_" << walkerLengths[i][j];
        int id=plist.add(sstr.str());
        myIndex=std::min(myIndex,id);
        //app_log()<<"Observables named "<<sstr.str()<< " at " << id <<endl;
      }
    app_log()<<"TrialDMCCorrection::Observables [" << myIndex << ", " << myIndex+nc << ")" << endl;
  }
}

/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1581 $   $Date: 2007-01-04 10:02:14 -0600 (Thu, 04 Jan 2007) $
 * $Id: trialDMCcorrection.h 1581 2007-01-04 16:02:14Z jnkim $ 
 ***************************************************************************/


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
#ifndef QMCPLUSPLUS_TRIALENERGY_H
#define QMCPLUSPLUS_TRIALENERGY__H
#include "Particle/ParticleSet.h"
#include "Particle/WalkerSetRef.h"
#include "QMCHamiltonians/QMCHamiltonianBase.h"




namespace qmcplusplus {


  struct TrialEnergy: public QMCHamiltonianBase {
    vector<int> walkerLengths;
    vector<double> Values;
    int nValues,historyIndex,numT ;

    /** constructor
     *
     * Pressure operators need to be re-evaluated during optimization.
     */
    TrialEnergy(ParticleSet& P) { 
      UpdateMode.set(OPTIMIZABLE,1);
    }
    ///destructor
    ~TrialEnergy() { }

    void resetTargetParticleSet(ParticleSet& P) {
      }

    inline Return_t 
    evaluate(ParticleSet& P) {
      Value = tWalker->Properties(TRIALENERGY);
      tWalker->addPropertyHistoryPoint(historyIndex,Value);

      vector<double>::iterator Vit=Values.begin();
        for(int i=0;i<numT;i++,Vit++){
            (*Vit) = tWalker->PropertyHistory[historyIndex][walkerLengths[i]-1];
        }
      return 0.0;
    }

    inline Return_t 
    evaluate(ParticleSet& P, vector<NonLocalData>& Txy) {
      return evaluate(P);
    }
    
    /** implements the virtual function.
     * 
     * Nothing is done but should check the mass
     */
    
    bool put(xmlNodePtr cur,ParticleSet& P ) {
      putContent(walkerLengths, cur);
      numT=walkerLengths.size();
      app_log()<<"   "<<numT<<" E_T values will be calculated at step numbers:"<<endl;
      app_log()<<"      ";
      for(int nm=0;nm<walkerLengths.size();nm++) app_log()<<walkerLengths[nm]<<"  ";
      app_log()<<endl;
      nValues=walkerLengths.back();
      historyIndex = P.addPropertyHistory(nValues);
      Values.resize(numT);
      return true;
    }
    bool put(xmlNodePtr cur) {
      return true;
    }

    bool get(std::ostream& os) const {
      os << "TrialEnergy";
      return true;
    }

    QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi)
    {
      TrialEnergy* myClone = new TrialEnergy(qp);
      ///Do I need to do this?
      myClone->numT=numT;
      myClone->walkerLengths=walkerLengths;
      myClone->nValues=nValues;
      myClone->myIndex=myIndex;
      myClone->historyIndex=historyIndex;
      myClone->Values=Values;
      return myClone;
    }
    
    void addObservables(PropertySetType& plist)
    {
      bool first(true);
        for(int j=0;j<walkerLengths.size();j++){
          std::stringstream sstr;
          sstr << "E_T_" << walkerLengths[j];
          app_log()<<"Observables named "<<sstr.str()<<endl;
          if (first) {
            myIndex = plist.add(sstr.str());
            first=false;
          } else {
            plist.add(sstr.str());
          }
      }
    }

    void setObservables(PropertySetType& plist)
    {
      for (int i=0;i<numT;i++) plist[myIndex+i]=Values[i];
    }

  };
}
#endif

/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1581 $   $Date: 2007-01-04 10:02:14 -0600 (Thu, 04 Jan 2007) $
 * $Id: Pressure.h 1581 2007-01-04 16:02:14Z jnkim $ 
 ***************************************************************************/


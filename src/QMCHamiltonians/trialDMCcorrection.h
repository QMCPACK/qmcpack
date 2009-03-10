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
//   Materials Computation Center, UIUC
//   Department of Physics, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef QMCPLUSPLUS_TRIALDMCCORRECTION_H
#define QMCPLUSPLUS_TRIALDMCCORRECTION_H
#include "QMCHamiltonians/QMCHamiltonianBase.h"

namespace qmcplusplus {

  class QMCHamiltonian;

  struct TrialDMCCorrection: public QMCHamiltonianBase 
  {
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
      //UpdateMode.set(OPTIMIZABLE,1);
    }
    
    ///destructor
    ~TrialDMCCorrection() { }

    void resetTargetParticleSet(ParticleSet& P) { }

    inline Return_t rejectedMove(ParticleSet& P) 
    {
      vector<double>::iterator Vit=Values.begin();
      vector<double>::iterator Vit2=EValues.begin();

      for(int i=0;i<nObservables;i++){
        //           app_log()<<"Obs#"<<i;
        (*Vit)=0.0;
        int k=0;
        int DMindex = tWalker->PHindex[Pindices[i]]-1;
        for(int j=0;j<tWalker->PropertyHistory[Pindices[i]].size();j++){
          if (DMindex<0) DMindex=tWalker->PropertyHistory[Pindices[i]].size()-1;
          (*Vit) += tWalker->PropertyHistory[Pindices[i]][DMindex];
          DMindex--;
          if(j==walkerLengths[i][k]){
            double Tsum=(*Vit);
            (*Vit2)=Tsum* (tWalker->Properties(LOCALENERGY));
            Vit++;
            Vit2++;
            if (Vit != Values.end())
            {
              (*Vit)=Tsum;
            }
            k++;
          }
          //             app_log()<<"  "<<tWalker->PropertyHistory[Pindices[i]][walkerLengths[i][j]-1];
        }
        // 	  app_log()<<endl;
      }
      //       }
      // 	double* wFWval = tWalker->getPropertyBase();
      std::copy(Values.begin(),Values.end(),tWalker->getPropertyBase()+FirstHamiltonian+myIndex );
      std::copy(EValues.begin(),EValues.end(),tWalker->getPropertyBase()+FirstHamiltonian+myIndex+nValues);

      return 0.0;
  }
    
  inline Return_t evaluate(ParticleSet& P) 
  {

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
        int DMindex = tWalker->PHindex[Pindices[i]]-1;
        for(int j=0;j<tWalker->PropertyHistory[Pindices[i]].size();j++){
          if (DMindex<0) DMindex=tWalker->PropertyHistory[Pindices[i]].size()-1;
          (*Vit) += tWalker->PropertyHistory[Pindices[i]][DMindex];
          DMindex--;
          if(j==walkerLengths[i][k]){
            double Tsum=(*Vit);
            (*Vit2)=Tsum* (tWalker->Properties(LOCALENERGY));
            Vit++;
            if (Vit != Values.end())
            {
              (*Vit)=Tsum;
              (*Vit2)=Tsum* (tWalker->Properties(LOCALENERGY));
              Vit2++;
            }
            k++;
          }
          //             app_log()<<"  "<<tWalker->PropertyHistory[Pindices[i]][walkerLengths[i][j]-1];
        }
        //    app_log()<<endl;
      }
    //       }
    // 	double* wFWval = tWalker->getPropertyBase();
    std::copy(Values.begin(),Values.end(),tWalker->getPropertyBase()+FirstHamiltonian+myIndex );
    std::copy(EValues.begin(),EValues.end(),tWalker->getPropertyBase()+FirstHamiltonian+myIndex+nValues);
    // 	wFWval += FirstHamiltonian;

    return 0.0;
  }

  inline Return_t evaluate(ParticleSet& P, vector<NonLocalData>& Txy) 
  {
    return evaluate(P);
  }

  bool put(xmlNodePtr cur) {return true;}
      
  bool putSpecial(xmlNodePtr cur, QMCHamiltonian& h, ParticleSet& P );

  bool get(std::ostream& os) const {
    os << "TrialVCorrection";
    return true;
  }

  QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi);

  void addObservables(PropertySetType& plist); 
  void addObservables(PropertySetType& plist, BufferType& collectables);

  void setObservables(PropertySetType& plist)
  {
    std::copy(Values.begin(),Values.end(),plist.begin()+myIndex);
    std::copy(EValues.begin(),EValues.end(),plist.begin()+myIndex+nValues);
    //for (int i=0;i<nValues;i++) plist[myIndex+ i]=Values[i];
    //for (int i=0;i<nValues;i++) plist[myIndex+i+nValues]=EValues[i];
  }
    
  void setParticlePropertyList(PropertySetType& plist, int offset)
  {
    std::copy(Values.begin(),Values.end(),plist.begin()+myIndex+offset);
    std::copy(EValues.begin(),EValues.end(),plist.begin()+myIndex+nValues+offset);
    //for (int i=0;i<nValues;i++) plist[myIndex+i+offset]=Values[i];
    //for (int i=0;i<nValues;i++) plist[myIndex+i+offset+nValues]=EValues[i];
  }
};
}
#endif

/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1581 $   $Date: 2007-01-04 10:02:14 -0600 (Thu, 04 Jan 2007) $
 * $Id: trialDMCcorrection.h 1581 2007-01-04 16:02:14Z jnkim $ 
 ***************************************************************************/


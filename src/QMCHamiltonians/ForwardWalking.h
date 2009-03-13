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
/** @file ForwardWalking.h
 * @brief Declarations of ForwardWalking
 */
#ifndef QMCPLUSPLUS_FORWARDWALKING_H
#define QMCPLUSPLUS_FORWARDWALKING_H
#include "QMCHamiltonians/QMCHamiltonianBase.h"

namespace qmcplusplus {

  class QMCHamiltonian;

  struct ForwardWalking: public QMCHamiltonianBase 
  {
    vector<int> Hindices;
    vector<int> Pindices;
    vector<vector<int> > walkerLengths;
    vector<double> Values;
    vector<string> Names;
    int blockT,nObservables,nValues,FirstHamiltonian;
    double count;

    /** constructor
     */
    ForwardWalking() 
    {
      UpdateMode.set(OPTIMIZABLE,1);
    }
    
    ///destructor
    ~ForwardWalking() { }

    void resetTargetParticleSet(ParticleSet& P) { }

    inline Return_t rejectedMove(ParticleSet& P) 
    {
      vector<Return_t>::iterator Vit=Values.begin();
      
      for(int i=0;i<nObservables;i++){
        int j=0; 
        int FWindex = tWalker->PHindex[Pindices[i]]; 
        while (j<walkerLengths[i].size())
        {
          int Cindex = FWindex - walkerLengths[i][j];
          if (Cindex< 0) Cindex += tWalker->PropertyHistory[Pindices[i]].size();
          (*Vit) = tWalker->PropertyHistory[Pindices[i]][Cindex];
          j++;
          Vit++;
        }
      }
      std::copy(Values.begin(),Values.end(),tWalker->getPropertyBase()+FirstHamiltonian+myIndex);
//       cerr<<"REJECTEDMOVE"<<endl;
//       Vit=Values.begin();
//       while(Vit!=Values.end())
//       {
//         cerr<<(*Vit)<<"  ";
//         Vit++;
//       }
//       cerr<<endl;
//       cerr<<"REJECTEDMOVE ARRAY"<<endl;
//       for(int dindex=0;dindex<tWalker->PropertyHistory.size();dindex++)
//       {
//         for(int f=0;f<20;f++) cerr<<tWalker->PropertyHistory[dindex][f]<<"  ";
//         cerr<<endl; 
//       }
//       cout<<tWalker->PropertyHistory[dindex].front()<<endl; 
      return 0.0;
    }

    inline Return_t evaluate(ParticleSet& P) 
    {
      for(int i=0;i<nObservables;i++)
        tWalker->addPropertyHistoryPoint(Pindices[i],  P.PropertyList[Hindices[i]]);

      vector<double>::iterator Vit=Values.begin();
//       for(int i=0;i<nObservables;i++)
//         for(int j=0;j<walkerLengths[i].size();j++,Vit++)
//           (*Vit) = tWalker->PropertyHistory[Pindices[i]][walkerLengths[i][j] ];
      for(int i=0;i<nObservables;i++){
        int j=0;
        int FWindex = tWalker->PHindex[Pindices[i]];
        while (j<walkerLengths[i].size())
        {
          int Cindex = FWindex - walkerLengths[i][j];
          if (Cindex< 0) Cindex += tWalker->PropertyHistory[Pindices[i]].size();
          (*Vit) = tWalker->PropertyHistory[Pindices[i]][Cindex];
          j++;
          Vit++;
        }
      }
      std::copy(Values.begin(),Values.end(),tWalker->getPropertyBase()+FirstHamiltonian+myIndex);
//       cerr<<"ACCEPTED VALUES"<<endl;
//       Vit=Values.begin();
//       while(Vit!=Values.end())
//       {
//         cerr<<(*Vit)<<"  ";
//         Vit++;
//       }
//       cerr<<endl;
//       cerr<<"ACCEPTED ARRAY"<<endl;
//       for(int dindex=0;dindex<tWalker->PropertyHistory.size();dindex++)
//       {
//         for(int f=0;f<20;f++) cerr<<tWalker->PropertyHistory[dindex][f]<<"  ";
//         cerr<<endl; 
//       }
      return 0.0;
    }

    inline Return_t evaluate(ParticleSet& P, vector<NonLocalData>& Txy) 
    {
      return evaluate(P);
    }

    bool put(xmlNodePtr cur) {return true;}

    ///rename it to avoid conflicts with put
    bool putSpecial(xmlNodePtr cur, QMCHamiltonian& h, ParticleSet& P);

    bool get(std::ostream& os) const {
      os << "ForwardWalking";
      return true;
    }

    QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi);

    void addObservables(PropertySetType& plist);

    void addObservables(PropertySetType& plist, BufferType& collectables);

    void setObservables(PropertySetType& plist)
    {
      std::copy(Values.begin(),Values.end(),plist.begin()+myIndex);
    }
    
    void setParticlePropertyList(PropertySetType& plist, int offset)
    {
      std::copy(Values.begin(),Values.end(),plist.begin()+myIndex+offset);
    }
  };
}
#endif

/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1581 $   $Date: 2007-01-04 10:02:14 -0600 (Thu, 04 Jan 2007) $
 * $Id: ForwardWalking.h 1581 2007-01-04 16:02:14Z jnkim $ 
 ***************************************************************************/

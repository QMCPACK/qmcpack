//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#ifndef QMCPLUSPLUS_TRIALDMCCORRECTION_H
#define QMCPLUSPLUS_TRIALDMCCORRECTION_H
#include "QMCHamiltonians/QMCHamiltonianBase.h"

namespace qmcplusplus
{

class QMCHamiltonian;

struct TrialDMCCorrection: public QMCHamiltonianBase
{
  std::vector<int> Hindices;
  std::vector<int> Pindices;
  std::vector<std::vector<int> > walkerLengths;
  std::vector<double> Values,EValues,FWValues;
  std::vector<std::string> Names;
  int resum,CountIndex,nObservables,nValues,FirstHamiltonian;

  double count;


  /** constructor
   *
   * Pressure operators need to be re-evaluated during optimization.
   */
  TrialDMCCorrection()
  {
    //UpdateMode.set(OPTIMIZABLE,1);
  }

  ///destructor
  ~TrialDMCCorrection() { }

  void resetTargetParticleSet(ParticleSet& P) { }

  inline Return_t rejectedMove(ParticleSet& P)
  {
    for (int i=0; i<nObservables; i++)
    {
      int lastindex = tWalker->PHindex[Pindices[i]]-1;
      if (lastindex<0)
        lastindex +=walkerLengths[i][2];//tWalker->PropertyHistory[Pindices[i]].size();
      tWalker->addPropertyHistoryPoint(Pindices[i],  tWalker->PropertyHistory[Pindices[i]][lastindex]  );
    }
    calculate(P);
    return 0.0;
  }

  inline void calculate(ParticleSet& P)
  {
    std::vector<double>::iterator Vit=Values.begin();
    std::vector<double>::iterator Vit2=EValues.begin();
    std::vector<double>::iterator Vit3=FWValues.begin();
    Return_t LocEn = tWalker->Properties(LOCALENERGY);
    std::vector<std::vector<Return_t> >& restrict walkerProps = (tWalker->PropertyHistory);
    std::vector<int>& restrict walkerPHindex = tWalker->PHindex;
    if (walkerProps[CountIndex][0] >= resum)
    {
      walkerProps[CountIndex][0]=0;
      for (int i=0; i<nObservables; i++)
      {
        walkerProps[Pindices[i]+1][0] = 0.0;
        int DMindex = walkerPHindex[Pindices[i]]-2;
        int k=0;
        for(int j=0; j<walkerLengths[i][1] ; )
        {
          if (DMindex<0)
            DMindex+=walkerLengths[i][2];
          walkerProps[Pindices[i]+1][j] += walkerProps[Pindices[i]][DMindex];
          DMindex--;
          k++;
          if(k==walkerLengths[i][0])
          {
            double Tsum=walkerProps[Pindices[i]+1][j];
            (*Vit)=Tsum;
            (*Vit2)=Tsum*LocEn ;
            (*Vit3)=walkerProps[Pindices[i]][DMindex+1];
            Vit++;
            Vit2++;
            Vit3++;
            k=0;
            j++;
            if (j<walkerLengths[i][1])
            {
              walkerProps[Pindices[i]+1][j]=Tsum;
            }
          }
        }
      }
    }
    else
    {
      walkerProps[CountIndex][0] +=1;
      for (int i=0; i<nObservables; i++)
      {
        int DMindex = walkerPHindex[Pindices[i]]-2;
        int FWindex = walkerPHindex[Pindices[i]]-1;
        if (DMindex<0)
          DMindex+=walkerLengths[i][2];
        Return_t newVal = walkerProps[Pindices[i]][DMindex];
        for (int j=0; j<walkerLengths[i][1]; j++)
        {
          walkerProps[Pindices[i]+1][j] += newVal;
          FWindex -= walkerLengths[i][0];
          if (FWindex<0)
            FWindex+=walkerLengths[i][2];
          (*Vit3) = walkerProps[Pindices[i]][FWindex];
          Vit3++;
        }
        int hin = DMindex;
        for (int j=0; j<walkerLengths[i][1]; j++)
        {
          hin -= walkerLengths[i][0];
          if (hin<0)
            hin+=walkerLengths[i][2];
          walkerProps[Pindices[i]+1][j] -=  walkerProps[Pindices[i]][hin];
          (*Vit)=walkerProps[Pindices[i]+1][j];/*tWalker->PropertyHistory[Pindices[i]+1][j];*/
          (*Vit2)=LocEn*walkerProps[Pindices[i]+1][j];
          Vit++;
          Vit2++;
        }
        //         for(int j=0;j<tWalker->PropertyHistory[Pindices[i]].size();j++){
        //           if (DMindex<0) DMindex=tWalker->PropertyHistory[Pindices[i]].size()-1;
        //           (*Vit) += tWalker->PropertyHistory[Pindices[i]][DMindex];
        //           DMindex--;
        //           if(j==walkerLengths[i][k]){
        //             double Tsum=(*Vit);
        //             (*Vit2)=Tsum* (tWalker->Properties(LOCALENERGY));
        //             Vit++;
        //             Vit2++;
        //             if (Vit != Values.end())
        //             {
        //               (*Vit)=Tsum;
        //             }
        //             k++;
        //           }
        //           //             app_log()<<"  "<<tWalker->PropertyHistory[Pindices[i]][walkerLengths[i][j]-1];
        //         }
        //    app_log()<< std::endl;
      }
    }
    copy(Values.begin(),Values.end(),tWalker->getPropertyBase()+FirstHamiltonian+myIndex);
    copy(EValues.begin(),EValues.end(),tWalker->getPropertyBase()+FirstHamiltonian+myIndex+nValues);
    copy(FWValues.begin(),FWValues.end(),tWalker->getPropertyBase()+FirstHamiltonian+myIndex+nValues+nValues);
  }

  inline Return_t evaluate(ParticleSet& P)
  {
    for (int i=0; i<nObservables; i++)
    {
      tWalker->addPropertyHistoryPoint(Pindices[i],  P.PropertyList[Hindices[i]]);
    }
    calculate(P);
    return 0.0;
  }

  inline Return_t evaluate(ParticleSet& P, std::vector<NonLocalData>& Txy)
  {
    return evaluate(P);
  }

  bool put(xmlNodePtr cur)
  {
    return true;
  }

  bool putSpecial(xmlNodePtr cur, QMCHamiltonian& h, ParticleSet& P);

  bool get(std::ostream& os) const
  {
    os << "TrialVCorrection";
    return true;
  }

  QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi);

  void addObservables(PropertySetType& plist);
  void addObservables(PropertySetType& plist, BufferType& collectables);

  void setObservables(PropertySetType& plist)
  {
    copy(Values.begin(),Values.end(),plist.begin()+myIndex);
    copy(EValues.begin(),EValues.end(),plist.begin()+myIndex+nValues);
    copy(FWValues.begin(),FWValues.end(),plist.begin()+myIndex+nValues+nValues);
    //for (int i=0;i<nValues;i++) plist[myIndex+ i]=Values[i];
    //for (int i=0;i<nValues;i++) plist[myIndex+i+nValues]=EValues[i];
  }

  void setParticlePropertyList(PropertySetType& plist, int offset)
  {
    copy(Values.begin(),Values.end(),plist.begin()+myIndex+offset);
    copy(EValues.begin(),EValues.end(),plist.begin()+myIndex+nValues+offset);
    copy(FWValues.begin(),FWValues.end(),plist.begin()+myIndex+nValues+nValues+offset);
    //for (int i=0;i<nValues;i++) plist[myIndex+i+offset]=Values[i];
    //for (int i=0;i<nValues;i++) plist[myIndex+i+offset+nValues]=EValues[i];
  }
};
}
#endif



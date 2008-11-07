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
#ifndef QMCPLUSPLUS_DMCPRESSURE_CORR_H
#define QMCPLUSPLUS_DMCPRESSURE_CORR_H
#include "Particle/ParticleSet.h"
#include "Particle/WalkerSetRef.h"
#include "QMCHamiltonians/QMCHamiltonianBase.h"
#include "ParticleBase/ParticleAttribOps.h"
#include "OhmmsData/ParameterSet.h"



namespace qmcplusplus {

  /** @ingroup hamiltonian
   @brief Evaluate the Bare Pressure.
   P=/frac{2T+V}{d* /Omega}
   where d is the dimension of space and /Omega is the volume.
  **/

  struct DMCPressureCorr: public QMCHamiltonianBase {
    double pNorm, Value0, Value1, Value2;
    int phlen;
    int ELVindex, Vindex, Eindex;
    Walker<Return_t, ParticleSet::ParticleGradient_t>* tWalker;

    /** constructor
     *
     * Pressure operators need to be re-evaluated during optimization.
     */
    DMCPressureCorr(ParticleSet& P, int phLen) {
      UpdateMode.set(OPTIMIZABLE,1);
      pNorm = 1.0/(P.Lattice.DIM*P.Lattice.Volume);
//       app_log()<<"Tau "<<Tau<<endl;
      phlen=phLen;
      ELVindex = P.addPropertyHistory(phlen);
      Vindex = P.addPropertyHistory(phlen);
      int tpl= max(phlen,1000);
      Eindex = P.addPropertyHistory(tpl);
      app_log()<<"INITIAL  indices :"<<ELVindex<<"  "<<Vindex<<"  "<<Eindex<<endl;
    }
    
    DMCPressureCorr(ParticleSet& P) {
      UpdateMode.set(OPTIMIZABLE,1);
      pNorm = 1.0/(P.Lattice.DIM*P.Lattice.Volume);
    }
    
    ///destructor
    ~DMCPressureCorr() { }

    void resetTargetParticleSet(ParticleSet& P) {
      pNorm = 1.0/(P.Lattice.DIM*P.Lattice.Volume);
    }
    
    void setHistories(Walker<Return_t, ParticleSet::ParticleGradient_t>& ThisWalker){
      (tWalker)= &(ThisWalker);
//       app_log()<<"  Plist for ThisWalker is len:"<<ThisWalker.PropertyHistory.size()<<endl;
    }

    inline Return_t 
    evaluate(ParticleSet& P) {
//       app_log()<<"  Plist for P is len:"<<P.PropertyHistory.size()<<endl;
//       app_log()<<"  Plist for twalker is len:"<<tWalker->PropertyHistory.size()<<endl;
      
//       app_log()<<"  indices :"<<ELVindex<<"  "<<Vindex<<"  "<<Eindex<<endl;
      double eloc = P.PropertyList[LOCALENERGY];
      double vloc = P.PropertyList[LOCALPOTENTIAL];
      tWalker->addPropertyHistoryPoint(ELVindex,eloc*vloc);
      tWalker->addPropertyHistoryPoint(Vindex,vloc);
      tWalker->addPropertyHistoryPoint(Eindex,eloc);
      double ELVsum = tWalker->getPropertyHistorySum(ELVindex,phlen);
      double Vsum = tWalker->getPropertyHistorySum(Vindex,phlen);
      double Ebar = tWalker->getPropertyHistoryAvg(Eindex);
      
      int AC(0);
      double tm=1.0;
      for (int i=1;((i<tWalker->PropertyHistory[Eindex].size())&&(tm>0.0));i++){
        tm=tWalker->PropertyHistory[Eindex][i]-Ebar;
        AC+=1;
      }

      Value0=Vsum*pNorm*Tau;
      Value1=ELVsum*pNorm*Tau;
      Value2=AC;
      return 0.0;
    }

    inline Return_t 
    evaluate(ParticleSet& P, vector<NonLocalData>& Txy) {
      return evaluate(P);
    }

    bool put(xmlNodePtr cur ) {return true;}

    bool get(std::ostream& os) const {
      os << "DMCPressureCorr";
      return true;
    }

    QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi)
    {
      DMCPressureCorr* myClone = new DMCPressureCorr(qp);
      myClone->ELVindex =ELVindex;
      myClone->Vindex = Vindex;
      myClone->Eindex = Eindex;
      return myClone;
    }
    
    void addObservables(PropertySetType& plist)
    {
      myIndex=plist.add("SumPotDMC");
      plist.add("ElSumPotDMC");
      plist.add("decorr");

    }

    void setObservables(PropertySetType& plist)
    {
      plist[myIndex]=Value0;
      plist[myIndex+1]=Value1;
      plist[myIndex+2]=Value2;
    }
    void setParticlePropertyList(PropertySetType& plist, int offset)
    {
      plist[myIndex+offset]=Value0;
      plist[myIndex+1+offset]=Value1;
      plist[myIndex+2+offset]=Value2;
      
    }
  };
}
#endif

/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1581 $   $Date: 2007-01-04 10:02:14 -0600 (Thu, 04 Jan 2007) $
 * $Id: Pressure.h 1581 2007-01-04 16:02:14Z jnkim $ 
 ***************************************************************************/


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
    }
    
    ///destructor
    ~DMCPressureCorr() { }

    void resetTargetParticleSet(ParticleSet& P) {
      pNorm = 1.0/(P.Lattice.DIM*P.Lattice.Volume);
//       app_log()<<"Tau "<<Tau<<endl;
//       assert(phlen==P.phLength);
    }

    inline Return_t 
    evaluate(ParticleSet& P) {
//       assert(phlen==P.phLength);
      
      double eloc = P.PropertyList[LOCALENERGY];
      double vloc = P.PropertyList[LOCALPOTENTIAL];
      P.addPropertyHistoryPoint(ELVindex,eloc*vloc);
      P.addPropertyHistoryPoint(Vindex,vloc);
      P.addPropertyHistoryPoint(Eindex,eloc);
      double ELVsum = P.getPropertyHistorySum(ELVindex,phlen);
      double Vsum = P.getPropertyHistorySum(Vindex,phlen);
      double Ebar = P.getPropertyHistoryAvg(Eindex);
      
      int AC(0.0);
      double tm=1.0;
      for (int i=1;((i<phlen)&&(tm>0.0));i++){
        tm=P.PropertyHistory[Eindex][i]-Ebar;
        AC+=1;
      }
//       app_log()<<"DATAS "<<Vsum<<"  "<<ELVsum<<"  "<<AC<<endl;
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
      return new DMCPressureCorr(qp,phlen);
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

  };
}
#endif

/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1581 $   $Date: 2007-01-04 10:02:14 -0600 (Thu, 04 Jan 2007) $
 * $Id: Pressure.h 1581 2007-01-04 16:02:14Z jnkim $ 
 ***************************************************************************/


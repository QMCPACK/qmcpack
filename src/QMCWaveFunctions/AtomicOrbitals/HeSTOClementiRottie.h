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
#ifndef OHMMS_QMC_HEPRESETHF_H
#define OHMMS_QMC_HEPRESETHF_H
#include "QMCWaveFunctions/OrbitalBuilderBase.h"

namespace ohmmsqmc {

  /**class HePresetHF
   *@brief a specialized class that implments He-HF orbital
   *Used to debug the package.
   */
  struct HePresetHF: public QMCTraits {
    
    enum {N=5};
    TinyVector<RealType,N> C;
    TinyVector<RealType,N> Z, ZZ;
    HePresetHF(){
      C[0]=0.76838;C[1]=0.22346;C[2]=0.04082;C[3]=-0.00994;C[4]=0.00230;
      Z[0]=1.41714;Z[1]=2.37682;Z[2]=4.39628;Z[3]=6.52699;Z[4]=7.94252;
      const RealType fourpi = 4.0*(4.0*atan(1.0));
      for(int i=0; i<N; i++) {
	//C[i] *= sqrt(pow(2.0*Z[i],3.0)/(2.0*fourpi));
	C[i] *= sqrt(pow(2.0*Z[i],3.0)/2)/fourpi;
      }
      for(int i=0; i<N; i++) ZZ[i] = Z[i]*Z[i];
    }
    
    inline void reset() { }
    
    inline void resizeByWalkers(int nw) { }
    
    template<class VM, class GM>
    inline void 
    evaluate(const ParticleSet& P, int first, int last,
	     VM& logdet, GM& dlogdet, VM& d2logdet) {
      RealType r = myTable->r(first);
      
      RealType rinv = myTable->rinv(first);
      PosType dr = myTable->dr(first);
      RealType rinv2 = rinv*rinv;
      RealType chi = 0.0, d2chi = 0.0;
      PosType dchi;
      for(int i=0; i<C.size(); i++) {
	RealType u = C[i]*exp(-Z[i]*r);
	RealType du = -u*Z[i]*rinv; // 1/r du/dr 
	RealType d2u = u*ZZ[i]; // d2u/dr2
	chi += u;
	dchi += du*dr;
	d2chi += (d2u+2.0*du); 
      }
      logdet(0,0) =chi;
      dlogdet(0,0) = dchi;
      d2logdet(0,0) = d2chi; 
    }

    template<class VM, class GM>
    inline void 
    evaluate(const WalkerSetRef& W, int first, int last,
	     vector<VM>& logdet, vector<GM>& dlogdet, vector<VM>& d2logdet) {

      int nptcl = last-first;
      for(int iw=0; iw<W.walkers(); iw++) {
	int nn = first;///first pair of the particle subset
	RealType r = myTable->r(iw,first);
	RealType rinv = myTable->rinv(iw,first);
	PosType dr = myTable->dr(iw,first);
	RealType rinv2 = rinv*rinv;
	RealType chi = 0.0, d2chi = 0.0;
	PosType dchi;
	for(int i=0; i<C.size(); i++) {
	  RealType u = C[i]*exp(-Z[i]*r);
	  RealType du = -u*Z[i]*rinv; // 1/r du/dr 
	  RealType d2u = u*ZZ[i]; // d2u/dr2
	  chi += u;
	  dchi += du*dr;
	  d2chi += (d2u+2.0*du); 
	}
	logdet[iw](0,0) =chi;
	dlogdet[iw](0,0) = dchi;
	d2logdet[iw](0,0) = d2chi; 
      }
    }

    void setTable(DistanceTableData* dtable) { myTable = dtable;}
    DistanceTableData* myTable;
  };


  struct HePresetHFBuilder: public OrbitalBuilderBase {

    HePresetHFBuilder(TrialWaveFunction& wfs,
		      ParticleSet& ions,
		      ParticleSet& els);
    bool put(xmlNodePtr cur);

  };

}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

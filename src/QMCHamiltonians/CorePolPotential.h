//////////////////////////////////////////////////////////////////
// (c) Copyright 2003  by Jeongnim Kim and Jordan Vincent
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
#ifndef OHMMS_QMC_COREPOLPOTENTIAL_H
#define OHMMS_QMC_COREPOLPOTENTIAL_H
#include <algo.h>
#include "Particle/ParticleSet.h"
#include "Particle/WalkerSetRef.h"
#include "Particle/DistanceTableData.h"
#include "QMCHamiltonians/QMCHamiltonianBase.h"

namespace ohmmsqmc {

  struct CPPelel: public QMCHamiltonianBase {
    int Centers;
    RealType alpha, r_b, r_binv;
    RealType C;
    DistanceTableData* d_table;
    vector<bool> CoreCoef;
    Matrix<PosType> CoreCoreDipole, ElCoreDipole;

    //index for attribute charge
    int iz = ions.Species.addAttribute("charge");

    CPPelel(ParticleSet& ions, ParticleSet& els): 
      d_table(NULL) { 
      
      d_table = DistanceTable::getTable(DistanceTable::add(ions,els));
      d_ii;
      Centers = ions.getTotalNum();
      alpha = 0.3558;
      C = -0.5*alpha;
      r_b = 0.7048;
      r_binv = 1.0/r_b;

      CoreCoef.resize(Centers);
      CoreCoreDipole.resize(Centers,Centers);
      ElCoreDipole.reisze(Centers,els.getTotalNum());

      for(int i=0; i<Centers; i++){
	string sname = ions.Species.speciesName[ions.GroupID[i]];
	if(sname == "Ge"){
	  cout << "Adding a core-electron potential for " << sname << endl;
	  CoreCoef[i] = true;
	}
	else CoreCoef[i] = false;
      }

      int nn=0;
      for(int i=0; i<Centers; i++) {
	for(int j=i+1;j<Centers; j++, nn++) {
	  RealType rinv3 = pow(d_ii->rinv(nn),3);
	  PosType dipole = rinv3*d_ii->dr(nn);
	  CoreCoreDipole(i,j) = dipole*ions.Species(iz,ions.GroupID[j]);//charge of j
	  CoreCoreDipole(j,i) = dipole*ions.Species(iz,ions.GroupID[i]);//charge of i
	}
      }
    }
    
    ~CPPelel() { }

    inline ValueType evaluate(ParticleSet& P) {

      RealType esum=0.0;
      int nn = P.getTotalNum();
      vector<RealType> fc(nn),rinv2(nn);
      vector<PosType> dr(nn);
      int nnoff = 0;
      for(int iat=0; iat<Centers; iat++) {
	if(CoreCoef[iat]) {
	  for(int j=0; j<nn; j++, nnoff++) {
	    dr[j] = d_table->dr(nnoff);
	    rinv2[j] = pow(d_table->rinv(nnoff),2); // 1/r^2
	    fc[j] = fcpp(d_table->r(nnoff)*r_binv); //*rinv2[nnj]*d_table->rinv(j); // f(r)/r^3
	    enum += pow(rinv2[j],2)*fc[j];
	  }

	  for(int j=0; j<nn; j++) {
            for(int k=j+1; k<nn; k++) {
              esum += dot(dr[j],dr[k])*rinv2[j]*rinv2[k]*fc[j]*fc[k];
	    }
	  }
	  for(int j=nnoff, nnj=0; nnj<nn-1; nnj++,j++) {
	    esum += CoreCoef[iat]*fc[nnj]*fc[nnj]*rinv2[nnj];
	    for(int k=j+1,nnk=nnj+1; nnk<nn; nnk++,k++) {
	      esum += fc[nnj]*fc[nnk]*dot(d_table->dr(j),d_table->dr(k));
	    }
	  }
	  for(int jat=iat+1; jat<Centers; jat++) {
	    int i=nnoff;
	    int j=d_table->M[jat];
	    for(int k=0; nnj=0; nni=0; k<nn; j++; i++; nni++; nnj++; k++){
	      esum += dot(CoreCoreDipole(iat,jat),d_table->dr(j))*fc[nnj];
	      esum += dot(CoreCoreDipole(jat,iat),d_table->dr(i))*fc[nni];
	  }
	}
	}
      } //iat


      return C*esum;
    }


    inline ValueType 
    evaluate(ParticleSet& P, RealType& x) {
      RealType esum=0.0;
      int nn = P.getTotalNum();
      vector<RealType> fc(nn),rinv2(nn);
      for(int iat=0; iat<Centers; iat++) {
	int nnoff = d_table->M[iat];
	for(int j=nnoff, nnj=0; nnj<nn; j++, nnj++) {
	  rinv2[nnj] = pow(d_table->rinv(j),2);
	  fc[nnj] = fcpp(d_table->r(j)*r_binv)*rinv2[nnj]*d_table->rinv(j);
	}
	for(int j=nnoff, nnj=0; nnj<nn-1; nnj++,j++) {
	  esum += CoreCoef[iat]*fc[nnj]*fc[nnj]*rinv2[nnj];
	  for(int k=j+1,nnk=nnj+1; nnk<nn; nnk++,k++) {
            esum += CoreCoef[iat]*fc[nnj]*fc[nnk]*dot(d_table->dr(j),d_table->dr(k));
	  }
	}
      }

      x = C*esum;
      return x;
    }

    inline RealType fcpp(RealType z) {
      return pow((1.0-exp(-1.0*z*z)),2.0);
    }


#ifdef WALKERSTORAGE_COLUMNMAJOR
    inline void 
    evaluate(WalkerSetRef& P, ValueVectorType& LE) {
      ValueVectorType e(P.walkers());
      for(int iw=0; iw<P.walkers(); iw++) {
	e=0.0;
	for(int iat=0; iat<Centers; iat++) {
	  for(int nn=d_table->M[iat]; nn<d_table->M[iat+1]; nn++) {
	    e[iw] += d_table->rinv(iw,nn);
	  }
	}
      }
      for(int iw=0; iw<P.walkers(); iw++) { LE[iw] += e[iw];} 
    }
#else
    inline void 
    evaluate(WalkerSetRef& P, ValueVectorType& LE) {
      for(int iw=0; iw<P.walkers(); iw++) {
	RealType e=0.0;
	for(int iat=0; iat<Centers; iat++) {
	  RealType esub = 0.0;
	  for(int nn=d_table->M[iat]; nn<d_table->M[iat+1]; nn++) {
	    esub += d_table->rinv(iw,nn);
	  }
	  e += esub; ///multiply z
	}
	LE[iw] += e;
      }
    }
#endif
  };
}
#endif


//  struct CPPcorecore: public QMCHamiltonianBase {

//    RealType alpha;
//    RealType C;
//    DistanceTableData* d_table;
    
//    CPPcorecore(ParticleSet& ref): d_sum(0.0) { 

//      alpha = 0.3558; 
//      C = -0.5*alpha;

//      IndexType iii = DistanceTable::add(ref);
//      DistanceTableData* d_ii = DistanceTable::getTable(iii);
//      d_ii->create(1);
     
//      int charge = ref.Species.addAttribute("charge");
//      int nat = ref.getTotalNum();
//      vector<RealType> Z(nat);
//      for(int iat=0; iat<nat;iat++) {
//        Z[iat] = ref.Species(charge,ref.GroupID[iat]);
//      }
//      d_ii->evaluate(ref);
//      d_sum = 0.0;
//      for(int i=0; i< nat; i++) {
//        PosType fcc = 0.0;
//        for(int nn=d_ii->M[i]; nn<d_ii->M[i+1]; nn++) {
//           fcc += Z[ref.GroupID[d_ii->J[nn]]]*pow(d_ii->rinv(nn),3.0)*d_ii->dr(nn);
// 	}
//         d_sum += C*dot(fcc,fcc);
//       }
//       LOGMSG("Energy of ion-ion interaction " << d_sum)
//    }
   
//    ~CPPcorecore() { }

//    inline ValueType evaluate(ParticleSet& P) { return d_sum; }

//    inline ValueType evaluate(ParticleSet& P, RealType& x) {
//      ///only for consistency
//      return x=d_sum;
//    }
   
//  };

//  struct CPPelcore: public QMCHamiltonianBase {
//    RealType alpha;
//    int Centers;
//    vector<RealType> Z;
//    DistanceTableData* d_table;
//    DistanceTableData* d_ii;
    
//    CPPelcore(ParticleSet& ions, ParticleSet& els): d_table(NULL) { 
     
//      alpha = 0.3558; 
//      IndexType iii = DistanceTable::add(ions);
//      DistanceTableData* d_ii = DistanceTable::getTable(iii);
//      d_ii->create(1);

//      d_table = DistanceTable::getTable(DistanceTable::add(ions,els));
//      int iz = ions.Species.addAttribute("charge");
//      Centers = ions.getTotalNum();
//      Z.resize(Centers);
     
//      int charge = ref.Species.addAttribute("charge");
//      int nat = ref.getTotalNum();
//      vector<RealType> Z(nat);
//      vector<PosType> fcc(nat);
//      for(int iat=0; iat<nat;iat++) {
//        Z[iat] = ref.Species(charge,ref.GroupID[iat]);
//      }
//      d_ii->evaluate(ref);
//      d_sum = 0.0;
//      for(int i=0; i<nat; i++) {
//        for(int nn=d_ii->M[i]; nn<d_ii->M[i+1]; nn++) {
//           fcc[i] += Z[ref.GroupID[d_ii->J[nn]]]*pow(d_ii->rinv(nn),3.0)*d_ii->dr(nn);
// 	}
//       }
//    }
   
//    ~CPPelcore() { }

//    inline ValueType evaluate(ParticleSet& P) {

//      d_ie->evaluate(P);
//      for(int iat=0; iat<Centers; iat++) {
//        int eln = 0;
// 	for(int nn=d_ie->M[i]; nn<d_ie->M[iat+1]; nn++; eln++) {
// 	  esum+=alpha*fcc[iat]*pow(d_ie->rinv[nn],3.0)*dot(ions.R[iat],P.R[eln]);
// 	}
//      }
//      return esum;
//    }

//  inline ValueType evaluate(ParticleSet& P, RealType& x) {

//      d_ie->evaluate(P);
//      for(int iat=0; iat<Centers; iat++) {
//        int eln = 0;
// 	for(int nn=d_ie->M[i]; nn<d_ie->M[iat+1]; nn++; eln++) {
// 	  esum+=alpha*fcc[iat]*pow(d_ie->rinv[nn],3.0)*dot(ions.R[iat],P.R[eln]);
// 	}
//      }
//      x = esum;
//      return esum;
//    }

//  };



/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/


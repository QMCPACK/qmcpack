//////////////////////////////////////////////////////////////////
// (c) Copyright 2003  by Jeongnim Kim and Dyutiman Das
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
#ifndef OHMMS_QMC_WOSPOTENTIAL_H
#define OHMMS_QMC_WOSPOTENTIAL_H
#include <algo.h>
#include <vector>
#include "Particle/ParticleBase.h"
#include "Particle/WalkerSetRef.h"
#include "Particle/DistanceTableData.h"
#include "QMCHamiltonians/QMCHamiltonianBase.h"
#include "QMCHamiltonians/WOS/Domain.h"
#include "QMCHamiltonians/WOS/HeteroStructure.h"

namespace ohmmsqmc {

  struct WOSPotential: public QMCHamiltonianBase {
    
    typedef double RealType;
    
    DistanceTableData* d_table;
    HeteroStructure* QDot;
    
    
    /// Constructor
    WOSPotential(ParticleBase& els,
		 double delta,
		 Grid3D* aGrid3D,
		 const char* rhofile,
		 const char* vfile,
		 std::vector<double>& dielectric,
		 std::vector<double>& BandOffset,
		 std::vector<double>& allIntervals) : 
      d_table(NULL), QDot(NULL) { 

      d_table = DistanceTable::getTable(DistanceTable::add(els));

      QDot = new HeteroStructure(delta,
				 aGrid3D,
				 rhofile,
				 vfile,
				 dielectric,
				 BandOffset,
				 allIntervals);

    }
    
    ~WOSPotential() { }

    
    inline ValueType evaluate(ParticleBase& P) {

      RealType PE = 0.0; RealType PEsq = 0.0;

      RealType vHartree = 0.0;

      /// create a domain
      Domain domain;

      int itermax = 1;
      //      Random.init(0,1,1356823);

      for(int iat = 0; iat < P.getTotalNum(); iat++){

	for(int iter = 0; iter < itermax ; iter++){

	  RealType vbare = 0.0;
	  RealType Gee = 0.0;
	  
	/// start runner from particle position
	  domain.runner = P.R[iat];
	  domain.in_device = true;
	  //cout << "runner: " << domain.runner << endl;

	/// dielectric at runner position
	  QDot->epsilon(domain);
	  //cout << "eps " << domain.eps_d << endl;

	/// initialise weight_bc
	  QDot->weight_init();
	  //cout << "wt: " << QDot->weight_bc << endl;

	/// create spherical domain within layer
	  QDot->MakeLimitedSphere( domain );
	  //cout << "radius: " << domain.radius << endl;

	  vbare = QDot->Image_Contribution( iat, domain ,P );
	  //	  vbare += QDot->Domain_Contribution( domain, P );
	  //cout << "vbare: " << vbare << endl;

	  QDot->sample_point( domain );
	  //cout << "r: " << domain.runner << endl;

	  QDot->MakeMaximumSphere( domain );
	  //cout << "radius " << domain.radius << endl;
	
	  vbare += QDot->Domain_Contribution( domain, P );
	  //cout << "vbare: " << vbare << endl;

	  vbare += QDot->passage( domain );
	  //cout << "vbare: " << vbare << '\t' << domain.in_device << endl;

	  while( domain.in_device ) {

	    QDot->epsilon( domain );
	    //cout << "eps" << domain.eps_d << endl;

	    QDot->sample_point( domain );
	    //cout << "r: " << domain.runner << endl;

	    QDot->MakeMaximumSphere( domain );
	    //cout << "d " << domain.radius << endl;

	    vbare += QDot->Domain_Contribution( domain, P );
	    //cout << "vbare: " << vbare << endl;

	    vbare += QDot->passage( domain );
	    //cout << "vbare: " << vbare << '\t' << domain.in_device << endl;

	  }

	  PE += vbare;
	  PEsq += vbare*vbare;
	  //	  cout << iter << '\t' << vbare << endl;
	  //cout << "poten:: " << vbare << '\t' << QDot->weight_bc << endl;
	}
      
	PE /= 1.0*itermax; PEsq /= 1.0*itermax;

	//	cout <<  P.R[0][0] << '\t' << -PE << '\t' << sqrt(fabs(PEsq-PE*PE)/(1.0*itermax)) << endl;
      //      //cout << " I am here:: very good " << P.R[0] << endl;
	//cout << endl << endl << endl << endl ;

      }

      return -PE*0.036749309;
    } 
    

    /*      RealType esum = 0.0;
	    int nn = 0;
	    for(int iat=0; iat<P.getTotalNum(); iat++) {
	    for(int jat=iat+1; jat<P.getTotalNum(); jat++) {
	    esum += d_table->r(nn++);
	    }*/
    //for(int nn=d_table->M[iat]; nn<d_table->M[iat+1]; nn++)
    //esum+=d_table->r(nn); 
  
    

  

    inline ValueType
    evaluate(ParticleBase& P, RealType& x){
      return x = evaluate(P);
    } 
  
#ifdef USE_FASTWALKER
    inline void 
    evaluate(WalkerSetRef& W, ValueVectorType& LE) {
    }
#else
    inline void 
    evaluate(WalkerSetRef& W, ValueVectorType& LE) {
    }
#endif
  };
}
#endif

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/


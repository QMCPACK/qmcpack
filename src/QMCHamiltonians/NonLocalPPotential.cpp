//////////////////////////////////////////////////////////////////
// (c) Copyright 2004 by Jeongnim Kim and Jordan Vincent
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
//#include "Particle/MCWalkerConfiguration.h"
#include "QMCHamiltonians/NonLocalPPotential.h"
#include "Particle/DistanceTableData.h"
#include "Particle/DistanceTable.h"
#include "Numerics/OhmmsBlas.h"
#include "Numerics/HDFNumericAttrib.h"
#include "Utilities/SimpleParser.h"
#include "Utilities/OhmmsInfo.h"
#include "Utilities/RandomGenerator.h"

namespace ohmmsqmc {

  NonLocalPPotential::RadialPotentialSet::~RadialPotentialSet() {
    for(int ig=0; ig<lgrid_m.size(); ig++)  delete lgrid_m[ig];
    for(int ip=0; ip<lpp_m.size(); ip++)    delete lpp_m[ip];
    for(int ig=0; ig<nlgrid_m.size(); ig++) delete nlgrid_m[ig];
    for(int ip=0; ip<nlpp_m.size(); ip++)   delete nlpp_m[ip];
  }

  NonLocalPPotential::ValueType 
  NonLocalPPotential::RadialPotentialSet::evaluate(
      ParticleSet& W, DistanceTableData* d_table, int iat, 
      TrialWaveFunction& Psi){

    RealType esum=0.0;

    int nchannel=nlpp_m.size();
    int nknot=sgridxyz_m.size();

    randomize_grid();

    //int iel=0;
    for(int nn=d_table->M[iat],iel=0; nn<d_table->M[iat+1]; nn++,iel++){

      if(d_table->r(nn)>Rmax) continue;

      register RealType r(d_table->r(nn));
      register RealType rinv(d_table->rinv(nn));
      register PosType  dr(d_table->dr(nn));
      // Compute ratio of wave functions
      for (int j=0; j < nknot ; j++){ 
	PosType deltar(r*rrotsgrid_m[j]-dr);
	W.makeMove(iel,deltar);
	psiratio[j]=Psi.ratio(W,iel)*sgridweight_m[j];
      }
      // Compute radial potential
      for(int ip=0;ip< nchannel; ip++){
	vrad[ip]=nlpp_m[ip]->evaluate(r,rinv)*wgt_angpp_m[ip];
      }
      // Compute spherical harmonics on grid
      for (int j=0, jl=0; j<nknot ; j++){ 
	RealType zz=dot(dr,rrotsgrid_m[j])*rinv;
	lpol[0]=1.0;
	RealType lpolprev=0.0;
	// Forming the Legendre polynomials
	for (int l=0 ; l< lmax ; l++){
	  lpol[l+1]=(2*l+1)*zz*lpol[l]-l*lpolprev;
	  lpol[l+1]/=(l+1);
	  lpolprev=lpol[l];
	}
	//for (int l=0; l < nchannel; l++)
	//  Amat[l*nknot+j] = lpol[ angpp_m[l] ];
        for(int l=0; l <nchannel; l++,jl++) Amat[jl]=lpol[ angpp_m[l] ];
      } 

      BLAS::gemv(nknot, nchannel, &Amat[0], &psiratio[0], &wvec[0]);
      esum += BLAS::dot(nchannel, &vrad[0], &wvec[0]);
      ////////////////////////////////////
      //Original implmentation by S. C.
      //const char TRANS('T');
      //const int ione=1;
      //const RealType one=1.0;
      //const RealType zero=0.0;
      //dgemv(TRANS,nknot,nchannel,one,&Amat[0],nknot,&psiratio[0],ione,zero,&wvec[0],ione);
      //esum += ddot(nchannel,&vrad[0],ione,&wvec[0],ione);
      ////////////////////////////////////
      //iel++;
    }   /* end loop over electron */
    return esum;
  }

  /** constructor
   *\param ions the positions of the ions
   *\param els the positions of the electrons
   *\param psi trial wavefunction
   *
   * For each ion-type, an ASCII file "*.psf" must
   be provided.  The "*.psf" must contain two columns, 
   the first column being the grid, the second being
   the potential on the grid. 
   */

  NonLocalPPotential::NonLocalPPotential(ParticleSet& ions, ParticleSet& els,
      TrialWaveFunction& psi):
    Centers(ions.GroupID), d_table(NULL), Psi(psi)
  { 

    d_table = DistanceTable::getTable(DistanceTable::add(ions,els));

    const SpeciesSet& Species(ions.getSpeciesSet());
    for(int i=0; i<Species.getTotalNum();i++) {
      vector<RealType> grid_temp, pp_temp;
      string species(Species.speciesName[i]);
      string fname = species+".psf";
      ifstream fin(fname.c_str(),ios_base::in);
      if(!fin){
	ERRORMSG("Could not open file " << fname)
	  exit(-1);
      }      
      XMLReport("Reading a file for the PseudoPotential for species " << species);
      // Add a new center to the list.
      PP.push_back(new RadialPotentialSet);
      // Read Number of potentials (local and non) for this atom
      int npotentials;
      fin >> npotentials;
      RealType r, f1;
      int numnonloc=0, lmax=-1;
      RealType rmax(0.0);
      for (int ij=0; ij<npotentials; ij++){
	int angmom,npoints;
	fin >> angmom >> npoints;
        if(grid_temp.size()<npoints) grid_temp.resize(npoints);
        if(pp_temp.size()<npoints) pp_temp.resize(npoints);
	for (int j=0; j<npoints; j++){
          fin >> grid_temp[j] >> pp_temp[j];
	  //fin >> r ; grid_temp.push_back(r);
	  //fin >> f1; pp_temp.push_back(f1);
	}
	GridType *agrid = new NumericalGrid<ValueType>(grid_temp);
	LocalPotentialType *app = new OneDimCubicSpline<ValueType>(agrid,pp_temp);
	int imin = 0;
	RealType yprime_i = ((*app)(imin+1)-(*app)(imin))/app->dr(imin);
        LOGMSG("NonPP l=" << angmom << " deriv= " << yprime_i)
	app->spline(imin,yprime_i,app->size()-1,0.0);
	if(angmom < 0)
	  PP[i]->add(agrid,app);
	else {
	  PP[i]->add(angmom,agrid,app);
	  numnonloc++;
	}
	lmax=std::max(lmax,angmom);
	rmax=std::max(rmax,agrid->rmax());
      }
      //cout << npotentials << " potentials read" << endl;
      PP[i]->lmax=lmax; PP[i]->Rmax=rmax;
      fin.close();
      int numsgridpts=0;
      // if NonLocal look for file containing the spherical grid
      if(numnonloc>0){
	string fname = species+".sgr";
	ifstream fin(fname.c_str(),ios_base::in);
	if(!fin){
	  ERRORMSG("Could not open file " << fname)
	  exit(-1);
	}
	PosType xyz;
	ValueType weight;
	while(fin >> xyz >> weight){
	  PP[i]->addknot(xyz,weight);
	  numsgridpts++;
	}
	//cout << "Spherical grid : " << numsgridpts << " points" <<endl;
      }

      PP[i]->resize_warrays(numsgridpts,numnonloc,lmax);
    }//centers
  }

  ///destructor
  NonLocalPPotential::~NonLocalPPotential() { 
    for(int pp=0; pp<PP.size(); pp++) delete PP[pp];
  }

  ///random rotation of the spherical grid
  void NonLocalPPotential::RadialPotentialSet::randomize_grid(){
    const RealType twopi(6.28318530718);
    RealType phi(twopi*Random()),psi(twopi*Random()),cth(Random()-0.5),
	     sph(sin(phi)),cph(cos(phi)),sth(sqrt(1-cth*cth)),sps(sin(psi)),
	     cps(cos(psi));
    Tensor<double,3> rmat( cph*cth*cps-sph*sps, sph*cth*cps+cph*sps,-sth*cps,
                          -cph*cth*sps-sph*cps,-sph*cth*sps+cph*cps, sth*sps,
			   cph*sth,             sph*sth,             cth     );
    vector<PosType>::iterator it(sgridxyz_m.begin());
    vector<PosType>::iterator it_end(sgridxyz_m.end());
    vector<PosType>::iterator jt(rrotsgrid_m.begin());
    while(it != it_end) {*jt = dot(rmat,*it); ++it; ++jt;}
  }
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/


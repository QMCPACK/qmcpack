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
#include "Particle/MCWalkerConfiguration.h"
#include "Particle/DistanceTableData.h"
#include "Particle/DistanceTable.h"
#include "QMCHamiltonians/QMCHamiltonianBase.h"
#include "QMCHamiltonians/NonLocalPPotential.h"
#include "Numerics/OneDimGridBase.h"
#include "Numerics/OneDimGridFunctor.h"
#include "Numerics/OneDimCubicSpline.h"
#include "Numerics/Blasf.h"
#include "Numerics/HDFNumericAttrib.h"
#include "OhmmsPETE/Tensor.h"
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

  /*!
   *\param ions the positions of the ions
   *\param els the positions of the electrons
   *\brief the constructor
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

    maxsgridpts=0, maxnonloc=0; maxangmom=0;
    for(int i=0; i<ions.Species.getTotalNum();i++) {
      vector<ValueType> grid_temp, pp_temp;
      string species = ions.Species.speciesName[i];
      string fname = species+".psf";
      ifstream fin(fname.c_str(),ios_base::in);
      if(!fin){
	ERRORMSG("Could not open file " << fname)
	  exit(-1);
      }      
      XMLReport("Reading a file for the PseudoPotential for species " << species);
      // Add a new center to the list.
      PP.push_back(new RadialPotentialSet);
      // Read Number of channels for this atom
      int npotentials;
      fin >> npotentials;
      ValueType r, f1;
      int angmom,npoints;
      int numnonloc=0;
      RealType rmax(0.0);
      for (int ij=0; ij<npotentials; ij++){
	fin >> npoints >> angmom;
	for (int j=0; j<npoints; j++){
	  fin >> r ; grid_temp.push_back(r);
	  fin >> f1; pp_temp.push_back(f1);
	}
	GridType *agrid = new NumericalGrid<ValueType>(grid_temp);
	LocalPotentialType *app = new OneDimCubicSpline<ValueType>(agrid,pp_temp);
	int imin = 1;
	RealType yprime_i = ((*app)(imin+1)-(*app)(imin))/app->dr(imin);
	app->spline(imin,yprime_i,app->size()-1,0.0);
	if(angmom < 0)
	  PP[i]->add(agrid,app);
	else {
	  PP[i]->add(angmom,agrid,app);
	  numnonloc++;
	}
	rmax=std::max(rmax,agrid->rmax());
	//PP[i]->Rmax=agrid->rmax();
      }
      PP[i]->Rmax=rmax;
      if(numnonloc>maxnonloc)maxnonloc=numnonloc;
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
	  numnonloc++;
	}
	if(numsgridpts>maxsgridpts) maxsgridpts=numsgridpts;
      }
      else 
	PP[i]->lmax=-1;
      if(PP[i]->lmax > maxangmom) maxangmom=PP[i]->lmax;
    }//centers
    // If there are NonLocal components => allocate working arrays
    if(maxnonloc>0){
      psiratio = new ValueType[maxsgridpts];
      vrad = new ValueType[maxnonloc];
      wvec = new ValueType[maxnonloc];
      Amat = new ValueType[maxnonloc*maxsgridpts];
      lpol = new ValueType[maxangmom];
    }
  }

  ///destructor
  NonLocalPPotential::~NonLocalPPotential() { 
    for(int pp=0; pp<PP.size(); pp++) delete PP[pp];
    if(maxnonloc>0){
      delete [] psiratio ;
      delete [] vrad ;
      delete [] wvec ;
      delete [] Amat ;
      delete [] lpol ;
    }
  }

  NonLocalPPotential::ValueType 
  NonLocalPPotential::evaluate(ParticleSet& W) {

    MCWalkerConfiguration& Wref=dynamic_cast<MCWalkerConfiguration&>(W);

    const ValueType oneoverfourpi=0.0795774715;
    PosType knot,deltar;
    RealType esum=0.0;
    int iel;
    //loop over all the ions
    for(int iat=0; iat<Centers.size(); iat++) {
      int ispecies = Centers[iat];
      // evaluate the local parts
      RadialPotentialSet *pp=PP[ispecies];

      esum += pp->evaluate(d_table,iat);
      int nchannel = pp->nlpp_m.size();
      if(nchannel==0) continue;
      int nknot = pp->sgridxyz_m.size();
      iel=0;
      pp->randomize_grid();

      register RealType rcut(pp->Rmax);

      for(int nn=d_table->M[iat]; nn<d_table->M[iat+1]; nn++){

	if(d_table->r(nn)>rcut) continue;

	register RealType r(d_table->r(nn));
	register RealType rinv(d_table->rinv(nn));
	register PosType dr(d_table->dr(nn));

	// Compute ratio of wave functions
	for (int j=0; j < nknot ; j++){ 
	  knot = pp->getknot(j);
	  deltar = dr - r * knot;
	  Wref.makeMove(iel,deltar);
	  psiratio[j]=Psi.ratio(Wref,iel) * pp->sgridweight_m[j];
	}
	// Compute radial potential
	for(int ig=0; ig< nchannel; ig++)
	  pp->nlgrid_m[ig]->index(r);
	for(int ip=0;ip< nchannel; ip++){
	  vrad[ip]=pp->nlpp_m[ip]->evaluate(r,rinv);
	  vrad[ip]*=(2*pp->angpp_m[ip]+1);
	}
	// Compute spherical harmonics on grid
	for (int j=0; j<nknot ; j++){ 
	  RealType zz= dot(dr, pp->getknot(j));
	  lpol[0]=1.0;
	  lpol[1]=zz;
	  // Forming the Legendre polynomial
	  for (int l=1 ; l< pp->lmax ; l++){
	    lpol[l+1]=(2*l+1)*zz*lpol[l]-l*lpol[l-1];
	    lpol[l+1]/=(l+1);
	  }
	  for (int l=0; l < nchannel; l++)
	    Amat[l*nknot+j] = lpol[ pp->angpp_m[l] ];
	} 
	const char TRANS('T');
	const int ione=1;
	const RealType one=1.0;
	const RealType zero=0.0;
        dgemv(TRANS,nknot,nchannel,one,Amat,nknot,psiratio,ione,zero,
	    wvec,ione);
	esum += ddot(nchannel,vrad,ione,wvec,ione)*oneoverfourpi;
	iel++;
      }   /* end loop over electron */
    }
    return esum;
  }


  void NonLocalPPotential::RadialPotentialSet::randomize_grid(){
    const double twopi=6.28318531;
    double phi,psi,sph,cph,sth,cth,sps,cps;
    Tensor<double,3> rmat;
    /* Rotation matrix with Euler angles defined as: 
       1) counterclockwise rotation around z (phi)
       2) clockwise rotation around y (theta)
       3) counterclockwise rotation around z (psi) */
    phi=twopi*Random(); sph=sin(phi); cph=cos(phi);
    psi=twopi*Random(); sps=sin(psi); cps=cos(psi);
    cth=2*(Random()-0.5); sth=sqrt(1-cth*cth);
    rmat(0,0)= cph*cth*cps-sph*sps;
    rmat(0,1)= sph*cth*cps+cph*sps;
    rmat(0,2)=-sth*cps;
    rmat(1,0)=-cph*cth*sps-sph*cps;
    rmat(1,1)=-sph*cth*sps+cph*cps;
    rmat(1,2)= sth*sps;
    rmat(2,0)= cph*sth;
    rmat(2,1)= sph*sth;
    rmat(2,2)= cth;
    for (vector<PosType>::iterator it=sgridxyz_m.begin(); 
	it != sgridxyz_m.end(); ++it) *it=dot(rmat,*it);
  }
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/


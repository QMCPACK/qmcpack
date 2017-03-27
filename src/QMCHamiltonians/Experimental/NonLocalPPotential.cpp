//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
//#include "Particle/MCWalkerConfiguration.h"
#include "QMCHamiltonians/NonLocalPPotential.h"
//#include "Numerics/HDFNumericAttrib.h"
#include "Particle/DistanceTableData.h"
#include "Particle/DistanceTable.h"
#include "Utilities/SimpleParser.h"
#include "Utilities/RandomGenerator.h"

namespace qmcplusplus
{

NonLocalPPotential::RadialPotentialSet::~RadialPotentialSet()
{
  if(lpp_m)
    delete lpp_m;
  if(lgrid_m)
    delete lgrid_m;
  //This has to be fixed
  //for(int ig=0; ig<nlgrid_m.size(); ig++) delete nlgrid_m[ig];
  for(int ip=0; ip<nlpp_m.size(); ip++)
    delete nlpp_m[ip];
}

void NonLocalPPotential::RadialPotentialSet::resize_warrays(int n,int m,int l)
{
  psiratio.resize(n);
  vrad.resize(m);
  wvec.resize(m);
  Amat.resize(n*m);
  lpol.resize(l+1,1.0);
  rrotsgrid_m.resize(n);
  HasNonLocalPP = nlpp_m.size()>0;
  if(HasNonLocalPP)
  {
    nchannel=nlpp_m.size();
    nknot=sgridxyz_m.size();
    //This is just to check
    //for(int nl=1; nl<nlpp_m.size(); nl++) nlpp_m[nl]->setGridManager(false);
    Lfactor1.resize(lmax);
    Lfactor2.resize(lmax);
    for(int nl=0; nl<lmax; nl++)
    {
      Lfactor1[nl]=static_cast<RealType>(2*nl+1);
      Lfactor2[nl]=1.0e0/static_cast<RealType>(nl+1);
    }
    LOGMSG(" Found Non-Local Potential. Number of channels=" << 1)
  }
}

NonLocalPPotential::ValueType
NonLocalPPotential::RadialPotentialSet::evaluate(ParticleSet& W,
    DistanceTableData* d_table, int iat, TrialWaveFunction& Psi, bool randomize)
{
  RealType esum=0.0;
  if(HasNonLocalPP)
    randomize_grid(*(W.Sphere[iat]),randomize);
  //int iel=0;
  for(int nn=d_table->M[iat],iel=0; nn<d_table->M[iat+1]; nn++,iel++)
  {
    register RealType r(d_table->r(nn));
    esum += lpp_m->splint(r);
    if(r>Rmax)
      continue;
    register RealType rinv(d_table->rinv(nn));
    register PosType  dr(d_table->dr(nn));
    // Compute ratio of wave functions
    for (int j=0; j < nknot ; j++)
    {
      PosType deltar(r*rrotsgrid_m[j]-dr);
      W.makeMove(iel,deltar);
      psiratio[j]=Psi.ratio(W,iel)*sgridweight_m[j];
      W.rejectMove(iel);
    }
    // Compute radial potential
    for(int ip=0; ip< nchannel; ip++)
    {
      vrad[ip]=nlpp_m[ip]->splint(r)*wgt_angpp_m[ip];
    }
    // Compute spherical harmonics on grid
    for (int j=0, jl=0; j<nknot ; j++)
    {
      RealType zz=dot(dr,rrotsgrid_m[j])*rinv;
      // Forming the Legendre polynomials
      lpol[0]=1.0;
      RealType lpolprev=0.0;
      for (int l=0 ; l< lmax ; l++)
      {
        //Not a big difference
        //lpol[l+1]=(2*l+1)*zz*lpol[l]-l*lpolprev;
        //lpol[l+1]/=(l+1);
        lpol[l+1]=Lfactor1[l]*zz*lpol[l]-l*lpolprev;
        lpol[l+1]*=Lfactor2[l];
        lpolprev=lpol[l];
      }
      for(int l=0; l <nchannel; l++,jl++)
        Amat[jl]=lpol[ angpp_m[l] ];
    }
    if(nchannel==1)
    {
      esum += vrad[0]*BLAS::dot(nknot, &Amat[0],&psiratio[0]);
    }
    else
    {
      BLAS::gemv(nknot, nchannel, &Amat[0], &psiratio[0], &wvec[0]);
      esum += BLAS::dot(nchannel, &vrad[0], &wvec[0]);
    }
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

void NonLocalPPotential::resetTargetParticleSet(ParticleSet& P)
{
  d_table = DistanceTable::add(IonConfig,P);
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
 If there is any non-local Pseudo potential term, UpdateMode is
 set so that the values are re-evaluated during the optimizations.
 */
NonLocalPPotential::NonLocalPPotential(ParticleSet& ions, ParticleSet& els,
                                       TrialWaveFunction& psi):
  IonConfig(ions), Centers(ions.GroupID), d_table(0), Psi(&psi)
{
  d_table = DistanceTable::add(ions,els);
  els.resizeSphere(ions.getTotalNum());
  const SpeciesSet& Species(ions.getSpeciesSet());
  for(int i=0; i<Species.getTotalNum(); i++)
  {
    std::vector<RealType> grid_temp, pp_temp;
    std::string species(Species.speciesName[i]);
    std::string fname = species+".psf";
    std::ifstream fin(fname.c_str(),std::ios_base::in);
    if(!fin)
    {
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
    GridType* _nlgrid=0;
    GridType* _lgrid=0;
    for (int ij=0; ij<npotentials; ij++)
    {
      int angmom,npoints;
      fin >> angmom >> npoints;
      if(grid_temp.size()<npoints)
        grid_temp.resize(npoints);
      if(pp_temp.size()<npoints)
        pp_temp.resize(npoints);
      for (int j=0; j<npoints; j++)
      {
        fin >> grid_temp[j] >> pp_temp[j];
        //fin >> r ; grid_temp.push_back(r);
        //fin >> f1; pp_temp.push_back(f1);
      }
      GridType *agrid = new NumericalGrid<ValueType>(grid_temp);
      LocalPotentialType *app = new OneDimCubicSpline<ValueType>(agrid,pp_temp);
      int imin = 0;
      RealType yprime_i = ((*app)(imin+1)-(*app)(imin))/app->dr(imin);
      if(angmom < 0)
      {
        LOGMSG("NonPP l=" << angmom << " deriv= " << yprime_i)
        app->spline(imin,yprime_i,app->size()-1,0.0);
        PP[i]->add(agrid,app);
      }
      else
      {
        LOGMSG("NonPP l=" << angmom << " deriv= " << yprime_i)
        app->spline(imin,yprime_i,app->size()-1,0.0);
        //app->spline(imin,0.0,app->size()-1,0.0);
        PP[i]->add(angmom,agrid,app);
        numnonloc++;
        lmax=std::max(lmax,angmom);
        rmax=std::max(rmax,agrid->rmax());
      }
    }
    //cout << npotentials << " potentials read" << std::endl;
    PP[i]->lmax=lmax;
    PP[i]->Rmax=rmax;
    LOGMSG("   Maximum cutoff of NonPP " << rmax)
    fin.close();
    int numsgridpts=0;
    // if NonLocal look for file containing the spherical grid
    if(numnonloc>0)
    {
      std::string fname = species+".sgr";
      std::ifstream fin(fname.c_str(),std::ios_base::in);
      if(!fin)
      {
        ERRORMSG("Could not open file " << fname)
        exit(-1);
      }
      PosType xyz;
      ValueType weight;
      while(fin >> xyz >> weight)
      {
        PP[i]->addknot(xyz,weight);
        numsgridpts++;
      }
      //cout << "Spherical grid : " << numsgridpts << " points" << std::endl;
    }
    PP[i]->resize_warrays(numsgridpts,numnonloc,lmax);
    ////Has non-local pseudo potentials.
    //if(numnonloc) {
    //  UpdateMode.set(OPTIMIZABLE,1);
    //  UpdateMode.set(RATIOUPDATE,1);
    //}
  }//centers
  for(int ic=0; ic<ions.getTotalNum(); ic++)
  {
    int nsg=PP[ions.GroupID[ic]]->nknot;
    if(nsg)
      els.Sphere[ic]->resize(nsg);
  }
}

///destructor
NonLocalPPotential::~NonLocalPPotential()
{
  for(int pp=0; pp<PP.size(); pp++)
    delete PP[pp];
}

NonLocalPPotential::Return_t
NonLocalPPotential::evaluate(ParticleSet& P)
{
  Value=0.0;
  //loop over all the ions
  for(int iat=0; iat<Centers.size(); iat++)
  {
    Value += PP[Centers[iat]]->evaluate(P,d_table,iat,*Psi,UpdateMode[PRIMARY]);
  }
  return Value;
}

QMCHamiltonianBase* NonLocalPPotential::clone()
{
  ERRORMSG("Incomplete implementation of NonLocalPPotential::clone() ")
  return this;
}
}

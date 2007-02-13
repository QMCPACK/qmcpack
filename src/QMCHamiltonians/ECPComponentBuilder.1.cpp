//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim
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
#include "QMCHamiltonians/ECPComponentBuilder.h"
#include "Numerics/GaussianTimesRN.h"
#include "Numerics/Transform2GridFunctor.h"
#include "Utilities/SimpleParser.h"

namespace qmcplusplus {

  void ECPComponentBuilder::addSemiLocal(xmlNodePtr cur) {

    GridType* grid_semilocal=0;
    RealType rmax= pp_nonloc->Rmax;
    cur=cur->children;
    while(cur != NULL) {
      string cname((const char*)cur->name);
      if(cname == "grid") {
        grid_semilocal=createGrid(cur);
        rmax=grid_semilocal->rmax();
      } else if(cname == "vps") {
        //should be able to overwrite rmax
        int l=angMon[(const char*)xmlGetProp(cur,(const xmlChar*)"l")];
        Lmax = std::max(l,Lmax);
        xmlNodePtr cur1=cur->children;
        while(cur1 != NULL) {
          if(xmlStrEqual(cur1->name,(const xmlChar*)"basisGroup")) {
            pp_nonloc->add(l,createVr(cur1,grid_semilocal));
          }
          cur1=cur1->next;
        }
        NumNonLocal++;
      }
      cur=cur->next;
    }

    pp_nonloc->lmax=Lmax;
    pp_nonloc->Rmax=rmax;

    //Debug
    //ofstream fout("C.semilocal.dat");
    //fout.setf(std::ios::scientific, std::ios::floatfield);
    //fout.precision(12);
    //int ig=0;
    //while(ig<grid_semilocal->size()) {
    //  double r=(*grid_semilocal)[ig++];
    //  fout << setw(20) << r;
    //  for(int i=0; i<pp_nonloc->nlpp_m.size(); i++) {
    //    fout << setw(20) << pp_nonloc->nlpp_m[i]->operator()(ig);
    //  }
    //  fout << endl;
    //}
  }

  ECPComponentBuilder::RadialPotentialType*
  ECPComponentBuilder::createVr(xmlNodePtr cur, GridType* agrid) {
    //todo rcut should be reset if necessary
    typedef GaussianTimesRN<RealType> InFuncType;
    InFuncType a;
    a.putBasisGroup(cur);

    int ng=agrid->size();
    if(agrid->GridTag != LINEAR_1DGRID)
    {
      RealType ri=0.0;
      RealType rf=agrid->rmax();
      agrid = new LinearGrid<RealType>;
      agrid->set(ri,rf,ng);
      app_log() << "  Reset the grid for SemiLocal component to a LinearGrid. " << endl;
    }

    std::vector<RealType> v(ng);
    for(int ig=0; ig<ng; ig++) {
      v[ig]=a.f((*agrid)[ig]);
    }
    RadialPotentialType *app=new RadialPotentialType(agrid,v);
    app->spline();

    //ofstream fout("C.nonlocal.dat");
    //fout.setf(std::ios::scientific, std::ios::floatfield);
    //fout.precision(12);
    //RealType r=0;
    //while(r<agrid->rmax())
    //{
    //  RealType nvr=app->f(r);
    //  RealType avr= a.f(r);
    //  fout << setw(20) << r << setw(20) << nvr<< " " << avr << " " << nvr-avr << endl;
    //  r+=0.01;
    //}
    return app;
  }

  /** build a Local Pseudopotential
   *
   * For numerical stability, the radial function of the local pseudopotential is stored 
   * as \f$rV(r)/Z_{eff}/Z_{e}\f$ for the distance r and \f$Z_e = -1\f$. 
   * The coulomb factor $Z_{e}Z_{eff}/r$ is applied by LocalECPotential 
   * (see LocalECPotential::evaluate).
   */
  void ECPComponentBuilder::buildLocal(xmlNodePtr cur) {

    if(pp_loc) return; //something is wrong

    string vFormat("V");
    const xmlChar* vptr=xmlGetProp(cur,(const xmlChar*)"format");
    if(vptr != NULL)
    {
      vFormat=(const char*)vptr;
    }

    int vPowerCorrection=1;
    if(vFormat == "r*V")
    {
      app_log() << "  Local pseudopotential format = r*V" << endl;
      vPowerCorrection=0;
    } 
    else 
    {
      app_log() << "  Local pseudopotential format = V" << endl;
    }

    typedef GaussianTimesRN<RealType> InFuncType;
    GridType* grid_local=0;
    InFuncType vr;

    bool bareCoulomb=true;
    cur=cur->children;
    while(cur != NULL) {
      string cname((const char*)cur->name);
      if(cname == "grid") {
        grid_local=createGrid(cur,true); 
      } else if(cname == "basisGroup") {
        vr.putBasisGroup(cur,vPowerCorrection);
        bareCoulomb=false;
      }
      cur=cur->next;
    }

    if(grid_local->GridTag == CUSTOM_1DGRID)
    {
      app_error() << "  Custom grid is used. Need to recast to the linear grid" << endl;
      OHMMS::Controller->abort();
    } 
    else
    {
      vector<RealType> v;
      if(bareCoulomb) {
        app_log() << "   Bare Coulomb potential is used." << endl;
        grid_local->set(0.0,1.,3);
        v.resize(3);
        for(int ig=0; ig<3; ig++) v[ig]=1.0;
        pp_loc=new RadialPotentialType(grid_local,v);
        pp_loc->spline(0,0.0,2,0.0);
      } else {
        app_log() << "   Guassian basisGroup is used: base power " << vr.basePower << endl;
        const RealType eps=1e-12;//numeric_limits<RealType>::epsilon();//1e-12;
        RealType zinv=1.0/Zeff;
        RealType r=5.;
        bool ignore=true;
        int last=grid_local->size()-1;
        while(ignore&&last)
        {
          r=(*grid_local)[last];
          ignore=(abs(zinv*vr.f(r))<eps);
          --last;
        }
        if(last ==0)
        {
          app_error() << "  Illegal Local Pseudopotential " <<endl;
        }
        //Add the reset values here
        int ng=static_cast<int>(r/1e-3)+1;
        app_log() << "     Use a Linear Grid: [0,"<< r << "] Number of points = " << ng << endl;
        grid_local->set(0.0,r,ng);
        v.resize(ng);
        v[0]=1.0;
        for(int ig=1; ig<ng-1; ig++) {
          double r=(*grid_local)[ig];
          v[ig]=1.0-zinv*vr.f(r);
        }
        v[ng-1]=1.0;
        pp_loc=new RadialPotentialType(grid_local,v);
        pp_loc->spline(); //use the fixed conditions
      }
    }

    //ofstream fout("C.local.dat");
    //fout.setf(std::ios::scientific, std::ios::floatfield);
    //fout.precision(12);
    //RealType r=1e-6;
    //RealType zinv=1.0/Zeff;
    //while(r<5.0)
    //{
    //  RealType nvr=pp_loc->f(r);
    //  RealType avr= 1.0-zinv*r*vr.f(r);
    //  fout << setw(20) << r << setw(20) << nvr<< " " << nvr-avr << endl;
    //  r+=0.01;
    //}
  }

  ECPComponentBuilder::GridType* ECPComponentBuilder::createGrid(xmlNodePtr cur, bool useLinear)
  {
    GridType *agrid=0;
    RealType ri = 1e-5;
    RealType rf = 100.0;
    RealType ascale = -1.0e0;
    RealType astep = -1.0;
    //RealType astep = 1.25e-2;
    int npts = 1001;

    string gridType("log");
    OhmmsAttributeSet radAttrib;
    radAttrib.add(gridType,"type"); 
    radAttrib.add(npts,"npts"); 
    radAttrib.add(ri,"ri"); radAttrib.add(rf,"rf");
    radAttrib.add(ascale,"ascale"); radAttrib.add(astep,"astep");
    radAttrib.add(ascale,"scale"); radAttrib.add(astep,"step");
    radAttrib.put(cur);

    //overwrite the grid type to linear starting at 0.0
    if(useLinear)
    {
      gridType="linear";
      ri=0.0;
    }

    if(gridType == "log") {
      if(ascale>0.0) {
        agrid = new LogGridZero<RealType>;
        agrid->set(astep,ascale,npts);
      } else {
        if(ri<numeric_limits<RealType>::epsilon())
        {
          ri=numeric_limits<RealType>::epsilon();
        }
        agrid = new LogGrid<RealType>;
        agrid->set(ri,rf,npts);
      }
    } else if(gridType == "linear") {
      agrid = new LinearGrid<RealType>;
      if(astep>0.0)
      {
        npts = static_cast<int>((rf-ri)/astep)+1;
        app_log() << "   Linear grid overwrites npts = " << npts << " with step = " << astep << endl;
      }
      agrid->set(ri,rf,npts);
    }
    return agrid;
  }

  bool ECPComponentBuilder::parseCasino(const std::string& fname, xmlNodePtr cur)
  {

    RealType rmax=2.0;
    const xmlChar* rptr=xmlGetProp(cur,(const xmlChar*)"cutoff");
    if(rptr != NULL) rmax = atof((const char*)rptr);

    //first create with the input
    GridType *agrid=0;
    cur=cur->children;
    while(cur != NULL)
    {
      if(xmlStrEqual(cur->name,(const xmlChar*)"grid"))
      {
        agrid=createGrid(cur);
      }
      cur=cur->next;
    }

    const RealType d=5e-4;
    int ng=0;
    if(agrid == 0)
    {
      agrid = new LinearGrid<RealType>;
      ng=static_cast<int>(rmax/d)+1;
      agrid->set(0,rmax,ng);
    }
    else 
    {
      ng=agrid->size();
      rmax=agrid->rmax();
    }

    if(ng == 0)
    {
      app_error() << "Grid for the pseudopotentials is not defined." << endl;
      OHMMS::Controller->abort();
    }

    ifstream fin(fname.c_str(),ios_base::in);
    if(!fin)
    {
      app_error() << "Could not open file " << fname << endl;
      OHMMS::Controller->abort();
    }      

    if(pp_nonloc==0) pp_nonloc=new NonLocalECPComponent; 

    OhmmsAsciiParser aParser;
    int atomNumber;
    int npts;
    string eunits("rydberg");

    app_log() << "    ECPComponentBuilder::parseCasino" <<endl;
    aParser.skiplines(fin,1);//Header
    aParser.skiplines(fin,1);//Atomic number and pseudo-charge
    aParser.getValue(fin,atomNumber,Zeff);
    app_log() << "      Atomic number = " << atomNumber << "  Zeff = " << Zeff << endl;
    aParser.skiplines(fin,1);//Energy units (rydberg/hartree/ev):
    aParser.getValue(fin,eunits);
    app_log() << "      Unit of the potentials = " << eunits << endl;
    RealType Vprefactor = (eunits == "rydberg")?0.5:1.0;

    aParser.skiplines(fin,1);//Angular momentum of local component (0=s,1=p,2=d..)
    aParser.getValue(fin,Lmax);
    aParser.skiplines(fin,1);//NLRULE override (1) VMC/DMC (2) config gen (0 ==> input/default value)
    aParser.skiplines(fin,1);//0 0, not sure what to do yet
    aParser.skiplines(fin,1);//Number of grid points
    aParser.getValue(fin,npts);
    app_log() << "      Input Grid size = " << npts << endl;
    vector<RealType> temp(npts);

    aParser.skiplines(fin,1);//R(i) in atomic units
    aParser.getValues(fin,temp.begin(),temp.end());

    const RealType rLower=1e-4;
    //pick the lower and upper bound of the input grid
    vector<RealType>::iterator upper=upper_bound(temp.begin(),temp.end(),rmax);
    vector<RealType>::iterator lower=lower_bound(temp.begin(),temp.end(),rLower);
    int nUpper=temp.size()-1;
    int nLower=0;
    if(upper != temp.end()) {
      nUpper=upper-temp.begin();
      nLower=lower-temp.begin();
      ++upper;
    }

    //build a numerical grid of [nUpper,nLower]
    int ngKept=nUpper-nLower+1;
    //vector<RealType> grid_data(ngKept);
    Matrix<RealType> vnn(Lmax+1,npts);
    //std::copy(temp.begin()+nLower, temp.begin()+nUpper+1, grid_data.begin());
    for(int l=0; l<=Lmax; l++)
    {
      aParser.skiplines(fin,1);
      aParser.getValues(fin,vnn[l],vnn[l]+npts);
    }

    RealType vfac=-Vprefactor/Zeff;


    //temp is used to create a temporary data to spline
    vector<RealType> newP(ng),newPin(ngKept);

    for(int l=0; l<Lmax; l++)
    {
      //store corrected non-local part in Hartree
      for(int i=nLower,j=0; i<=nUpper; i++,j++) 
        newPin[j]=Vprefactor*(vnn(l,i)-vnn(Lmax,i));

      OneDimLinearSpline<RealType> inFunc;
      inFunc.assign(lower,upper,newPin.begin(),newPin.end());
      inFunc.spline();

      for(int i=1; i<ng-1; i++)
      {
        RealType r((*agrid)[i]);
        newP[i]=inFunc.splintNG(r)/r;
      }
      newP[0]=newP[1];
      newP[ng-1]=0.0;

      RadialPotentialType *app = new RadialPotentialType(agrid,newP);
      app->spline();
      pp_nonloc->add(l,app);
    }

    Lmax--;
    pp_nonloc->lmax=Lmax;
    pp_nonloc->Rmax=rmax;
    NumNonLocal=Lmax;

    {
      int locL=Lmax+1;
      newPin[0]=0.0;
      for(int i=nLower+1,j=1; j<ngKept-1; i++,j++) 
        newPin[j]= vfac*vnn(locL,i);
      newPin[ngKept-1]=1.0;

      //Time to build local
      OneDimLinearSpline<RealType> inFunc;
      inFunc.assign(lower, ++upper, newPin.begin(),newPin.end());
      inFunc.spline();

      newP[0]=0.0;
      for(int i=1; i<ng-1; i++)
      {
        newP[i]=inFunc.splintNG((*agrid)[i]);
      }
      newP[ng-1]=1.0;

      pp_loc = new RadialPotentialType(agrid,newP);
      pp_loc->spline();
    }

    SetQuadratureRule(Nrule);

    app_log() << "    Non-local pseudopotential parameters" <<endl;
    pp_nonloc->print(app_log());

    return true;
  }


} // namespace qmcPlusPlus
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1551 $   $Date: 2006-12-02 09:32:17 -0600 (Sat, 02 Dec 2006) $
 * $Id: ECPComponentBuilder.cpp 1551 2006-12-02 15:32:17Z jnkim $ 
 ***************************************************************************/

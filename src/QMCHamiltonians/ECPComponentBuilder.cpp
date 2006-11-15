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
#include "QMCHamiltonians/Ylm.h"
#include "QMCHamiltonians/FSAtomPseudoPot.h"
#include "Utilities/IteratorUtility.h"
#include "Utilities/SimpleParser.h"
#include <cmath>

namespace qmcplusplus {

  ECPComponentBuilder::ECPComponentBuilder(const string& aname):
    NumNonLocal(0), Lmax(0), Zeff(0), Species(aname), Nrule(4),
    grid_inp(0), pp_loc(0), pp_nonloc(0) {
    angMon["s"]=0; angMon["p"]=1; angMon["d"]=2; angMon["f"]=3; angMon["g"]=4;
    angMon["0"]=0; angMon["1"]=1; angMon["2"]=2; angMon["3"]=3; angMon["4"]=4;
  }

  bool ECPComponentBuilder::parse(const string& fname) {

    // build an XML tree from a the file;
    xmlDocPtr m_doc = xmlParseFile(fname.c_str());
    if (m_doc == NULL) {
      xmlFreeDoc(m_doc);
      app_error() <<"File " << fname << " is invalid" << endl;
      OHMMS::Controller->abort();
    }    
    // Check the document is of the right kind
    xmlNodePtr cur = xmlDocGetRootElement(m_doc);
    if (cur == NULL) {
      xmlFreeDoc(m_doc);
      app_error() << "Empty document" << endl;
      OHMMS::Controller->abort();
    }

    bool success=put(cur);
    xmlFreeDoc(m_doc);

    return success;
  }

  bool ECPComponentBuilder::put(xmlNodePtr cur) {
    int nk=0;
    //vector<RealType> kpts;
    vector<xmlNodePtr> semiPtr;
    cur=cur->children;
    while(cur != NULL) {
      string cname((const char*)cur->name);
      if(cname == "header") 
      {
        Zeff = atoi((const char*)xmlGetProp(cur,(const xmlChar*)"zval"));
      } 
      else if(cname == "grid")
      {
        if(grid_inp) {
          app_warning() << "   Only one global grid is accepted. Ignore this grid" << endl;
        }
        else 
        {
          grid_inp = createGrid(cur);
        }
      }
      else if(cname == "semilocal") 
      {
        semiPtr.push_back(cur);//save the pointer
      } 
      else if(cname == "local") 
      {
        buildLocal(cur);
      }
      // else if(cname == "sphericalGrid") 
      // {
      //  nk=atoi((const char*)xmlGetProp(cur,(const xmlChar*)"size"));
      //  kpts.resize(nk*4);
      //  putContent(kpts,cur);
      // } 
      cur=cur->next;
    }

    if(semiPtr.size()) {
      if(pp_nonloc==0) pp_nonloc=new NonLocalECPComponent; 
      if(pp_loc) 
      {
        for(int i=0; i<semiPtr.size(); i++) addSemiLocal(semiPtr[i]);
      }
      else 
        buildSemiLocalAndLocal(semiPtr);
    } 

    if(pp_nonloc) {
      SetQuadratureRule(Nrule);
      pp_nonloc->print(app_log());
//       if(nk == 0) { //assign default
//         RealType w=1.0/12.0;
//         pp_nonloc->addknot(PosType( 0.000000, 0.525731, 0.850651),w);
//         pp_nonloc->addknot(PosType( 0.000000,-0.525731,-0.850651),w);
//         pp_nonloc->addknot(PosType( 0.000000,-0.525731, 0.850651),w);
//         pp_nonloc->addknot(PosType( 0.000000, 0.525731,-0.850651),w);
//         pp_nonloc->addknot(PosType(-0.850651, 0.000000, 0.525731),w);
//         pp_nonloc->addknot(PosType( 0.850651, 0.000000,-0.525731),w);
//         pp_nonloc->addknot(PosType( 0.850651, 0.000000, 0.525731),w);
//         pp_nonloc->addknot(PosType(-0.850651, 0.000000,-0.525731),w);
//         pp_nonloc->addknot(PosType(-0.525731, 0.850651, 0.000000),w);
//         pp_nonloc->addknot(PosType( 0.525731,-0.850651, 0.000000),w);
//         pp_nonloc->addknot(PosType( 0.525731, 0.850651, 0.000000),w);
//         pp_nonloc->addknot(PosType(-0.525731,-0.850651, 0.000000),w);
//         pp_nonloc->resize_warrays(12,NumNonLocal,Lmax);
//       } 
//       else {
//         for(int ik=0, kk=0; ik<nk; ik++, kk+=4) 
//           pp_nonloc->addknot(PosType(kpts[kk],kpts[kk+1],kpts[kk+2]),kpts[kk+3]);
//         pp_nonloc->resize_warrays(nk,NumNonLocal,Lmax);
//       }
    }

    return true;
  }

  void ECPComponentBuilder::SetQuadratureRule(int rule) {
    using namespace std;

    int nk;
    RealType w;
    typedef enum {SINGLE, TETRA, OCTA, ICOSA} SymmType;
    SymmType symmetry;
    int lexact;
    RealType A, B, C, D;
    A = B = C = D = 0.0;
    switch (rule) {
    case 1:
      nk = 1;   symmetry = SINGLE; lexact = 0;
      A = 1.0;
      break;
    case 2:
      nk = 4;   symmetry = TETRA;  lexact = 2;
      A=0.25;
      break;
    case 3:
      nk = 6;   symmetry = OCTA;   lexact = 3;
      A=1.0/6.0;
      break;
    case 4:
      nk = 12;  symmetry = ICOSA;  lexact = 5;
      A = 1.0/12.0;
      B = 1.0/12.0;
      break;
    case 5:
      nk = 18;  symmetry = OCTA;   lexact = 5;
      A = 1.0/30.0; 
      B = 1.0/15.0;
      break;
    case 6:
      nk = 26;  symmetry = OCTA;   lexact = 7;
      A = 1.0  / 21.0;
      B = 4.0  / 105.0;
      C = 27.0 / 840.0;
      break;
    case 7:
      nk = 50;  symmetry = OCTA;   lexact = 1;
      A = 4.0/315.0;
      B = 64.0/2835.0;
      C = 27.0/1280.0;
      D = 14641.0/725760.0;
      break;
    default:
      ERRORMSG("Unrecognized spherical quadrature rule " << rule << ".");
      abort();
    }

    // First, build a_i, b_i, and c_i points
    vector<PosType> a, b, c, d;
    RealType p = 1.0/std::sqrt(2.0);
    RealType q = 1.0/std::sqrt(3.0);
    RealType r = 1.0/std::sqrt(11.0);
    RealType s = 3.0/std::sqrt(11.0);

    if (symmetry == SINGLE) {
      a.push_back (PosType(1.0, 0.0, 0.0));
    }
    else if (symmetry == TETRA) {
      a.push_back(PosType( q, q, q));
      a.push_back(PosType( q,-q,-q));
      a.push_back(PosType(-q, q,-q));
      a.push_back(PosType(-q,-q, q));
    }
    else if (symmetry == OCTA) {
      a.push_back(PosType( 1.0, 0.0, 0.0));
      a.push_back(PosType(-1.0, 0.0, 0.0));
      a.push_back(PosType( 0.0, 1.0, 0.0));
      a.push_back(PosType( 0.0,-1.0, 0.0));
      a.push_back(PosType( 0.0, 0.0, 1.0));
      a.push_back(PosType( 0.0, 0.0,-1.0));

      b.push_back(PosType(   p,   p, 0.0));
      b.push_back(PosType(   p,  -p, 0.0));
      b.push_back(PosType(  -p,   p, 0.0));
      b.push_back(PosType(  -p,  -p, 0.0));
      b.push_back(PosType(   p, 0.0,   p));
      b.push_back(PosType(   p, 0.0,  -p));
      b.push_back(PosType(  -p, 0.0,   p));
      b.push_back(PosType(  -p, 0.0,  -p));
      b.push_back(PosType( 0.0,   p,   p));
      b.push_back(PosType( 0.0,   p,  -p));
      b.push_back(PosType( 0.0,  -p,   p));
      b.push_back(PosType( 0.0,  -p,  -p));

      c.push_back(PosType(   q,   q,   q));
      c.push_back(PosType(   q,   q,  -q));
      c.push_back(PosType(   q,  -q,   q));
      c.push_back(PosType(   q,  -q,  -q));
      c.push_back(PosType(  -q,   q,   q));
      c.push_back(PosType(  -q,   q,  -q));
      c.push_back(PosType(  -q,  -q,   q));
      c.push_back(PosType(  -q,  -q,  -q));

      d.push_back(PosType(   r,   r,   s));
      d.push_back(PosType(   r,   r,  -s));
      d.push_back(PosType(   r,  -r,   s));
      d.push_back(PosType(   r,  -r,  -s));
      d.push_back(PosType(  -r,   r,   s));
      d.push_back(PosType(  -r,   r,  -s));
      d.push_back(PosType(  -r,  -r,   s));
      d.push_back(PosType(  -r,  -r,  -s));

      d.push_back(PosType(   r,   s,   r));
      d.push_back(PosType(   r,   s,  -r));
      d.push_back(PosType(   r,  -s,   r));
      d.push_back(PosType(   r,  -s,  -r));
      d.push_back(PosType(  -r,   s,   r));
      d.push_back(PosType(  -r,   s,  -r));
      d.push_back(PosType(  -r,  -s,   r));
      d.push_back(PosType(  -r,  -s,  -r));

      d.push_back(PosType(   s,   r,   r));
      d.push_back(PosType(   s,   r,  -r));
      d.push_back(PosType(   s,  -r,   r));
      d.push_back(PosType(   s,  -r,  -r));
      d.push_back(PosType(  -s,   r,   r));
      d.push_back(PosType(  -s,   r,  -r));
      d.push_back(PosType(  -s,  -r,   r));
      d.push_back(PosType(  -s,  -r,  -r));
    }
    else if (symmetry == ICOSA) {
      double t, p;  // theta and phi
      // a points
      t = 0.0;  p=0.0;
      a.push_back(PosType(std::cos(t),std::sin(t)*std::cos(p),std::sin(t)*std::sin(p)));
      t = M_PI;  p=0.0;
      a.push_back(PosType (std::cos(t),std::sin(t)*std::cos(p),std::sin(t)*std::sin(p)));
      // b points
      for (int k=0; k<5; k++) {
	t = std::atan(2.0);          p = (double)(2*k+0)*M_PI/5.0;
	b.push_back(PosType (std::cos(t), std::sin(t)*std::cos(p), std::sin(t)*std::sin(p)));
	t = M_PI-std::atan(2.0);     p = (double)(2*k+1)*M_PI/5.0;
	b.push_back(PosType (std::cos(t), std::sin(t)*std::cos(p), std::sin(t)*std::sin(p)));
      }
      // c points
      double t1 = std::acos ((2.0+std::sqrt(5.0)) / std::sqrt(15.0+6.0*std::sqrt(5.0)));
      double t2 = std::acos (      1.0            / std::sqrt(15.0+6.0*std::sqrt(5.0)));
      for (int k=0; k<5; k++) {
	t = t1; p = (double)(2*k+1)*M_PI/5.0;
	c.push_back(PosType (std::cos(t), std::sin(t)*std::cos(p), std::sin(t)*std::sin(p)));
	t = t2; p = (double)(2*k+1)*M_PI/5.0;
	c.push_back(PosType (std::cos(t), std::sin(t)*std::cos(p), std::sin(t)*std::sin(p)));
	t = M_PI - t1; p = (double)(2*k+0)*M_PI/5.0;
	c.push_back(PosType (std::cos(t), std::sin(t)*std::cos(p), std::sin(t)*std::sin(p)));
	t = M_PI - t2; p = (double)(2*k+0)*M_PI/5.0;
	c.push_back(PosType (std::cos(t), std::sin(t)*std::cos(p), std::sin(t)*std::sin(p)));
      }
    }
    // Now, construct rule
    if (std::fabs(A) > 1.0e-10) 
      for (int i=0; i<a.size(); i++)
	pp_nonloc->addknot(a[i], A);
    if (std::fabs(B) > 1.0e-10) 
      for (int i=0; i<b.size(); i++)
	pp_nonloc->addknot(b[i], B);
    if (std::fabs(C) > 1.0e-10) 
      for (int i=0; i<c.size(); i++)
	pp_nonloc->addknot(c[i], C);
    if (std::fabs(D) > 1.0e-10) 
      for (int i=0; i<d.size(); i++)
	pp_nonloc->addknot(d[i], D);

    // Allocate storage for wave function ratios
    pp_nonloc->resize_warrays(nk,NumNonLocal,Lmax);

    // Finally, check the rule for correctness
    assert (pp_nonloc->sgridxyz_m.size() == nk);
    assert (pp_nonloc->sgridweight_m.size() == nk);
    double wSum = 0.0;
    for (int k=0; k < nk; k++) {
      PosType r = pp_nonloc->sgridxyz_m[k];
      double nrm = dot(r,r);
      assert (std::fabs(nrm-1.0) < 1.0e-14);
      wSum += pp_nonloc->sgridweight_m[k];
      //cout << pp_nonloc->sgridxyz_m[k] << " " << pp_nonloc->sgridweight_m[k] << endl;
    }
    assert (std::fabs(wSum - 1.0) < 1.0e-14);
    // Check the quadrature rule
    CheckQuadratureRule(lexact);

  }

//   double ECPComponentBuilder::AssociatedLegendre(int l, int m, double x)
//   {
//     int am = abs(m);
//     // Value for positive m
//     double posVal;
//     if (l == 0)
//       return 1.0;
//     else if (l == 1) {
//       if (am == 0)
// 	return x;
//       else
// 	posVal = -std::sqrt(1.0-x*x);
//     }
//     else if (l == am) {
//       double sign = ((l%2)==1) ? -1.0 : 1.0;
//       double dfact = 1.0;
//       for (int i=(2*l-1); i>0; i-=2)
// 	dfact *= (double)i;
//       double sqt = std::sqrt(1-x*x);
//       double powval = 1.0;
//       for (int i=0; i<l; i++)
// 	powval *= sqt;
//       posVal = sign*dfact*pval;
//     }
//     else if (l == (am+1)) {
//       posVal = x*(2.0*l+1.0)*AssociatedLegendre(l, l-1, x);
//     }
//     else {
//       // Start with values for am=l and am=l-1
//       double Plm2=AssociatedLegendre(l, 

//     }
	
//     if (m < 0) {
//       double sign = ((m%2)==1) ? -1.0 : 1.0;
//       double mfact = 1.0;
//       double pfact = 1.0;
//       for (int i=(l-am); i>0; i--)
// 	mfact *= (double)i;
//       for (int i=(l+am); i>0; i--)
// 	pfact *= (double)i;
//     }

//   }

//   complex<double> ECPComponentBuilder::Ylm (int l, int m, PosType r)
//   {
//     double costheta, phi;
//     costheta = r[0];
//     phi   = std::atan2(r[2],r[1]);
//     int lmm = l - m;
//     int lpm = l + m;
//     double mfact = 1.0;
//     double pfact = 1.0;
//     for (int i=lmm; i>0; i--)
//       mfact *= (double)i;
//     for (int i=lpm; i>0; i--)
//       pfact *= (double)i;
//     double prefactor = std::sqrt ((double)(2*l+1)*mfact/(4.0*M_PI*pfact));


//     double Pl;
//     if (l == 0)
//       Pl = 1.0;
//     else if (l == 1)
//       Pl = costheta;
//     else {
//       double Plm2 = 1.0;
//       double Plm1; 
//       if (= costheta;
//       for (int i=2; i<=l; i++) {
// 	double dl = (double)i;
// 	double dm = (double)m;
// 	Pl = (1.0/(dl-m))*((2.0*dl-1.0)*costheta*Plm1 - (dl+dm-1.0)*Plm2);
// 	Plm2 = Plm1;
// 	Plm1 = Pl;
//       }
//     }
//     complex<double> e2imphi (std::cos(m*phi), std::sin(m*phi));
//     return prefactor * Pl * e2imphi;
//   }



  void ECPComponentBuilder::CheckQuadratureRule(int lexact)
  {
    vector<PosType> &grid = pp_nonloc->sgridxyz_m;
    vector<RealType> &w = pp_nonloc->sgridweight_m;
    for (int l1=0; l1<=lexact; l1++) 
      for (int l2=0; l2 <= (lexact-l1); l2++) 
	for (int m1=-l1; m1<=l1; m1++)
	  for (int m2=-l2; m2<=l2; m2++) {
	    complex<double> sum(0.0, 0.0);
	    for (int k=0; k<grid.size(); k++) {
	      complex<double> v1 = Ylm(l1, m1, grid[k]);
	      complex<double> v2 = Ylm(l2, m2, grid[k]);
	      sum += 4.0*M_PI*w[k] * conj(v1)*v2;
	    }
	    double re = real (sum);
	    double im = imag (sum);
	    if ((l1==l2) && (m1==m2)) 
	      re -= 1.0;
	    if ((std::fabs(im) > 1.0e-14) || (std::fabs(re) > 1.0e-14)) {
	      app_error() << "Broken spherical quadrature for " << grid.size() << "-point rule.\n" << endl;
              OHMMS::Controller->abort();
	    }
// 	    fprintf (stderr, "(l1,m1,l2m,m2) = (%2d,%2d,%2d,%2d)  sum = (%20.16f %20.16f)\n",
// 	     l1, m1, l2, m2, real(sum), imag(sum));
	    
	  }
  }


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
        grid_local=createGrid(cur); 
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

  void ECPComponentBuilder::buildSemiLocalAndLocal(vector<xmlNodePtr>& semiPtr) {
    if(grid_inp==0)
    {
      app_error() << "  Pseudopotential file does not defined a global grid. vps/grid is disabled." << endl;
      OHMMS::Controller->abort();
    }

    //this is abinit/siesta format
    // There should only be one semilocal tag
    bool is_r_times_V(true);

    if (semiPtr.size()> 1) {
      ERRORMSG("We have more than one semilocal sections in the PP xml file.");
      OHMMS::Controller->abort();
    }

    RealType rmax = 0.0;
    //attributes: initailize by defaults
    string units("hartree");
    string format("r*V");
    string lloc;
    int ndown=1; 
    int nup=0;
    int nrule=4;//default quadrature
    int Llocal = -1;

    OhmmsAttributeSet aAttrib;
    aAttrib.add(units,"units");
    aAttrib.add(format,"format");
    aAttrib.add(ndown,"npots-down");
    aAttrib.add(nup,"npots-up");
    aAttrib.add(Llocal,"l-local");
    aAttrib.add(nrule,"nrule");

    xmlNodePtr cur_semilocal = semiPtr[0];
    aAttrib.put(cur_semilocal);

    RealType Vprefactor=1.0;
    if (units == "rydberg") 
      Vprefactor = 0.5;
    else if (units == "hartree")
      Vprefactor = 1.0;
    else {
      ERRORMSG("Unrecognized units """ << units << """ in PP file.");
      OHMMS::Controller->abort();
    }

    if (format == "r*V") 
      is_r_times_V = true;
    else if (format == "V")
      is_r_times_V = false;
    else {
      ERRORMSG("Unrecognized format """ << format << """ in PP file.");
      OHMMS::Controller->abort();
    }

    // Read which quadrature rule to use
    //const char *nrule = (const char*)xmlGetProp(cur_semilocal, (const xmlChar*)"nrule");
    //Nrule = (nrule==NULL) ? 4 : (int)atol(nrule);
    //cerr << "nrule = " << Nrule << endl;

    // We cannot construct the potentials as we construct them since
    // we may not know which one is local yet.
    typedef FSAtomPseudoPot<RealType> InFuncType;
    vector<InFuncType*> Vls;
    int iLocal=-1;
    // Now read vps sections
    xmlNodePtr cur_vps = cur_semilocal->children;
    while (cur_vps != NULL) {
      string vname ((const char*)cur_vps->name);
      if (vname == "vps") 
      {
        OhmmsAttributeSet aAttrib;
        string lstr("s");
        RealType rc=-1.0;
        aAttrib.add(lstr,"l");
        aAttrib.add(rc,"cutoff");
        aAttrib.put(cur_vps);

        InFuncType* avps=new InFuncType(angMon[lstr],rc,grid_inp);

        avps->put(cur_vps);
        Lmax=std::max(Lmax,avps->AngL);
        rmax=std::max(rmax,rc);

        //multiply by R
        if(!is_r_times_V) avps->convert2RV();

        Vls.push_back(avps);
      }
      cur_vps = cur_vps->next;
    }

    if (Llocal == -1) Llocal = Lmax;    

    Lmax=0;
    // Find which channel is the local one
    for (int i=0; i<Vls.size(); i++)
    {
      int l(Vls[i]->AngL);
      if (l == Llocal)
        iLocal = i;
      else
        Lmax = std::max(l,Lmax);
    }

    pp_loc = Vls[iLocal]->getLocalPot(-Vprefactor/Zeff);

    // Now construct the radial potentials
    NumNonLocal=0;
    for (int i=0; i<Vls.size(); i++) {
      if (i == iLocal) continue;
      RadialPotentialType* newpot = Vls[i]->getNonLocalPot(*Vls[iLocal],Vprefactor);
      pp_nonloc->add(Vls[i]->AngL,newpot);
      NumNonLocal++;
    }

    delete_iter(Vls.begin(),Vls.end());
    pp_nonloc->lmax=Lmax;
    pp_nonloc->Rmax=rmax;
  }

  ECPComponentBuilder::GridType* ECPComponentBuilder::createGrid(xmlNodePtr cur)
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

  bool ECPComponentBuilder::parseCasino(std::string& fname, RealType rmax)
  {
      ifstream fin(fname.c_str(),ios_base::in);
      if(!fin){
	ERRORMSG("Could not open file " << fname)
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

      LinearGrid<RealType> *agrid=new LinearGrid<RealType>;
      int ng=static_cast<int>(rmax/1.0e-3)+1;
      agrid->set(0,rmax,ng);

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
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

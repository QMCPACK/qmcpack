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
#include "QMCHamiltonians/ECPComponentBuilder.h"
#include "Numerics/GaussianTimesRN.h"
#include "Numerics/Transform2GridFunctor.h"
#include "QMCFactory/OneDimGridFactory.h"
#include "QMCHamiltonians/Ylm.h"
#include <cmath>

namespace qmcplusplus {

  ECPComponentBuilder::ECPComponentBuilder(const string& aname):
    NumNonLocal(0), Lmax(0), Zeff(0), Species(aname), Nrule(4),
    pp_loc(0), pp_nonloc(0) {
    angMon["s"]=0; angMon["p"]=1; angMon["d"]=2; angMon["f"]=3; angMon["g"]=4;
    angMon["0"]=0; angMon["1"]=1; angMon["2"]=2; angMon["3"]=3; angMon["4"]=4;
  }

  bool ECPComponentBuilder::parse(const string& fname) {

    // build an XML tree from a the file;
    xmlDocPtr m_doc = xmlParseFile(fname.c_str());
    if (m_doc == NULL) {
      ERRORMSG("File " << fname << " is invalid")
        xmlFreeDoc(m_doc);
    }    
    // Check the document is of the right kind
    xmlNodePtr cur = xmlDocGetRootElement(m_doc);
    if (cur == NULL) {
      ERRORMSG("Empty document");
      xmlFreeDoc(m_doc);
    }

    bool success=put(cur);
    xmlFreeDoc(m_doc);

    return success;
  }

  bool ECPComponentBuilder::put(xmlNodePtr cur) {
    int nk=0;
    vector<RealType> kpts;
    vector<xmlNodePtr> semiPtr;
    cur=cur->children;
    while(cur != NULL) {
      string cname((const char*)cur->name);
      if(cname == "header") {
        Zeff = atoi((const char*)xmlGetProp(cur,(const xmlChar*)"zval"));
      } else if(cname == "semilocal") {
        semiPtr.push_back(cur);//save the pointer
      } else if(cname == "local") {
        buildLocal(cur);
      } else if(cname == "sphericalGrid") {
        nk=atoi((const char*)xmlGetProp(cur,(const xmlChar*)"size"));
        kpts.resize(nk*4);
        putContent(kpts,cur);
      }
      cur=cur->next;
    }

    if(semiPtr.size()) {
      if(pp_nonloc==0) pp_nonloc=new NonLocalECPComponent; 
      if(pp_loc) {
        for(int i=0; i<semiPtr.size(); i++) addSemiLocal(semiPtr[i]);
      } 
      else 
        buildSemiLocalAndLocal(semiPtr);
 
      app_log() << "    Maximum angular momentum of NonLocalECP " << pp_nonloc->lmax << endl;
      app_log() << "    Maximum cutoff of NonLocalECP " << pp_nonloc->Rmax << endl;
    } 

    if(pp_nonloc) {
      pp_nonloc->lmax=Lmax;
      SetQuadratureRule(Nrule);
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
      double nrm = r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
      assert (std::fabs(nrm-1.0) < 1.0e-14);
      wSum += pp_nonloc->sgridweight_m[k];
//       fprintf (stderr, "%14.10f %14.10f %14.10f %14.10f\n",
// 	       pp_nonloc->sgridxyz_m[k][0],
// 	       pp_nonloc->sgridxyz_m[k][1],
// 	       pp_nonloc->sgridxyz_m[k][2],
// 	       pp_nonloc->sgridweight_m[k]);
    }
    assert (std::fabs(wSum - 1.0) < 1.0e-14);
    // Check the quadrature rule
    CheckQuadratureRule (lexact);

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
	      ERRORMSG("Broken spherical quadrature for " << grid.size()
		       << "-point rule.\n");
	      	      abort();
	    }
// 	    fprintf (stderr, "(l1,m1,l2m,m2) = (%2d,%2d,%2d,%2d)  sum = (%20.16f %20.16f)\n",
// 	     l1, m1, l2, m2, real(sum), imag(sum));
	    
	  }
    cerr << "Passed harmonic check for " << grid.size() << "-point quadrature rule.\n";
  }


  void ECPComponentBuilder::addSemiLocal(xmlNodePtr cur) {

    GridType* grid_semilocal=0;
    RealType rmax= pp_nonloc->Rmax;
    cur=cur->children;
    while(cur != NULL) {
      string cname((const char*)cur->name);
      if(cname == "grid") {
        grid_semilocal=OneDimGridFactory::createGrid(cur);
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
    RadialPotentialType *app=new RadialPotentialType(agrid);
    Transform2GridFunctor<InFuncType,RadialPotentialType> transform(a,*app);
    transform.generate(agrid->rmin(),agrid->rmax(), agrid->size());
    //int ig=0;
    //while(ig<agrid->size()) {
    //  (*app)(ig)=a.f((*agrid)[ig]);
    //  ++ig;
    //}
    //app->spline(0,0.0,app->size()-1,0.0);
    return app;
  }

  void ECPComponentBuilder::buildLocal(xmlNodePtr cur) {

    if(pp_loc) return; //something is wrong

    typedef GaussianTimesRN<RealType> InFuncType;
    GridType* grid_local=0;
    InFuncType vr;

    bool bareCoulomb=true;
    cur=cur->children;
    while(cur != NULL) {
      string cname((const char*)cur->name);
      if(cname == "grid") {
        grid_local=OneDimGridFactory::createGrid(cur);
      } else if(cname == "basisGroup") {
        vr.putBasisGroup(cur);
        bareCoulomb=false;
      }
      cur=cur->next;
    }

    int ng=grid_local->size();
    vector<RealType> v(ng);
    if(bareCoulomb) {
      for(int ig=0; ig<ng; ig++) {
        v[ig]=1.0/(*grid_local)[ig];
      }
    } else {
      RealType zinv=1.0/Zeff;
      for(int ig=0; ig<ng; ig++) {
        double r=(*grid_local)[ig];
        v[ig]=1.0/r-zinv*vr.f(r);
      }
    }
    //if(bareCoulomb) {
    //  for(int ig=0; ig<ng; ig++) {
    //    v[ig]=-Zeff/(*grid_local)[ig];
    //  }
    //} else {
    //  for(int ig=0; ig<ng; ig++) {
    //    double r=(*grid_local)[ig];
    //    v[ig]=vr.f(r)-Zeff/r;
    //  }
    //}

    pp_loc=new RadialPotentialType(grid_local,v);
    int imin = 0;
    RealType yprime_i = ((*pp_loc)(imin+1)-(*pp_loc)(imin))/pp_loc->dr(imin);
    pp_loc->spline(imin,yprime_i,ng-1,0.0);

    //cout << " Effective Z = " << Zeff << endl;
    //ofstream fout("C.local.dat");
    //fout.setf(std::ios::scientific, std::ios::floatfield);
    //fout.precision(12);
    //int ig=0;
    //while(ig<grid_semilocal->size()) {
    //  double r=(*grid_semilocal)[ig++];
    //  fout << setw(20) << r << setw(20) << localpp->f(r)-Zeff/r<< endl;
    //}
  }

  void ECPComponentBuilder::buildSemiLocalAndLocal(vector<xmlNodePtr>& semiPtr) {
    //this is abinit/siesta format
    // There should only be one semilocal tag
    bool is_r_times_V(false);

    if (semiPtr.size() != 1) {
      ERRORMSG("We have more than one semilocal sections in the PP xml file.");
      return;
    }
    RealType rmax = 0.0;
    xmlNodePtr cur_semilocal = semiPtr[0];

    // First, check the attributes of the semilocal section
    string units = (const char*) xmlGetProp(cur_semilocal, (const xmlChar*)"units");
    RealType Vprefactor;
    if (units == "rydberg")
      Vprefactor = 0.5;
    else if (units == "hartree")
      Vprefactor = 1.0;
    else {
      ERRORMSG("Unrecognized units """ << units << """ in PP file.");
      return;
    }

    // Read "format" type
    string format = (const char*)xmlGetProp(cur_semilocal, (const xmlChar*)"format");
    if (format == "r*V") 
      is_r_times_V = true;
    else if (format == "V")
      is_r_times_V = false;
    else {
      ERRORMSG("Unrecognized format """ << format << """ in PP file.");
      return;
    }

    // Read number of potentials
    int ndown = atol((const char*)xmlGetProp(cur_semilocal, (const xmlChar*)"npots-down"));
    int nup   = atol((const char*)xmlGetProp(cur_semilocal, (const xmlChar*)"npots-up"));
    
    // Read which channels is the local one
    string lloc = (const char *)xmlGetProp(cur_semilocal, (const xmlChar*)"l-local");
    if (lloc != "")
      Llocal = (int)atol(lloc.c_str());
    else
      Llocal = -1;

    // Read which quadrature rule to use
    const char *nrule = (const char*)xmlGetProp(cur_semilocal, (const xmlChar*)"nrule");
    Nrule = (nrule==NULL) ? 4 : (int)atol(nrule);
    cerr << "nrule = " << Nrule << endl;

    // Now read vps sections
    xmlNodePtr cur_vps = cur_semilocal->children;
    // We cannot construct the potentials as we construct them since
    // we may not know which one is local yet.
    vector<vector<double> > Vls;
    vector<GridType*> Grids;
    vector<int> Ls;
    while (cur_vps != NULL) {
      string vname ((const char*)cur_vps->name);
      if (vname == "vps") {
	int l=angMon[(const char*)xmlGetProp(cur_vps,(const xmlChar*)"l")];
	double rc = atof ((const char*)xmlGetProp(cur_vps, (const xmlChar*)"cutoff"));
	rmax = max(rmax, rc);

	xmlNodePtr cur_radfunc = cur_vps->children;
	if ((string)(const char*)cur_radfunc->name == "text")
	  cur_radfunc = cur_radfunc->next;
	assert ((string)(const char*)cur_radfunc->name == "radfunc");

	xmlNodePtr cur_data = cur_radfunc->children;
	GridType *grid_vps(NULL);
	vector<RealType> vr;
	bool addedPot = false;
	while (cur_data != NULL) {
	  if ((string) (const char*)cur_data->name == "grid")
	    grid_vps = OneDimGridFactory::createGrid(cur_data);
	  else if ((string) (const char*)cur_data->name == "data") {
	    // Read the numerical data
	    putContent(vr, cur_data);
	  }
	  else if ((string) (const char*)cur_data->name == "text")
	    ; // do nothing
	  else {
	    ERRORMSG( "Unrecognized section """ << cur_data->name << " in vps in reading PP.");
	    return;
	  }
	  cur_data = cur_data->next;
	}
	if ((vr.size()!=0) && (grid_vps != NULL) && !addedPot) {
	  if (is_r_times_V) {
	    for (int i=0; i<vr.size(); i++)
	      vr[i] /= (*grid_vps)[i];
	    // simple extrapolation
	    vr[0] = 2.0*vr[1] - vr[2];
	  }
	  for (int i=0; i<vr.size(); i++)
	    vr[i] *= Vprefactor;
	  Vls.push_back (vr);
	  Grids.push_back(grid_vps);
	  Ls.push_back(l);
	  addedPot = true;
	}
      }
      cur_vps = cur_vps->next;
    }

    if (Llocal == -1)
      Llocal = Lmax;    
    // Find which channel is the local one
    int iLocal = 0;
    for (int i=0; i<Ls.size(); i++)
      if (Ls[i] == Llocal)
	iLocal = i;
      else
	Lmax = std::max(Ls[i],Lmax);

    // Now construct the radial potentials
    for (int i=0; i<Vls.size(); i++) {
      if (i == iLocal) {  
	// Construct the local channel
	vector<double> vloc(Vls[i].size());
	for (int j=0; j<vloc.size(); j++)
	  vloc[j] = Vls[i][j]/(-Zeff);
	pp_loc = new RadialPotentialType (Grids[i], vloc);
	int imin = 0; int ng = vloc.size();
	RealType yprime_i = ((*pp_loc)(imin+1)-(*pp_loc)(imin))/pp_loc->dr(imin);
	pp_loc->spline(imin,yprime_i,ng-1,0.0);
	cerr << "Added local potential for l=" << Ls[i] << endl;
      }
      else {
	// Construct the nonlocal channels
	for (int j=0; j<Vls[i].size(); j++)
	  Vls[i][j] -= Vls[iLocal][j];
	RadialPotentialType *newpot = new RadialPotentialType(Grids[i], Vls[i]);
	int ng = Vls[i].size();
	RealType yprime_i = (Vls[i][1]-Vls[i][0])/Grids[i]->dr(0);	
	newpot->spline(0, yprime_i, ng-1, 0.0);
	pp_nonloc->add (Ls[i], newpot); 
	cerr << "Added nonlocal potential for l=" << Ls[i] << endl;
	NumNonLocal++;
      }
    }
    pp_nonloc->lmax=Lmax;
    pp_nonloc->Rmax=rmax;
  }

} // namespace qmcPlusPlus
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

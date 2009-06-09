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
#include "Message/CommOperators.h"
#include <cmath>

namespace qmcplusplus {

  ECPComponentBuilder::ECPComponentBuilder(const string& aname, Communicate* c):
    MPIObjectBase(c),
    RcutMax(-1), NumNonLocal(0), Lmax(0), Zeff(0), Species(aname), Nrule(4),
    grid_global(0),pp_loc(0), pp_nonloc(0) {
    angMon["s"]=0; angMon["p"]=1; angMon["d"]=2; angMon["f"]=3; angMon["g"]=4;
    angMon["0"]=0; angMon["1"]=1; angMon["2"]=2; angMon["3"]=3; angMon["4"]=4;
  }

  bool ECPComponentBuilder::parse(const string& fname, xmlNodePtr cur) 
  {
    const xmlChar* rptr=xmlGetProp(cur,(const xmlChar*)"cutoff");
    if(rptr != NULL) RcutMax = atof((const char*)rptr);

    int length=0;
    char* cbuffer=0;

    ifstream *fin=0;
    if(myComm->rank()==0)
    {
      fin = new ifstream(fname.c_str());
      // if (!fin->is_open()) 
      // 	APP_ABORT("Could not open pseudopotential file \"" + fname + "\".\n");
      fin->seekg (0, ios::end);
      length = fin->tellg();
      fin->seekg (0, ios::beg);
    }

    myComm->bcast(length);
    cbuffer = new char[length];
    if(myComm->rank()==0) fin->read (cbuffer,length);
    myComm->bcast(cbuffer,length);
    xmlDocPtr m_doc = xmlReadMemory(cbuffer,length,NULL,NULL,0);

    if(fin)  delete fin;
    if(cbuffer) delete [] cbuffer;

    // build an XML tree from a the file;
    //xmlDocPtr m_doc = xmlParseFile(fname.c_str());
    if (m_doc == NULL) {
      xmlFreeDoc(m_doc);
      APP_ABORT("ECPComponentBuilder::parse xml file "+fname+" is invalid");
    }    
    // Check the document is of the right kind
    cur = xmlDocGetRootElement(m_doc);
    if (cur == NULL) {
      xmlFreeDoc(m_doc);
      APP_ABORT("Empty document");
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
        //capture the global grid
        grid_global = createGrid(cur);
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
      app_log() << "    Non-local pseudopotential parameters" <<endl;
      pp_nonloc->print(app_log());
      app_log() << "    Maximum cutoff radius " << pp_nonloc->Rmax << endl;

    }
    return true;
  }

  void ECPComponentBuilder::printECPTable() 
  {
    char fname[12];
    sprintf(fname,"%s.pp.dat",Species.c_str());
    ofstream fout(fname);
    fout.setf(std::ios::scientific, std::ios::floatfield);
    fout.precision(12);
    int nl=pp_nonloc?pp_nonloc->nlpp_m.size():0;
    RealType d=1.7e-2;
    RealType rt=0.13*d;
    if(nl)
    {
      fout << "#  Lmax = " << Lmax+1 << " nonlocal L channels" << nl << endl;
      fout << "#  Units = bohr hartree " << endl;
      fout << "#  r  -r*V/zeff   Vnl ... " << endl;
      while(rt<5)
      {
        fout << rt << setw(25) << pp_loc->splint(rt);
        for(int l=0; l<nl; l++)
          fout << setw(25) << pp_nonloc->nlpp_m[l]->splint(rt);
        fout << endl; 
        rt+=d;
      }
    }
    else
    {
      fout << "#  Units = bohr hartree " << endl;
      fout << "#  r  -r*V/zeff " << endl;
      while(rt<5)
      {
        fout << rt << setw(25) << pp_loc->splint(rt) << endl;;
        rt+=d;
      }
    }
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
      nk = 50;  symmetry = OCTA;   lexact = 11;
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
              APP_ABORT("Give up");
	    }
// 	    fprintf (stderr, "(l1,m1,l2m,m2) = (%2d,%2d,%2d,%2d)  sum = (%20.16f %20.16f)\n",
// 	     l1, m1, l2, m2, real(sum), imag(sum));
	    
	  }
  }

} // namespace qmcPlusPlus
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

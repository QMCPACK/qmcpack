//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#ifndef QMCPLUSPLUS_BESSELZERO_FUNCTOR_H
#define QMCPLUSPLUS_BESSELZERO_FUNCTOR_H
#include "Numerics/OptimizableFunctorBase.h"
#include "Utilities/ProgressReportEngine.h"
#include "OhmmsData/AttributeSet.h"
#include <cstdio>

namespace qmcplusplus
{

template<class T>
struct BesselZero: public OptimizableFunctorBase
{

  typedef real_type value_type;
  int NumParams;
  std::vector<real_type> Parameters;
  real_type R_B,R_Binv;
  real_type C_0,C_0inv;
  std::string pr;

  ///constructor
  BesselZero():R_B(1.0)
  {
  }

  OptimizableFunctorBase* makeClone() const
  {
    return new BesselZero(*this);
  }

  void resize(int n)
  {
    NumParams = n;
    Parameters.resize(n,0);
  }

  void reset()
  {
    if (pr=="yes")
      print();
  }

  inline real_type evaluate(real_type r)
  {
    real_type u=0;
    real_type A(C_0);
    real_type Ainv(C_0inv);
    real_type AX(A*r);
    real_type Xinv(1.0/r);
    real_type AXinv(Ainv*Xinv);
    for(int i=0; i<NumParams; i++)
    {
      real_type j(i+1);
      real_type siniAX = std::sin(j*AX);
      real_type iinv(1.0/j);
      u += Parameters[i]*AXinv*iinv*siniAX;
    }
    return u;
  }

  inline real_type evaluate(real_type r, real_type rinv)
  {
    return evaluate(r);
  }

  inline void evaluateAll(real_type r, real_type rinv)
  {
    real_type du(0),d2u(0);
    evaluate(r,du,d2u);
  }

  inline real_type
  evaluate(real_type r, real_type& dudr, real_type& d2udr2)
  {
    real_type u(0);
    sumB(r,u,dudr,d2udr2);
    return u;
  }


  inline real_type
  evaluate(real_type r, real_type& dudr, real_type& d2udr2, real_type &d3udr3)
  {
    real_type u(0);
    sumB(r,u,dudr,d2udr2);
    return u;
  }

  inline void sumB(real_type X, real_type& u, real_type& du, real_type& d2u)
  {
    u=0;
    du=0;
    d2u=0;
    real_type A(C_0);
    real_type Ainv(C_0inv);
    real_type AX(A*X);
    real_type AX2(AX*AX);
    real_type Xinv(1.0/X);
    real_type AXinv(Ainv*Xinv);
    real_type AX2inv(Ainv*Xinv*Xinv);
    real_type AX3inv(Ainv*Xinv*Xinv*Xinv);
    for(int i=0; i<NumParams; i++)
    {
      real_type j(i+1.0);
      real_type cosjAX = std::cos(j*AX);
      real_type sinjAX = std::sin(j*AX);
      real_type jinv(1.0/j);
      u += Parameters[i]*AXinv*jinv*sinjAX;
      du += Parameters[i]*(j*AX*cosjAX - sinjAX)*AX2inv*jinv;
      d2u += Parameters[i]*(-2.0*j*AX*cosjAX + (2.0-j*j*AX2)*sinjAX)*AX3inv*jinv;
    }
  }


  inline bool
  evaluateDerivatives (real_type r, std::vector<TinyVector<real_type,3> >& derivs)
  {
    real_type A(C_0);
    real_type Ainv(C_0inv);
    real_type X(r);
    real_type Xinv(1.0/r);
    real_type AX(A*X);
    real_type AX2(AX*AX);
    real_type AXinv(Ainv*Xinv);
    real_type AX2inv(Ainv*Xinv*Xinv);
    real_type AX3inv(Ainv*Xinv*Xinv*Xinv);
    for(int i=0; i<NumParams; i++)
    {
      real_type j(i+1.0);
      real_type cosjAX = std::cos(j*AX);
      real_type sinjAX = std::sin(j*AX);
      real_type jinv(1.0/j);
      derivs[i][0] += AXinv*jinv*sinjAX;
      derivs[i][1] += (j*AX*cosjAX - sinjAX)*AX2inv*jinv;
      derivs[i][2] += (-2.0*j*AX*cosjAX + (2.0-j*j*AX2)*sinjAX)*AX3inv*jinv;
    }
    return true;
  }

  inline real_type f(real_type r)
  {
    real_type du, d2u;
    return evaluate (r, du, d2u);
  }

  inline real_type df(real_type r)
  {
    real_type du, d2u;
    evaluate (r, du, d2u);
    return du;
  }

  bool put(xmlNodePtr cur)
  {
    ReportEngine PRE("BesselZero","put(xmlNodePtr)");
    //CuspValue = -1.0e10;
    NumParams = 0;
    pr = "no";
    OhmmsAttributeSet rAttrib;
    rAttrib.add(NumParams,   "size");
    rAttrib.add(R_B,   "RB");
    rAttrib.add(pr,   "print");
    rAttrib.put(cur);
    app_log()<<" R_B is set to: "<<R_B<< std::endl;
    R_Binv = 1.0/R_B;
    C_0 = 3.1415926535897932384626433832795028841968*R_Binv;
    C_0inv = 1.0/C_0;
    if (NumParams == 0)
    {
      PRE.error("You must specify a positive number of parameters for the BesselZero jastrow function.",true);
    }
    app_log() << " size = " << NumParams << " parameters " << std::endl;
    resize (NumParams);
    // Now read coefficents
    xmlNodePtr xmlCoefs = cur->xmlChildrenNode;
    while (xmlCoefs != NULL)
    {
      std::string cname((const char*)xmlCoefs->name);
      if (cname == "coefficients")
      {
        std::string type("0"), id("0");
        OhmmsAttributeSet cAttrib;
        cAttrib.add(id, "id");
        cAttrib.add(type, "type");
        cAttrib.put(xmlCoefs);
        if (type != "Array")
        {
          PRE.error( "Unknown correlation type " + type + " in BesselZero." + "Resetting to \"Array\"");
          xmlNewProp (xmlCoefs, (const xmlChar*) "type", (const xmlChar*) "Array");
        }
        std::vector<real_type> params;
        putContent(params, xmlCoefs);
        if (params.size() == NumParams)
          Parameters = params;
        // Setup parameter names
        for (int i=0; i< NumParams; i++)
        {
          std::stringstream sstr;
          sstr << id << "_" << i;
          myVars.insert(sstr.str(),Parameters[i],true);
        }
        app_log() << "Parameter     Name      Value\n";
        myVars.print(app_log());
      }
      xmlCoefs = xmlCoefs->next;
    }
    reset();
    return true;
  }

  void checkInVariables(opt_variables_type& active)
  {
    active.insertFrom(myVars);
  }

  void checkOutVariables(const opt_variables_type& active)
  {
    myVars.getIndex(active);
  }

  void resetParameters(const opt_variables_type& active)
  {
    for(int i=0; i<Parameters.size(); ++i)
    {
      int loc=myVars.where(i);
      if(loc>=0)
        Parameters[i]=myVars[i]=active[loc];
    }
    reset();
  }


  void print()
  {
//       std::string fname = (elementType != "") ? elementType : pairType;
    std::string fname = "BesselZero.dat";
//       //cerr << "Writing " << fname << " file.\n";
    FILE *fout = fopen (fname.c_str(), "w");
    for (double r=0.0; r<R_B*4; r+=0.01)
      fprintf (fout, "%8.3f %16.10f\n", r, evaluate(r));
    fclose(fout);
  }


  void print(std::ostream& os)
  {
//       int n=100;
//       T d=cutoff_radius/100.,r=0;
//       T u,du,d2du;
//       for(int i=0; i<n; ++i)
//       {
//         u=evaluate(r,du,d2du);
//         os << std::setw(22) << r << std::setw(22) << u << std::setw(22) << du
//           << std::setw(22) << d2du << std::endl;
//         r+=d;
//       }
  }
};
}
#endif

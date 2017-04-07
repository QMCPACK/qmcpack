//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jordan E. Vincent, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jordan E. Vincent, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#ifndef QMCPLUSPLUS_CHEBYSHEV_JASTROWFUNCTION_H
#define QMCPLUSPLUS_CHEBYSHEV_JASTROWFUNCTION_H
#include "OhmmsData/libxmldefs.h"
#include "Optimize/VarList.h"

/**class Pade functional
 *@brief \f[ u(r) = \frac{a*r}{1+b*r} \f]
 * Prototype of the template parameter of TwoBodyJastrow and OneBodyJastrow
 */
template<class T>
struct ChebyshevFunctor
{

  ///coefficients
  T B, L;
  int lmax;
  std::vector<T> alpha;

  std::vector<T> TT, dTT, d2TT;
  ///constructor
  ChebyshevFunctor() { }

  /**
   *@brief reset the internal variables.
   */
  inline void reset() { }

  void reset(std::vector<T> a)
  {
    alpha = a;
  }

  /**@param r the distance
     @param dudr return value  \f$ du/dr = a/(1+br)^2 \f$
     @param d2udr2 return value  \f$ d^2u/dr^2 = -2ab/(1+br)^3 \f$
     @return \f$ u(r) = a*r/(1+b*r) \f$
  */
  inline T evaluate(T r, T& dudr, T& d2udr2)
  {
    T rbar = (2.0*r-L)/L;
    T drbardr = 2.0/L;
    TT[0] = 1.0;
    dTT[0] = 0.0;
    d2TT[0] = 0.0;
    TT[1] = rbar;
    dTT[1] = r;
    d2TT[1] = 0.0;
    T Tsum = alpha[0]*TT[0] + alpha[1]*TT[1];
    T dTsum = alpha[1]*dTT[1]*drbardr;
    T d2Tsum = 0.0;
    for(int l=1; l<=lmax-1; l++)
    {
      TT[l+1] = 2*rbar*TT[l]-TT[l-1];
      Tsum += alpha[l+1]*TT[l+1];
      dTT[l+1] = 2*TT[l]+2*r*dTT[l]-d2TT[l-1];
      dTsum += alpha[l+1]*dTT[l+1]*drbardr;
      d2TT[l+1] = 4*dTT[l]+2*r*d2TT[l]-d2TT[l-1];
      d2Tsum += alpha[l+1]*d2TT[l+1]*drbardr*drbardr;
    }
    T rsq = r*r;
    T r_minus_L = r-L;
    T r_minus_Lsq = r_minus_L*r_minus_L;
    T r_half_L = 0.5*L+r;
    T fact1 = r_minus_Lsq*rsq;
    T fact2 = 2.0*(r_minus_L*rsq+r_minus_Lsq*r);
    T fact3 = 2.0*(rsq+4.0*r*r_minus_L+r_minus_Lsq);
    T u = fact1*Tsum + B*r_minus_Lsq*r_half_L;
    dudr = fact2*Tsum + fact1*dTsum + 2.0*B*r_minus_L*r_half_L+B*r_minus_Lsq;
    d2udr2 = fact3*Tsum + 2.0*fact2*dTsum +
             2.0*fact1*dTsum+fact1*d2Tsum + 2.0*B*r_half_L+4.0*B*r_minus_L;
    return u;
  }

  /**@param cur current xmlNode from which the data members are reset
     @param vlist VarRegistry<T1> to which the Pade variables A and B
     are added for optimization
     @brief T1 is the type of VarRegistry, typically double.  Read
     in the Pade parameters from the xml input file.
  */
  template<class T1>
  void put(xmlNodePtr cur, VarRegistry<T1>& vlist)
  {
    std::string idalpha, idb;
    //jastrow[iab]->put(cur->xmlChildrenNode,wfs_ref.RealVars);
    xmlNodePtr tcur = cur->xmlChildrenNode;
    while(tcur != NULL)
    {
      //@todo Var -> <param(eter) role="opt"/>
      std::string cname((const char*)(tcur->name));
      if(cname == "parameter" || cname == "Var")
      {
        std::string aname((const char*)(xmlGetProp(tcur,(const xmlChar *)"name")));
        std::string idname((const char*)(xmlGetProp(tcur,(const xmlChar *)"id")));
        if(aname == "alpha")
        {
          idalpha = idname;
          putContent(alpha,tcur);
        }
        else
          if(aname == "L")
          {
            putContent(L,tcur);
          }
          else
            if(aname == "lmax")
            {
              putContent(lmax,tcur);
            }
            else
              if(aname == "B")
              {
                idb = idname;
                putContent(B,tcur);
              }
      }
      tcur = tcur->next;
    }
    vlist.add(idalpha,&alpha[0],alpha.size());
    vlist.add(idb,&B,1);
    XMLReport("Jastrow Parameters = (")
    for(int i=0; i<alpha.size(); i++)
      XMLReport(alpha[i] << ", ")
      XMLReport(")")
    }
};
#endif


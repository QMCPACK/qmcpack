//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
    
    
#ifndef QMCPLUSPLUS_GAUSSIANFUNCTOR
#define QMCPLUSPLUS_GAUSSIANFUNCTOR
#include "Numerics/OptimizableFunctorBase.h"
#include "QMCWaveFunctions/MolecularOrbitals/NGOBuilder.h"
#include <cmath>
#include <vector>

namespace qmcplusplus
{
/**  Linear combination of GTO.
 *   Actual evaluation is made from a NGO object build from the
 *   gaussian expansion.
 *   f(r) = \sum_i alp_i * exp(-expo_i*(r-R0_i)**2)
 */
//template<class T>
class GaussianFunctor: public OptimizableFunctorBase
{

public:

  int NumFuns;

  bool IsOptimizing;

  std::vector<real_type> R0;
  std::vector<real_type> alp;
  std::vector<real_type> expo;

  std::vector<std::string> R0_id;
  std::vector<std::string> alp_id;
  std::vector<std::string> expo_id;

  std::vector<bool> activeOptm;

  NGOrbital *NumOrb;

  GaussianFunctor():NumFuns(0),IsOptimizing(false) {}

  ~GaussianFunctor() {};

  void checkOutVariables(const opt_variables_type& active)
  {
    if(IsOptimizing)
      myVars.getIndex(active);
  }

  void checkInVariables(opt_variables_type& active)
  {
    if(IsOptimizing)
      active.insertFrom(myVars);
  }

  /** reset the optimizable variables
   * @param active list of active optimizable variables
   */
  void resetParameters(const opt_variables_type& active)
  {
    if(!IsOptimizing)
      return;
    int n=0,m=0;
    for(int i=0; i<NumFuns; i++)
    {
      if (activeOptm[n])
      {
        int ia=myVars.where(m);
        if(ia>-1)
          R0[i]=myVars[m++]=active[ia];
      }
      n++;
      if (activeOptm[n])
      {
        int ia=myVars.where(m);
        if(ia>-1)
          alp[i]=myVars[m++]=active[ia];
      }
      n++;
      if (activeOptm[n])
      {
        int ia=myVars.where(m);
        if(ia>-1)
          expo[i]=myVars[m++]=active[ia];
      }
      n++;
    }
  }

  /** reset function
   */
  void reset()
  {
    if(!IsOptimizing)
      return;
  }

  /** evaluate the value at r
   * @param r distance
   *
   * virtual function necessary for a transformation to a numerical functor
   */
  real_type f(real_type r)
  {
    real_type res=0;
    for(int i=0; i<NumFuns; i++)
    {
      real_type dr=r-R0[i];
      res += alp[i]*std::exp(-expo[i]*dr*dr);
    }
    return res;
  }

  /** evaluate the first derivate
   * @param r distance
   *
   * virtual function necessary for a transformation to a numerical functor
   */
  real_type df(real_type r)
  {
    real_type du=0.0;
    for(int i=0; i<NumFuns; i++)
    {
      real_type dr=r-R0[i];
      du -= 2.0*expo[i]*dr*alp[i]*std::exp(-expo[i]*dr*dr);
    }
    return du;
  }

  /** process xmlnode and registers variables to optimize
   * @param cur xmlNode for a functor
   */
  bool put(xmlNodePtr cur)
  {
    app_log() <<"GaussianFunctor::put() " << std::endl;
    real_type r(0.0),a(0.1),expon(1.0);
    cur = cur->xmlChildrenNode;
    bool renewed=false;
    while(cur != NULL)
    {
      std::string cname((const char*)(cur->name));
      if(cname == "parameter" || cname == "Var")
      {
        std::string id_in("BF");
        std::string optm("no");
        OhmmsAttributeSet rAttrib;
        rAttrib.add(id_in, "id");
        rAttrib.add(r, "r0");
        rAttrib.add(expon, "expo");
        rAttrib.add(a, "a");
        rAttrib.add(optm, "optimize");
        rAttrib.put(cur);
        R0_id.push_back(id_in+std::string("_R0"));
        alp_id.push_back(id_in+std::string("_alp"));
        expo_id.push_back(id_in+std::string("_expo"));
        R0.push_back(r);
        alp.push_back(a);
        expo.push_back(expon);
        NumFuns++;
        // optm R0 ?  disable for now
//         activeOptm.push_back(optm=="all" || optm=="yes" || optm=="1" || optm == "4" || optm =="5");
        activeOptm.push_back(false);
        // optm alp ?
        activeOptm.push_back(optm=="all" || optm=="yes" || optm=="2" || optm == "4" || optm =="6");
        // optm expo ?
        activeOptm.push_back(optm=="all" || optm=="yes" || optm=="3" || optm == "5" || optm =="6");
        if(activeOptm[activeOptm.size()-3])
          IsOptimizing=true;
        if(activeOptm[activeOptm.size()-2])
          IsOptimizing=true;
        if(activeOptm[activeOptm.size()-1])
          IsOptimizing=true;
        app_log() <<"Adding Gaussian funtion to Backflow Transformation.\n";
        app_log() <<"r0,alp,expon,optim_r0,optim_alp,optim_expon:"
                  <<a <<"  "
                  <<expon <<"  "
                  <<r <<"\n"
                  <<activeOptm[activeOptm.size()-3] <<" "
                  <<activeOptm[activeOptm.size()-2] <<" "
                  <<activeOptm[activeOptm.size()-1] << std::endl;
        if(activeOptm[activeOptm.size()-3])
          myVars.insert(R0_id.back(),R0.back(),true);
        if(activeOptm[activeOptm.size()-2])
          myVars.insert(alp_id.back(),alp.back(),true);
        if(activeOptm[activeOptm.size()-1])
          myVars.insert(expo_id.back(),expo.back(),true);
      }
      cur = cur->next;
    }
    return true;
  }

  /** empty virtual function to help builder classes
  */
  void setDensity(real_type n) { }

  inline bool evaluateDerivatives (real_type r, std::vector<qmcplusplus::TinyVector<real_type,3> >& derivs)
  {
    return false;
  }

  OptimizableFunctorBase* makeClone() const
  {
    GaussianFunctor* clone =  new GaussianFunctor();
    clone->NumFuns=NumFuns;
    clone->R0=R0;
    clone->alp=alp;
    clone->expo=expo;
    clone->R0_id=R0_id;
    clone->alp_id=alp_id;
    clone->expo_id=expo_id;
    clone->IsOptimizing=IsOptimizing;
    // NumOrb = fn.NumOrb->makeClone();
    return clone;
  }

  void add(real_type pos, real_type c, real_type exponent)
  {
    R0.push_back(pos);
    alp.push_back(c);
    expo.push_back(exponent);
    NumFuns++;
    //if(R0.size() != NumFuns)
    //  APP_ABORT("R0.size() != NumFuns");
  }

  void createNGO() {}

  inline real_type evaluate(const double r)
  {
    real_type res=0;
    for(int i=0; i<NumFuns; i++)
    {
      double dr=r-R0[i];
      res += alp[i]*std::exp(-expo[i]*dr*dr);
    }
    return res;
  }

  inline real_type evaluate(const double r, real_type& du, real_type& d2u)
  {
    real_type res=0;
    du=d2u=0.0;
    for(int i=0; i<NumFuns; i++)
    {
      double dr=r-R0[i];
      real_type u = alp[i]*std::exp(-expo[i]*dr*dr);
      real_type tmp = 2.0*expo[i]*dr;
      res+=u;
      du -= tmp*u;
      d2u += (-2.0*expo[i] + tmp*tmp)*u;
    }
    return res;
  }

};

}

#endif

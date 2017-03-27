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
    
    
#ifndef QMCPLUSPLUS_MCMILLANJ2_H
#define QMCPLUSPLUS_MCMILLANJ2_H
#include "Numerics/OptimizableFunctorBase.h"

namespace qmcplusplus
{
/** ModPade Jastrow functional
 *
 * \f[ u(r) = \frac{1}{2*A}\left[1-\exp(-A*r)\right], \f]
 */
template<class T>
struct McMillanJ2Functor: public OptimizableFunctorBase
{

  typedef typename OptimizableFunctorBase::real_type real_type;

  ///coefficients
  real_type A;
  real_type B;

  std::string ID_A;
  std::string ID_B;

  /** constructor
  * @param a A coefficient
  * @param samespin boolean to indicate if this function is for parallel spins
   */
  McMillanJ2Functor(real_type a=5.0, real_type b=4.9133):ID_A("McA"),ID_B("McB")
  {
    reset(a,b);
  }

  OptimizableFunctorBase* makeClone() const
  {
    return new McMillanJ2Functor<T>(*this);
  }
  inline void reset()
  {
    reset(A,B);
  }

  /** reset the internal variables.
  *
  * USE_resetParameters
   */
  void resetParameters(const opt_variables_type& active)
  {
    int ia=myVars.where(0);
    if(ia>-1)
      A=active[ia];
    int ib=myVars.where(1);
    if(ib>-1)
      B=active[ib];
    reset(A,B);
  }

  void checkInVariables(opt_variables_type& active)
  {
    active.insertFrom(myVars);
  }

  void checkOutVariables(const opt_variables_type& active)
  {
    myVars.getIndex(active);
  }


  /** reset the internal variables.
  *@param a New Jastrow parameter a
   */
  inline void reset(real_type a, real_type b)
  {
    A = a;
    B = b;
  }

  /** evaluate the value at r
  * @param r the distance
  * @return \f$ u(r) = \frac{A}{r}\left[1-\exp(-\frac{r}{F})\right]\f$
   */
  inline real_type evaluate(real_type r)
  {
    return std::pow(B/r,A);
  }

  /**@param r the distance
  @param dudr first derivative
  @param d2udr2 second derivative
  @return the value
   */
  inline real_type evaluate(real_type r, real_type& dudr, real_type& d2udr2)
  {
    real_type wert = std::pow(B/r,A);
    dudr = -A*wert/r;
    d2udr2= (A+1.)*A*wert/(r*r);
    return wert;
  }

  /** return a value at r
   */
  real_type f(real_type r)
  {
    return evaluate(r);
  }

  /** return a derivative at r
   */
  real_type df(real_type r)
  {
    return -A*std::pow(B/r,A)/r;
  }

  /** Read in the parameter from the xml input file.
  * @param cur current xmlNode from which the data members are reset
   */
  bool put(xmlNodePtr cur)
  {
    xmlNodePtr tcur = cur->xmlChildrenNode;
    while(tcur != NULL)
    {
      //@todo Var -> <param(eter) role="opt"/>
      std::string cname((const char*)(tcur->name));
      if(cname == "parameter" || cname == "Var")
      {
        std::string aname((const char*)(xmlGetProp(tcur,(const xmlChar *)"name")));
//            std::string idname((const char*)(xmlGetProp(tcur,(const xmlChar *)"id")));
        if(aname == "a")
        {
          putContent(A,tcur);
        }
        else
          if(aname == "b")
          {
            ID_B = (const char*)(xmlGetProp(tcur,(const xmlChar *)"id"));
            putContent(B,tcur);
          }
      }
      tcur = tcur->next;
    }
    myVars.insert(ID_A,A);
    myVars.insert(ID_B,B);
    reset(A,B);
    return true;
  }


};

//   template<class T>
//       struct ModMcMillanJ2Functor: public OptimizableFunctorBase
//   {
//
//     typedef typename OptimizableFunctorBase::real_type real_type;
//
//       ///coefficients
//     real_type A;
//     real_type B;
//     real_type RC;
//     real_type c0,c1,c2,c3,c4,c5,c6;
//
//     std::string ID_A,ID_B,ID_RC;
//
//       /** constructor
//      * @param a A coefficient
//      * @param samespin boolean to indicate if this function is for parallel spins
//        */
//     ModMcMillanJ2Functor(real_type a=4.9133, real_type b=5, real_type rc=1.0):ID_RC("Rcutoff"),ID_A("A"),ID_B("B")
//     {
//       A = a;
//       B = b;
//       RC = rc;
//       reset();
//     }
//
//
//     void checkInVariables(opt_variables_type& active)
//     {
//       active.insertFrom(myVars);
//     }
//
//     void checkOutVariables(const opt_variables_type& active)
//     {
//       myVars.getIndex(active);
//     }
//
//     OptimizableFunctorBase* makeClone() const
//     {
//       return new ModMcMillanJ2Functor<T>(*this);
//     }
//
//     void resetParameters(const opt_variables_type& active)
//     {
//       int ia=myVars.where(0); if(ia>-1) A=active[ia];
//       int ib=myVars.where(1); if(ib>-1) B=active[ib];
//       int ic=myVars.where(2); if(ic>-1) RC=active[ic];
//       reset();
//     }
//
//     inline void reset() {
//       c0 = std::pow(A,B);
//       c1 = -1.0*std::pow(A/RC,B);
//       c2 =  B*c0*std::pow(RC,-B-1);
//       c3 =  -0.5*(B+1)*B*c0*std::pow(RC,-B-2);
//
//       c4 = -B*c0;
//       c5 = 2.0*c3;
//       c6 = B*(B+1.0)*c0;
//     }
//
//
//       /** evaluate the value at r
//      * @param r the distance
//      * @return \f$ u(r) = \frac{A}{r}\left[1-\exp(-\frac{r}{F})\right]\f$
//        */
//     inline real_type evaluate(real_type r) {
//       if (r<RC) {
//         real_type r1 = (r-RC);
//         return c0*std::pow(r,-B)+c1 + r1*(c2+r1*c3) ;
//       }
//       else {return 0.0;};
//     }
//
//       /**@param r the distance
//     @param dudr first derivative
//     @param d2udr2 second derivative
//     @return the value
//        */
//     inline real_type evaluate(real_type r, real_type& dudr, real_type& d2udr2) {
//       if (r<RC){
//         real_type r1 = (r-RC);
//         real_type overr = 1.0/r;
//         real_type wert = std::pow(1.0/r,B);
//         dudr = c4*wert*overr + c2 + c5*r1;
//         d2udr2= c6*wert*overr*overr+ c5 ;
//         return c0*wert + c1 + r1*(c2+r1*c3) ;
//       } else {
//         dudr = 0.0;
//         d2udr2= 0.0;
//         return 0.0;
//       };
//     }
//
//       /** return a value at r
//        */
//     real_type f(real_type r) {
//       return evaluate(r);
//     }
//
//       /** return a derivative at r
//        */
//     real_type df(real_type r) {
//       if (r<RC) {return c4*std::pow(1.0/r,B+1) + c2 + c5*(r-RC) ;}
//       else {return 0.0;};
//     }
//
//       /** Read in the parameter from the xml input file.
//      * @param cur current xmlNode from which the data members are reset
//        */
//     bool put(xmlNodePtr cur) {
//       RC = cutoff_radius;
//       xmlNodePtr tcur = cur->xmlChildrenNode;
//       while(tcur != NULL) {
//           //@todo Var -> <param(eter) role="opt"/>
//         std::string cname((const char*)(tcur->name));
//         if(cname == "parameter" || cname == "Var") {
//           std::string aname((const char*)(xmlGetProp(tcur,(const xmlChar *)"name")));
//           std::string idname((const char*)(xmlGetProp(tcur,(const xmlChar *)"id")));
//           if(aname == "A") {
//             ID_A = idname;
//             putContent(A,tcur);
//           } else if(aname == "B") {
//             ID_B = idname;
//             putContent(B,tcur);
//           }else if(aname == "RC")
//           {
//             ID_RC = idname;
//             putContent(RC,tcur);
//           }
//         }
//         tcur = tcur->next;
//       }
//       myVars.insert(ID_A,A);
//       myVars.insert(ID_B,B);
//       myVars.insert(ID_RC,RC);
//       reset();
//       LOGMSG("  Modified McMillan Parameters ")
//           LOGMSG("    A (" << ID_A << ") = " << A  << "  B (" << ID_B << ") =  " << B<< "  R_c (" << ID_RC << ") =  " << RC)
//       return true;
//     }
//
//   };

template<class T>
struct ModMcMillanJ2Functor: public OptimizableFunctorBase
{

  ///f(r) = A*r^-5 - A*r_c^-5 + A/(B*r_c^6) * tanh(B*(r-r_c))
  typedef typename OptimizableFunctorBase::real_type real_type;

  ///coefficients
  real_type A;
  real_type B;
  real_type RC;
  real_type cA,c0,c1,c2,c3,c4,c5;

  std::string ID_A,ID_B,ID_RC;

  /** constructor
  * @param a A coefficient
  * @param samespin boolean to indicate if this function is for parallel spins
   */
  ModMcMillanJ2Functor(real_type a=0.0, real_type b=0.0, real_type rc=0.0):ID_RC("Rcutoff"),ID_A("A"),ID_B("B")
  {
  }


  void checkInVariables(opt_variables_type& active)
  {
    active.insertFrom(myVars);
  }

  void checkOutVariables(const opt_variables_type& active)
  {
    myVars.getIndex(active);
  }

  OptimizableFunctorBase* makeClone() const
  {
    return new ModMcMillanJ2Functor<T>(*this);
  }

  void resetParameters(const opt_variables_type& active)
  {
    int ia=myVars.where(0);
    if(ia>-1)
      A=active[ia];
    int ib=myVars.where(1);
    if(ib>-1)
      B=active[ib];
    int ic=myVars.where(2);
    if(ic>-1)
      RC=active[ic];
    reset();
  }

  inline void reset()
  {
    cA = std::pow(A,5);
    c0 = cA*std::pow(RC,-5);
    c1 = 5*cA*std::pow(RC,-6)/B;
    c2 = -5*cA;
    c3 = B*c1;
    c4 = 30*cA;
    c5 = -2*B*c3;
  }


  /** evaluate the value at r
  * @param r the distance
  * @return \f$ u(r) = \frac{A}{r}\left[1-\exp(-\frac{r}{F})\right]\f$
   */
  inline real_type evaluate(real_type r)
  {
    if (r<RC)
    {
      real_type r1 = (r-RC);
      real_type plZ = std::exp(B*r1);
      real_type miZ = 1.0/plZ;
      real_type tanhZ = (plZ-miZ)/(plZ+miZ);
      return cA*std::pow(r,-5)-c0 + c1*tanhZ;
    }
    else
    {
      return 0.0;
    };
  }

  /**@param r the distance
  @param dudr first derivative
  @param d2udr2 second derivative
  @return the value
   */
  inline real_type evaluate(real_type r, real_type& dudr, real_type& d2udr2)
  {
    if (r<RC)
    {
      real_type prt = cA*std::pow(r,-5);
      real_type r1 = (r-RC);
      real_type plZ = std::exp(B*r1);
      real_type miZ = 1.0/plZ;
      real_type tanhZ = (plZ-miZ)/(plZ+miZ);
      real_type sechZ = 2.0/(plZ+miZ);
      real_type sechZ2 = sechZ*sechZ;
      real_type ret = prt - c0 + c1*tanhZ;
      prt = prt/r;
      dudr = -5*prt + c3*sechZ2;
      prt = prt/r;
      d2udr2 = 30*prt + c5*sechZ2*tanhZ;
      return ret;
    }
    else
    {
      dudr = 0.0;
      d2udr2= 0.0;
      return 0.0;
    };
  }

  /** return a value at r
   */
  real_type f(real_type r)
  {
    return evaluate(r);
  }

  /** return a derivative at r
   */
  real_type df(real_type r)
  {
    if (r<RC)
    {
      real_type r1 = (r-RC);
      real_type plZ = std::exp(B*r1);
      real_type miZ = 1.0/plZ;
      real_type sechZ = 2.0/(plZ+miZ);
      real_type sechZ2 = sechZ*sechZ;
      return c2*std::pow(r,-6) + c3*sechZ2;
    }
    else
    {
      return 0.0;
    };
  }

  /** Read in the parameter from the xml input file.
  * @param cur current xmlNode from which the data members are reset
   */
  bool put(xmlNodePtr cur)
  {
    RC = cutoff_radius;
    xmlNodePtr tcur = cur->xmlChildrenNode;
    while(tcur != NULL)
    {
      //@todo Var -> <param(eter) role="opt"/>
      std::string cname((const char*)(tcur->name));
      if(cname == "parameter" || cname == "Var")
      {
        std::string aname((const char*)(xmlGetProp(tcur,(const xmlChar *)"name")));
        std::string idname((const char*)(xmlGetProp(tcur,(const xmlChar *)"id")));
        if(aname == "A")
        {
          ID_A = idname;
          putContent(A,tcur);
        }
        else
          if(aname == "B")
          {
            ID_B = idname;
            putContent(B,tcur);
          }
          else
            if(aname == "RC")
            {
              ID_RC = idname;
              putContent(RC,tcur);
            }
      }
      tcur = tcur->next;
    }
    myVars.insert(ID_A,A);
    myVars.insert(ID_B,B);
    myVars.insert(ID_RC,RC);
    reset();
    LOGMSG("  Modified McMillan Parameters ")
    LOGMSG("    A (" << ID_A << ") = " << A  << "  B (" << ID_B << ") =  " << B<< "  R_c (" << ID_RC << ") =  " << RC)
    return true;
  }

};

template<class T>
struct comboMcMillanJ2Functor: public OptimizableFunctorBase
{

  ///f(r) = A*r^-5 - A*r_c^-5 + A/(B*r_c^6) * tanh(B*(r-r_c))
  ///made to be used with another jastrow within a cutoff distance R0 due to divergent cusp.
  typedef typename OptimizableFunctorBase::real_type real_type;

  ///coefficients
  real_type A;
  real_type B;
  real_type RC;
  real_type R0;
  real_type cA,c0,c1,c2,c3,c4,c5;
  real_type b0,b1,b2,b3;

  std::string ID_A,ID_B,ID_RC,ID_R0;

  /** constructor
  * @param a A coefficient
  * @param samespin boolean to indicate if this function is for parallel spins
   */
  comboMcMillanJ2Functor(real_type a=0.0, real_type b=0.0, real_type r0=0.0 , real_type rc=0.0):ID_RC("Rcutoff"),ID_R0("Rinner"),ID_A("A"),ID_B("B")
  {
  }


  void checkInVariables(opt_variables_type& active)
  {
    active.insertFrom(myVars);
  }

  void checkOutVariables(const opt_variables_type& active)
  {
    myVars.getIndex(active);
  }

  OptimizableFunctorBase* makeClone() const
  {
    return new comboMcMillanJ2Functor<T>(*this);
  }

  void resetParameters(const opt_variables_type& active)
  {
    int ia=myVars.where(0);
    if(ia>-1)
      A=active[ia];
    int ib=myVars.where(1);
    if(ib>-1)
      B=active[ib];
    int ic=myVars.where(2);
    if(ic>-1)
      RC=active[ic];
    int id=myVars.where(3);
    if(id>-1)
      R0=active[id];
    reset();
  }

  inline void reset()
  {
    cA = std::pow(A,5);
    c0 = cA*std::pow(RC,-5);
    c1 = 5*cA*std::pow(RC,-6)/B;
    c2 = -5*cA;
    c3 = B*c1;
    c4 = 30*cA;
    c5 = -2*B*c3;
    real_type r1 = (R0-RC);
    real_type plZ = std::exp(B*r1);
    real_type miZ = 1.0/plZ;
    real_type tanhZ = (plZ-miZ)/(plZ+miZ);
    real_type sechZ = 2.0/(plZ+miZ);
    real_type sechZ2 = sechZ*sechZ;
    real_type btemp = cA/std::pow(R0,-7);
    b3 = 30.0*btemp + c5*sechZ2*tanhZ;
    b2 = 0.5*b3;
    b1 = -b3*R0 - 5.0*btemp*R0 + c3*sechZ2;
    b0 = -b2*R0*R0 -b1*R0 + btemp*R0*R0 - c0 +c1*tanhZ;
  }


  /** evaluate the value at r
  * @param r the distance
  * @return \f$ u(r) = \frac{A}{r}\left[1-\exp(-\frac{r}{F})\right]\f$
   */
  inline real_type evaluate(real_type r)
  {
    if ((r<RC)&&(r>R0))
    {
      real_type r1 = (r-RC);
      real_type plZ = std::exp(B*r1);
      real_type miZ = 1.0/plZ;
      real_type tanhZ = (plZ-miZ)/(plZ+miZ);
      return cA*std::pow(r,-5)-c0 + c1*tanhZ;
    }
    else
      if ((r<RC)&&(r<=R0))
      {
//         real_type r1 = (r-RC);
//         real_type plZ = std::exp(B*r1);
//         real_type miZ = 1.0/plZ;
//         real_type tanhZ = (plZ-miZ)/(plZ+miZ);
//         real_type sechZ = 2.0/(plZ+miZ);
//         real_type sechZ2 = sechZ*sechZ;
//         real_type ret = b0 + b1*r + b2*r*r - c0 + c1*tanhZ;
        real_type ret = b0 + b1*r + b2*r*r;
        return ret;
      }
      else
      {
        return 0.0;
      };
  }

  /**@param r the distance
  @param dudr first derivative
  @param d2udr2 second derivative
  @return the value
   */
  inline real_type evaluate(real_type r, real_type& dudr, real_type& d2udr2)
  {
    if ((r<RC)&&(r>R0))
    {
      real_type prt = cA*std::pow(r,-5);
      real_type r1 = (r-RC);
      real_type plZ = std::exp(B*r1);
      real_type miZ = 1.0/plZ;
      real_type tanhZ = (plZ-miZ)/(plZ+miZ);
      real_type sechZ = 2.0/(plZ+miZ);
      real_type sechZ2 = sechZ*sechZ;
      real_type ret = prt - c0 + c1*tanhZ;
      prt = prt/r;
      dudr = -5*prt + c3*sechZ2;
      prt = prt/r;
      d2udr2 = 30*prt + c5*sechZ2*tanhZ;
      return ret;
    }
    else
      if ((r<RC)&&(r<=R0))
      {
//         real_type r1 = (r-RC);
//         real_type plZ = std::exp(B*r1);
//         real_type miZ = 1.0/plZ;
//         real_type tanhZ = (plZ-miZ)/(plZ+miZ);
//         real_type sechZ = 2.0/(plZ+miZ);
//         real_type sechZ2 = sechZ*sechZ;
//         real_type ret = b0+b1*r+b2*r*r - c0 + c1*tanhZ;
//         dudr =  b1 + b3*r + c3*sechZ2;
//         d2udr2 =  b3 + c5*sechZ2*tanhZ;
        real_type ret = b0 + b1*r + b2*r*r ;
        dudr =  b1 + b3*r ;
        d2udr2 =  b3 ;
        return ret;
      }
      else
      {
        dudr = 0.0;
        d2udr2= 0.0;
        return 0.0;
      };
  }

  /** return a value at r
   */
  real_type f(real_type r)
  {
    return evaluate(r);
  }

  /** return a derivative at r
   */
  real_type df(real_type r)
  {
    if ((r<RC)&&(r>R0))
    {
      real_type r1 = (r-RC);
      real_type plZ = std::exp(B*r1);
      real_type miZ = 1.0/plZ;
      real_type sechZ = 2.0/(plZ+miZ);
      real_type sechZ2 = sechZ*sechZ;
      return c2*std::pow(r,-6) + c3*sechZ2;
    }
    else
      if ((r<RC)&&(r<R0))
      {
//         real_type r1 = (r-RC);
//         real_type plZ = std::exp(B*r1);
//         real_type miZ = 1.0/plZ;
//         real_type tanhZ = (plZ-miZ)/(plZ+miZ);
//         real_type sechZ = 2.0/(plZ+miZ);
//         real_type sechZ2 = sechZ*sechZ;
//         real_type ret = b1 + b3*r + c3*sechZ2;
        real_type ret = b1 + b3*r ;
        return ret;
      }
      else
      {
        return 0.0;
      };
  }

  /** Read in the parameter from the xml input file.
  * @param cur current xmlNode from which the data members are reset
   */
  bool put(xmlNodePtr cur)
  {
    RC = cutoff_radius;
    xmlNodePtr tcur = cur->xmlChildrenNode;
    while(tcur != NULL)
    {
      //@todo Var -> <param(eter) role="opt"/>
      std::string cname((const char*)(tcur->name));
      if(cname == "parameter" || cname == "Var")
      {
        std::string aname((const char*)(xmlGetProp(tcur,(const xmlChar *)"name")));
        std::string idname((const char*)(xmlGetProp(tcur,(const xmlChar *)"id")));
        if(aname == "A")
        {
          ID_A = idname;
          putContent(A,tcur);
        }
        else
          if(aname == "B")
          {
            ID_B = idname;
            putContent(B,tcur);
          }
          else
            if(aname == "RC")
            {
              ID_RC = idname;
              putContent(RC,tcur);
              if (RC==0.0)
                RC=cutoff_radius;
            }
            else
              if(aname == "R0")
              {
                ID_R0 = idname;
                putContent(R0,tcur);
              }
      }
      tcur = tcur->next;
    }
    myVars.insert(ID_A,A);
    myVars.insert(ID_B,B);
    myVars.insert(ID_RC,RC);
    myVars.insert(ID_R0,R0);
    reset();
    LOGMSG("  Combo McMillan Parameters ")
    LOGMSG("    A (" << ID_A << ") = " << A  << "  B (" << ID_B << ") =  " << B<< "  R_c (" << ID_RC << ") =  "<< RC <<"  R_0 (" << ID_R0 << ") =  " << R0)
    return true;
  }

};
}
#endif

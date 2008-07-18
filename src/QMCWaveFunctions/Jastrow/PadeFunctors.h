//////////////////////////////////////////////////////////////////
// (c) Copyright 2003  by Jeongnim Kim
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
/** @file PadeFunctors.h
 * @brief Functors which implement Pade functions
 */
#ifndef QMCPLUSPLUS_PADEFUNCTORS_H
#define QMCPLUSPLUS_PADEFUNCTORS_H
#include "Numerics/OptimizableFunctorBase.h"
#include <cmath>

namespace qmcplusplus {

  /** Implements a Pade Function \f$u[r]=A*r/(1+B*r)\f$
   * 
   * Similar to PadeJastrow with a scale.
   */
  template<class T>
    struct PadeFunctor:public OptimizableFunctorBase<T> {

      typedef typename OptimizableFunctorBase<T>::real_type real_type;
      typedef typename OptimizableFunctorBase<T>::OptimizableSetType OptimizableSetType;

      ///input A
      real_type A;
      ///input B
      real_type B0;
      ///input scaling, default=1.0
      real_type Scale;
      ///B=B0*Scale
      real_type B; 
      ///AB=A*B
      real_type AB;
      ///B2=2*B
      real_type B2;
      ///AoverB=A/B 
      real_type AoverB;
      ///id of A
      std::string ID_A;
      ///id of B
      std::string ID_B;

      ///default constructor
      PadeFunctor(): Scale(1.0),ID_A("0"),ID_B("0") { }

      OptimizableFunctorBase<T>* makeClone() const 
      {
        return new PadeFunctor<T>(*this);
      }

      ///constructor
      explicit PadeFunctor(real_type a, real_type b, real_type s=1.0): 
        Scale(s),ID_A("0"), ID_B("0") 
        {
          reset(a,b);
        }

      void reset(real_type a, real_type b, real_type s=1.0) {
        A=a; B0=b; Scale=s;
        B = B0*Scale;
        AB = A*B; B2=2.0*B;
        AoverB=A/B;
      }

      inline real_type evaluate(real_type r) {
        return A*r/(1.0+B*r);
      }

      inline real_type 
        evaluate(real_type r, real_type& dudr, real_type& d2udr2) {
          real_type u = 1.0/(1.0+B*r);
          dudr = A*u*u;
          d2udr2 = -B2*dudr*u;
          return A*u*r;
        }

      inline real_type f(real_type r) {
        return evaluate(r)-AoverB;
      }

      inline real_type df(real_type r) {
        real_type dudr,d2udr2;
        real_type res=evaluate(r,dudr,d2udr2);
        return dudr;
      }

      bool put(xmlNodePtr cur) {
        real_type Atemp(A),Btemp(B);
        cur = cur->xmlChildrenNode;
        while(cur != NULL) {
          //@todo Var -> <param(eter) role="opt"/>
          std::string cname((const char*)(cur->name));
          if(cname == "parameter" || cname == "Var") {
            const xmlChar* nptr=xmlGetProp(cur,(const xmlChar *)"name");
            const xmlChar* iptr=xmlGetProp(cur,(const xmlChar *)"id");
            if(nptr == NULL || iptr == NULL) return false;
            if(nptr[0] == 'A')
            { 
              ID_A = (const char*)iptr; 
              putContent(Atemp,cur);
            } else if(nptr[0] == 'B'){
              ID_B = (const char*)iptr;
              putContent(Btemp,cur);
            }
          }
          cur = cur->next;
        }
        reset(Atemp,Btemp);
        return true;
      }

      void addOptimizables(OptimizableSetType& vlist)
      {
        if(ID_A != "0") vlist[ID_A]=A;
        if(ID_B != "0") vlist[ID_B]=B0;
      }

      /** reset the internal variables.
       *
       * USE_resetParameters
       */
      void resetParameters(OptimizableSetType& optVariables) 
      {
        typename OptimizableSetType::iterator it_a(optVariables.find(ID_A));
        if(it_a != optVariables.end()) A=(*it_a).second;
        typename OptimizableSetType::iterator it_b(optVariables.find(ID_B));
        if(it_b != optVariables.end()) B0=(*it_b).second;
        B = B0*Scale;
        AB = A*B; 
        B2=2.0*B;
        AoverB=A/B;
      }
    };


  /** Pade function of \f[ u(r) = \frac{a*r+c*r^2}{1+b*r} \f]
   *
   * Prototype of the template parameter of TwoBodyJastrow and OneBodyJastrow
   */
  template<class T>
    struct Pade2ndOrderFunctor:public OptimizableFunctorBase<T> {

      typedef typename OptimizableFunctorBase<T>::real_type real_type;
      typedef typename OptimizableFunctorBase<T>::OptimizableSetType OptimizableSetType;
      ///coefficients
      real_type A, B, C, C2;
      ///id for A
      std::string ID_A;
      ///id for B
      std::string ID_B;
      ///id for C
      std::string ID_C;

      ///constructor
      Pade2ndOrderFunctor(real_type a=1.0, real_type b=1.0, real_type c=1.0) {reset(a,b,c);}

      //Pade2ndOrderFunctor(Pade2ndOrderFunctor<T>* func) {
      //  reset(1.0,1.0,1.0);
      //}
      OptimizableFunctorBase<T>* makeClone() const 
      {
        return new Pade2ndOrderFunctor<T>(*this);
      }

      /** reset the internal variables.
       *@param a Pade Jastrow parameter a 
       *@param b Pade Jastrow parameter b 
       *@param c Pade Jastrow parameter c 
       */
      void reset(real_type a, real_type b, real_type c) {
        A = a; B=b; C = c; C2 = 2.0*C;
      }

      /**@param r the distance
        @return \f$ u(r) = a*r/(1+b*r) \f$
        */
      inline real_type evaluate(real_type r) {
        real_type br(B*r);
        return (A+br)*r/(1.0+br);
      }

      /** evaluate the value at r
       * @param r the distance
       @param dudr return value  \f$ du/dr = a/(1+br)^2 \f$
       @param d2udr2 return value  \f$ d^2u/dr^2 = -2ab/(1+br)^3 \f$
       @return \f$ u(r) = a*r/(1+b*r) \f$
       */
      inline real_type evaluate(real_type r, real_type& dudr, real_type& d2udr2) {
        real_type u = 1.0/(1.0+B*r);
        real_type v = A*r+B*r*r;
        real_type w = A+C2*r;
        dudr = u*(w-B*u*v);
        d2udr2 = 2.0*u*(C-B*dudr);
        return u*v;
      }

      real_type f(real_type r) {
        return evaluate(r);
      }

      real_type df(real_type r) {
        real_type dudr,d2udr2;
        real_type res=evaluate(r,dudr,d2udr2);
        return dudr;
      }

      /** process input xml node
       * @param cur current xmlNode from which the data members are reset
       *
       * T1 is the type of VarRegistry, typically double.  
       * Read in the Pade parameters from the xml input file.
       */
      bool put(xmlNodePtr cur){
        real_type Atemp,Btemp, Ctemp;
        //jastrow[iab]->put(cur->xmlChildrenNode,wfs_ref.RealVars);
        xmlNodePtr tcur = cur->xmlChildrenNode;
        while(tcur != NULL) {
          //@todo Var -> <param(eter) role="opt"/>
          std::string cname((const char*)(tcur->name));
          if(cname == "parameter") {
            std::string aname((const char*)(xmlGetProp(tcur,(const xmlChar *)"name")));
            std::string idname((const char*)(xmlGetProp(tcur,(const xmlChar *)"id")));
            if(aname == "A") {
              ID_A = idname;
              putContent(Atemp,tcur);
            } else if(aname == "B"){
              ID_B = idname;
              putContent(Btemp,tcur);
            } else if(aname == "C") {
              ID_C = idname;
              putContent(Ctemp,tcur);
            }
          }
          tcur = tcur->next;
        }
        reset(Atemp,Btemp,Ctemp);
        //LOGMSG("Jastrow (A*r+C*r*r)/(1+Br) = (" << A << "," << B << "," << C << ")") 
          return true;
      }

      /** add optimizable variables to vlist
       * @param vlist VarRegistry<T1> to which the Pade variables A and B are added for optimization
       */
      void addOptimizables(OptimizableSetType& vlist)
      {
        vlist[ID_A]=A;
        vlist[ID_B]=B;
        vlist[ID_C]=C;
      }

      void resetParameters(OptimizableSetType& optVariables) 
      {
        typename OptimizableSetType::iterator it(optVariables.find(ID_A));
        if(it != optVariables.end()) A=(*it).second;
        it=optVariables.find(ID_B);
        if(it != optVariables.end()) B=(*it).second;
        it=optVariables.find(ID_C);
        if(it != optVariables.end()) C=(*it).second;
        C2 = 2.0*C;
      }
    };

  /** Pade functional of \f[ u(r) = \frac{a*f(r)}{1+b*f(r)} \f] with a scale function f(r)
   *
   * Prototype of the template parameter of TwoBodyJastrow and OneBodyJastrow
   */
  template<class T>
    struct ScaledPadeFunctor:public OptimizableFunctorBase<T> {

      typedef typename OptimizableFunctorBase<T>::real_type real_type;
      typedef typename OptimizableFunctorBase<T>::OptimizableSetType OptimizableSetType;

      ///coefficients
      real_type A, B, C; 
      real_type OneOverC, B2;

      ///constructor
      explicit ScaledPadeFunctor(real_type a=1.0, real_type b=1.0, real_type c=1.0) 
      {reset(a,b,c);}

      OptimizableFunctorBase<T>* makeClone() const 
      {
        return new ScaledPadeFunctor<T>(*this);
      }

      /** reset the internal variables.
       *@param a Pade Jastrow parameter a 
       *@param b Pade Jastrow parameter b 
       */
      void reset(real_type a, real_type b, real_type c) {
        A=a; B=b; C=c; 
        OneOverC=1.0/c;
        B2=2.0*b;
      }

      /** evaluate the value at r
       * @param r the distance
       * @return \f$ u(r_{eff}) = a*r_{eff}/(1+b*r_{eff}) \f$
       */
      inline real_type evaluate(real_type r) {
        real_type reff((1.0-std::exp(-C*r))*OneOverC);
        return A*reff/(1.0+B*reff);
      }

      /** evaluate the value, first derivative and second derivative
       * @param r the distance
       * @param dudr return value  \f$ du/dr\f$
       * @param d2udr2 return value  \f$ d^2u/dr^2 \f$
       * @return \f$ u(r_{eff}) = a*r_{eff}/(1+b*r_{eff}) \f$
       */
      inline real_type 
        evaluate(real_type r, real_type& dudr, real_type& d2udr2) {
          real_type reff_1(std::exp(-C*r));
          real_type reff((1.0-reff_1)*OneOverC);
          real_type u(1.0/(1.0+B*reff));
          real_type auu(A*u*u);
          dudr=reff_1*auu;
          //d2udr=auu*(-C*reff_1*reff_2-B2*reff_1*reff_1*u);
          d2udr2=-reff_1*auu*(C+B2*reff_1*u);
          return A*u*reff;
        }

      real_type f(real_type r) {
        return evaluate(r);
      }

      real_type df(real_type r) {
        real_type dudr,d2udr2;
        real_type res=evaluate(r,dudr,d2udr2);
        return dudr;
      }

      bool put(xmlNodePtr cur) {return true;}
      void addOptimizables(OptimizableSetType& vlist){}

      /** reset the internal variables.
       *
       * When RefPade is not 0, use RefPade->B to reset the values
       */
      inline void resetParameters(OptimizableSetType& optVariables) 
      {
        OneOverC=1.0/C;
        B2=2.0*B;
      }
    };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/


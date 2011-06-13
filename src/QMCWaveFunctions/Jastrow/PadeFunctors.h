//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
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
#include "OhmmsData/AttributeSet.h"
#include <cmath>
// #include <vector>
#include "OhmmsPETE/TinyVector.h"


namespace qmcplusplus {

  /** Implements a Pade Function \f$u[r]=A*r/(1+B*r)\f$
   * 
   * Similar to PadeJastrow with a scale.
   */
  template<class T>
    struct PadeFunctor:public OptimizableFunctorBase {
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
      ///constructor
      explicit PadeFunctor(real_type a, real_type b, real_type s=1.0): 
        A(a),B0(b),Scale(s)
      {
        reset();
      }


      /** set ID_A and ID_B
       * @param id_a ID of A
       * @param id_b ID of B
       * @param free_a if true, A is optimizable
       * @param free_b if true, B is optimizable
       */
      inline void setIDs(const std::string& id_a, const std::string& id_b, bool free_a=false, bool free_b=true)
      {
        ID_A=id_a;
        ID_B=id_b;
        myVars.clear();
        myVars.insert(ID_A,A,free_a);
        myVars.insert(ID_B,B0,free_b);
      }

      OptimizableFunctorBase* makeClone() const
      {
        return new PadeFunctor(*this);
      }

      void reset() {
        //A=a; B0=b; Scale=s;
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

      inline real_type 
      evaluate(real_type r, real_type& dudr, real_type& d2udr2, real_type d3udr3) {
          real_type u = 1.0/(1.0+B*r);
          dudr = A*u*u;
          d2udr2 = -B2*dudr*u;
	  std::cerr << "Third derivative not imlemented for Pade functor.\n";
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
      
      inline bool evaluateDerivatives (real_type r, std::vector<TinyVector<real_type,3> >& derivs)
    {
      real_type u = 1.0/(1.0+B*r);
      derivs[0][0]= r*u;
      derivs[1][0]= -A*r*r*u*u;
      derivs[0][1]= u*u;
      derivs[1][1]= -2.0*A*r*u*u*u;
      derivs[0][2]= -B2*u*u*u;
      derivs[1][2]= 2*A*(B2*r-1)*u*u*u*u;
      return true; 
    }

      bool put(xmlNodePtr cur) {
        real_type Atemp(A),Btemp(B0);
        cur = cur->xmlChildrenNode;
        bool renewed=false;
        while(cur != NULL) 
        {
          //@todo Var -> <param(eter) role="opt"/>
          std::string cname((const char*)(cur->name));
          if(cname == "parameter" || cname == "Var") 
          {
            std::string id_in("0");
            std::string p_name("B");
            OhmmsAttributeSet rAttrib;
            rAttrib.add(id_in, "id");
            rAttrib.add(p_name, "name"); 
            rAttrib.put(cur);
            if(p_name=="A")
            {
              ID_A = id_in;
              putContent(Atemp,cur);
              renewed=true;
            } else if(p_name == "B"){
              ID_B = id_in;
              putContent(Btemp,cur);
              renewed=true;
            }
          }
          cur = cur->next;
        }
        if(renewed)
        {
          A=Atemp;
          B0=Btemp;
          reset();
          myVars.clear();
          myVars.insert(ID_A,A, ID_A != "0");
          myVars.insert(ID_B,B0,ID_B != "0");
        }
        return true;
      }

      void checkInVariables(opt_variables_type& active)
      {
        active.insertFrom(myVars);
        //std::cout << "Checking variables by PadeFunctor" << endl;
        //myVars.print(std::cout);
      }

      void checkOutVariables(const opt_variables_type& active)
      {
        myVars.getIndex(active);
        //myVars.print(std::cout);
      }

      void resetParameters(const opt_variables_type& active) 
      {
        int ia=myVars.where(0); if(ia>-1) A=myVars[0]=active[ia];
        int ib=myVars.where(1); if(ib>-1) B0=myVars[1]=active[ib];
        reset();
        //B = B0*Scale;
        //AB = A*B; 
        //B2=2.0*B;
        //AoverB=A/B;
      }
    };


  /** Pade function of \f[ u(r) = \frac{a*r+c*r^2}{1+b*r} \f]
   *
   * Prototype of the template parameter of TwoBodyJastrow and OneBodyJastrow
   */
  template<class T>
    struct Pade2ndOrderFunctor:public OptimizableFunctorBase {

      ///coefficients
      real_type A, B, C, C2;
      ///id for A
      std::string ID_A;
      ///id for B
      std::string ID_B;
      ///id for C
      std::string ID_C;

      ///constructor
      Pade2ndOrderFunctor(real_type a=1.0, real_type b=1.0, real_type c=1.0): A(a),B(b),C(c)
      {
        reset();
      }

      OptimizableFunctorBase* makeClone() const
      {
        return new Pade2ndOrderFunctor(*this);
      }

      /** reset the internal variables.
       */
      void reset() {
       // A = a; B=b; C = c; 
        C2 = 2.0*C;
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
        real_type v = A*r+C*r*r;
        real_type w = A+C2*r;
        dudr = u*(w-B*u*v);
        d2udr2 = 2.0*u*u*u*(C-B*A);
        return u*v;
      }

    inline real_type evaluate(real_type r, real_type& dudr, 
			      real_type& d2udr2, real_type& d3udr3) {
        real_type u = 1.0/(1.0+B*r);
        real_type v = A*r+C*r*r;
        real_type w = A+C2*r;
        dudr = u*(w-B*u*v);
        d2udr2 = 2.0*u*u*u*(C-B*A);
	std::cerr << "Third derivative not imlemented for Pade functor.\n";
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
    
    inline bool evaluateDerivatives (real_type r, std::vector<TinyVector<real_type,3> >& derivs)
    {
      real_type u = 1.0/(1.0+B*r);
      real_type u2 = u*u;
      real_type u3 = u*u*u;
      real_type u4 = u*u*u*u;
      
      derivs[0][0]= r*u;
      derivs[1][0]= -r*r*(A+C*r)*u2;
      derivs[2][0]= r*r*u;
      
      derivs[0][1]= u2;
      derivs[1][1]= -r*(2.0*A+C*r*(3.0+B*r))*u3;
      derivs[2][1]= r*(2.0+B*r)*u2;
      
      derivs[0][2]= -2.0*B*u3;
      derivs[1][2]= -2.0*u4*(A-2.0*A*B*r+3.0*C*r);
      derivs[2][2]= 2.0*u3;
      return true; 
    }
      
      
      

      /** process input xml node
       * @param cur current xmlNode from which the data members are reset
       *
       * T1 is the type of VarRegistry, typically double.  
       * Read in the Pade parameters from the xml input file.
       */
      bool put(xmlNodePtr cur){
        real_type Atemp,Btemp, Ctemp;
        ID_A="pade2A";
        ID_B="pade2B";
        ID_C="pade2C";
        //jastrow[iab]->put(cur->xmlChildrenNode,wfs_ref.RealVars);
        xmlNodePtr tcur = cur->xmlChildrenNode;
        bool renewed=false;
        while(tcur != NULL) {
          std::string cname((const char*)(tcur->name));
          if(cname == "parameter" || cname == "Var") 
          {
            std::string id_in("0");
            std::string p_name("B");
            OhmmsAttributeSet rAttrib;
            rAttrib.add(id_in, "id");
            rAttrib.add(p_name, "name"); 
            rAttrib.put(tcur);
            if(p_name=="A")
            {
              if (id_in!="0") ID_A = id_in;
              putContent(Atemp,tcur);
              renewed=true;
            } else if(p_name == "B"){
              if (id_in!="0") ID_B = id_in;
              putContent(Btemp,tcur);
              renewed=true;
            } else if(p_name == "C"){
              if (id_in!="0") ID_C = id_in;
              putContent(Ctemp,tcur);
              renewed=true;
            }
          }

          tcur = tcur->next;
        }
        if (renewed)
        {
        A=Atemp; B=Btemp; C=Ctemp;
        reset();
        //these are always active 
        myVars.clear();
        myVars.insert(ID_A,A,ID_A!="pade2A");
        myVars.insert(ID_B,B,ID_B!="pade2B");
        myVars.insert(ID_C,C,ID_C!="pade2C");
        }
        //LOGMSG("Jastrow (A*r+C*r*r)/(1+Br) = (" << A << "," << B << "," << C << ")") 
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
        if (ID_A!="0") {int ia=myVars.where(0); if(ia>-1) A=myVars[0]=active[ia];}
        if (ID_B!="0") {int ib=myVars.where(1); if(ib>-1) B=myVars[1]=active[ib];}
        if (ID_C!="0") {int ic=myVars.where(2); if(ic>-1) C=myVars[2]=active[ic];}
        C2 = 2.0*C;
      }
    };

     /** Pade function of \f[ u(r) = \frac{a*r+b*r^2}{1+c*r+d*r^2} \f]
   *
   * Prototype of the template parameter of TwoBodyJastrow and OneBodyJastrow
   */
  template<class T>
    struct PadeTwo2ndOrderFunctor:public OptimizableFunctorBase {

      ///coefficients
      real_type A, B, C, D;
      ///id for A
      std::string ID_A;
      ///id for B
      std::string ID_B;
      ///id for C
      std::string ID_C;
      ///id for D
      std::string ID_D;
      
      ///constructor
      PadeTwo2ndOrderFunctor(real_type a=1.0, real_type b=1.0, real_type c=1.0, real_type d=1.0): A(a),B(b),C(c),D(d)
      {
        reset();
      }

      OptimizableFunctorBase* makeClone() const
      {
        return new PadeTwo2ndOrderFunctor(*this);
      }

      /** reset the internal variables.
       */
      void reset() {
       // A = a; B=b; C = c;  
      }

      inline real_type evaluate(real_type r) {
        real_type br(B*r);
        real_type dr(D*r);
        return (A*r+br*r)/(1.0+C*C*r+dr*dr);
      }

      /** evaluate the value at r
       * @param r the distance
       @param dudr return value 
       @param d2udr2 return value 
       */
      inline real_type evaluate(real_type r, real_type& dudr, real_type& d2udr2) 
      {
        real_type ar(A*r);
        real_type br(B*r);
        real_type cr(C*r);
        real_type dr(D*r);
        real_type bttm( 1.0/(1.0+C*cr+dr*dr) );
        dudr = (A - A*dr*dr + br*(2.0 + C*cr))*bttm*bttm;
        d2udr2 = -2.0*(A*(C*C + 3.0*dr*D - dr*dr*dr*D) + B*(-1.0 + dr*dr*(3.0 + C*cr)))*bttm*bttm*bttm;
        return (A*r+br*r)*bttm;
      }

    inline real_type evaluate(real_type r, real_type& dudr, 
            real_type& d2udr2, real_type& d3udr3) {
        d3udr3 = 0;
        return evaluate(r,dudr,d2udr2);
      }


      real_type f(real_type r) {
        return evaluate(r);
      }

      real_type df(real_type r) {
        real_type dudr,d2udr2;
        real_type res=evaluate(r,dudr,d2udr2);
        return dudr;
      }
    
    inline bool evaluateDerivatives (real_type r, std::vector<TinyVector<real_type,3> >& derivs)
    {
      real_type ar(A*r);
      real_type br(B*r);
      real_type cr(C*r);
      real_type dr(D*r);
      real_type r2(r*r);
      real_type dr2(D*r2);
      real_type d2r2(dr*dr);
      real_type c2(C*C);
      real_type c2r(c2*r);
      real_type bttm(1.0/(1.0+c2r+d2r2));
      real_type tp( A*r+br*r );
      
      real_type c2r2(cr*cr);
      real_type d4r4(d2r2*d2r2);
      real_type bttm2(bttm*bttm);
      real_type bttm3(bttm*bttm2);
      real_type bttm4(bttm2*bttm2);
      
      derivs[0][0]= r2*bttm;
      derivs[1][0]= -2.0*cr*tp*bttm2;
      derivs[2][0]= -2.0*dr2*tp*bttm2;

      derivs[0][1]= r*(2.0+c2r)*bttm2;
      derivs[1][1]= -2*cr*(A*(2.0 - 2.0*d2r2) + br*(3.0 + c2r - d2r2))*bttm3;
      derivs[2][1]= -2.0*dr2*(2.0*br*(2.0 + c2r) + A*(3.0 + c2r - d2r2))*bttm3;
      
      derivs[0][2]= (2.0 - 2.0*d2r2*(3.0+c2r))*bttm3;
      derivs[1][2]= 4.0*C*(A*(-1.0 + 2.0*c2r + 8.0*d2r2 - 3.0*d4r4) + br*(-3.0 - d4r4 + 2.0*d2r2 *(4.0 + c2r)))*bttm4;
      derivs[2][2]= -4.0*dr*(br*(6.0 + 4.0* c2r + c2*c2r2 - 6.0*d2r2 - 2.0*c2r*d2r2) + A*(3.0 + d4r4 - 2.0 *d2r2 *(4.0 + c2r))) *bttm4;
      return true; 
    }
      
      
      

      /** process input xml node
       * @param cur current xmlNode from which the data members are reset
       *
       * T1 is the type of VarRegistry, typically double.  
       * Read in the Pade parameters from the xml input file.
       */
      bool put(xmlNodePtr cur){
        std::string fcup("yes");
        OhmmsAttributeSet p;
        p.add(fcup,"fixcusp");
        p.put(cur);
        if (fcup=="true") fcup="yes";
        // if (fcup=="yes") app_log()<<" fixing cusp conditions"<<endl;


        real_type Atemp,Btemp, Ctemp, Dtemp;
        ID_A="pade2A";
        ID_B="pade2B";
        ID_C="pade2C";
        ID_D="pade2D";
        //jastrow[iab]->put(cur->xmlChildrenNode,wfs_ref.RealVars);
        xmlNodePtr tcur = cur->xmlChildrenNode;
        bool renewed=false;
        while(tcur != NULL) {
          std::string cname((const char*)(tcur->name));
          if(cname == "parameter" || cname == "Var") 
          {
            std::string id_in("0");
            std::string p_name("B");
            OhmmsAttributeSet rAttrib;
            rAttrib.add(id_in, "id");
            rAttrib.add(p_name, "name"); 
            rAttrib.put(tcur);
            if(p_name=="A")
            {
              if (id_in!="0") ID_A = id_in;
              putContent(Atemp,tcur);
              renewed=true;
            } else if(p_name == "B"){
              if (id_in!="0") ID_B = id_in;
              putContent(Btemp,tcur);
              renewed=true;
            } else if(p_name == "C"){
              if (id_in!="0") ID_C = id_in;
              putContent(Ctemp,tcur);
              renewed=true;
            } else if(p_name == "D"){
              if (id_in!="0") ID_D = id_in;
              putContent(Dtemp,tcur);
              renewed=true;
            }
          }

          tcur = tcur->next;
        }
        if (renewed)
        {
        A=Atemp; B=Btemp; C=Ctemp;  D=Dtemp;
        reset();
        //these are always active 
        myVars.clear();
        
        myVars.insert(ID_B,B,ID_B!="pade2B");
        myVars.insert(ID_C,C,ID_C!="pade2C");
        myVars.insert(ID_D,D,ID_D!="pade2D");
        //myVars.insert(ID_A,A,fcup!="yes");
        }
        //LOGMSG("Jastrow (A*r+C*r*r)/(1+Br) = (" << A << "," << B << "," << C << ")") 
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
        int ib=myVars.where(0); if(ib>-1) B=myVars[0]=active[ib];
        int ic=myVars.where(1); if(ic>-1) C=myVars[1]=active[ic];
        int id=myVars.where(2); if(id>-1) D=myVars[2]=active[id]; 
        //int ia=myVars.where(3); if(ia>-1) A=myVars[3]=active[ia];
      }
    };
    
  /** Pade functional of \f[ u(r) = \frac{a*f(r)}{1+b*f(r)} \f] with a scale function f(r)
   *
   * Prototype of the template parameter of TwoBodyJastrow and OneBodyJastrow
   */
  template<class T>
    struct ScaledPadeFunctor:public OptimizableFunctorBase {

      ///coefficients
      real_type A, B, C; 
      real_type OneOverC, B2;

      ///constructor
      explicit ScaledPadeFunctor(real_type a=1.0, real_type b=1.0, real_type c=1.0) :
        A(a),B(b),C(c)
      {reset();}

      OptimizableFunctorBase* makeClone() const
      {
        return new ScaledPadeFunctor(*this);
      }

      /** reset the internal variables.
       */
      void reset() {
        OneOverC=1.0/C;
        B2=2.0*B;
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

      inline real_type 
        evaluate(real_type r, real_type& dudr, real_type& d2udr2, real_type d3udr3) {
	  real_type reff_1(std::exp(-C*r));
	  real_type reff((1.0-reff_1)*OneOverC);
	  real_type u(1.0/(1.0+B*reff));
	  real_type auu(A*u*u);
	  dudr=reff_1*auu;
	  //d2udr=auu*(-C*reff_1*reff_2-B2*reff_1*reff_1*u);
	  d2udr2=-reff_1*auu*(C+B2*reff_1*u);
	  std::cerr << "Third derivative not imlemented for Pade functor.\n";
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

      void checkInVariables(opt_variables_type& active)
      {
        active.insertFrom(myVars);
      }

      void checkOutVariables(const opt_variables_type& active)
      {
        myVars.getIndex(active);
      }

      inline void resetParameters(const opt_variables_type& active) 
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


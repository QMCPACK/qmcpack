//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Eric Neuscamman, eneuscamman@berkeley.edu, University of California, Berkeley
//
// File created by: Eric Neuscamman, eneuscamman@berkeley.edu, University of California, Berkeley
//////////////////////////////////////////////////////////////////////////////////////


/** @file ShortRangeCuspFunctor.h
 * @brief Functor designed to encode short-ranged structure near a nuclear cusp
 */
#ifndef QMCPLUSPLUS_SHORTRANGECUSPFUNCTOR_H
#define QMCPLUSPLUS_SHORTRANGECUSPFUNCTOR_H
#include "OptimizableFunctorBase.h"
#include "OhmmsData/AttributeSet.h"
#include <cmath>
#include <stdexcept>
#include "OhmmsPETE/TinyVector.h"
#include "Utilities/ProgressReportEngine.h"
#include "ModernStringUtils.hpp"

namespace qmcplusplus
{

/*********************************************************************************************************
*
*  A functor intended for helping with short-ranged structure near electron-nuclear (en) cusps.
*
*  In QMCPACK, an en radial Jastrow is parameterized as exp(-U(r)), in which U(r) is the functor.
*  This class implements the functor
*
*    U(r) = -1.0 * exp(-r/R0) * ( A * R0 + sum_{k=0}^{N-1} B_k * (r/R0)^{k+2} / ( 1.0 + (r/R0)^{k+2} ) )
*
*  in which A, R0, and B_k are the variational parameters, although A is likely to be held fixed
*  as it enforces the en cusp condition ( dU/dr -> A as r -> 0 ).
*
*  The exp(-r/R0) keeps the functor short ranged.  For example, initial testing on the LiH molecule
*  suggests that in that case, a value of about 0.03 Bohr is desirable for the soft cutoff parameter R0.
*  With this R0 value, r values beyond just 0.5 Bohr already incur an exp(-r/R0) factor smaller than
*  10^{-7} and so the functor is quite short ranged.
*
*  The B_k variables add structure through an expansion in a set of sigmoidal functions.
*  The sigmoidal form (as opposed to a bare polynomial expansion) of
*
*                          (r/R0)^{k+2} / ( 1.0 + (r/R0)^{k+2} )
*
*  was chosen so that these do not compete with the exponential decay and so help keep things short ranged.
*  Note that the lowest power of (r/R0) in this expansion is 2, so it does not affect the cusp condition.
*
**********************************************************************************************************/
template<class T>
struct ShortRangeCuspFunctor : public OptimizableFunctorBase
{

  //*****************************************************************************//
  //*******************************  Member Data ********************************//
  //*****************************************************************************//

  ///true, if A is optimizable
  bool Opt_A;
  ///true, if R0 is optimizable
  bool Opt_R0;
  ///true, if B variables are optimizable
  bool Opt_B;
  ///variable that controls the cusp
  real_type A;
  ///the soft cutoff distance that controls how short ranged the functor is
  real_type R0;
  ///variables that add detail through an expansion in sigmoidal functions
  std::vector<real_type> B;
  ///id of A
  std::string ID_A;
  ///id of R0
  std::string ID_R0;
  ///id of B
  std::string ID_B;

  //*****************************************************************************//
  //*****************************  Member Functions *****************************//
  //*****************************************************************************//

  ///default constructor
  ShortRangeCuspFunctor() : Opt_A(false), Opt_R0(true), Opt_B(true), A(1.0), R0(0.06),
                            ID_A("string_not_set"), ID_R0("string_not_set"), ID_B("string_not_set")
  {
    cutoff_radius = 1.0e4; //some big range
    reset();
  }

  ///sets the cusp condition and disables optimization of the cusp-determining parameter
  void setCusp(real_type cusp) override
  {
    //throw std::runtime_error("ShortRangeCuspFunctor::setCusp was called");
    A     = cusp;
    Opt_A = false;
    reset();
  }

  ///clone the functor
  OptimizableFunctorBase* makeClone() const override { return new ShortRangeCuspFunctor(*this); }

  ///Implement the reset function, which was pure virtual in OptimizableFunctorBase, even though we don't need it
  void reset() override
  {
    //cutoff_radius = 1.0e4; //some big range
  }

  ///compute U(r) at a particular value of r
  inline real_type evaluate(real_type r) const
  {

    // the functor will be almost exactly zero at long range, so just return that if r is beyond the hard cutoff
    if (r >= cutoff_radius)
      return 0.0;

    // get the ratio of the distance and the soft cutoff distance
    const real_type s = r / R0;

    // sum up the sigmoidal function expansion
    real_type sig_sum = 0.0;
    real_type sn = s * s; // s^n
    for (int i = 0; i < B.size(); i++) {
      sig_sum += B[i] * sn / ( 1.0 + sn );
      sn *= s;  // update s^n
    }

    // return U(r)
    return -1.0 * std::exp(-s) * ( A * R0 + sig_sum );

  }

  ///compute U(r), dU/dr, and d^2U/dr^2 at a particular value of r
  inline real_type evaluate(real_type r, real_type& dudr, real_type& d2udr2) const
  {

    // the functor will be almost exactly zero at long range, so just return that if r is beyond the hard cutoff
    if (r >= cutoff_radius) {
      dudr = 0.0;
      d2udr2 = 0.0;
      return 0.0;
    }

    // get the ratio of the distance and the soft cutoff distance
    const real_type s = r / R0;

    // get the exponential factor
    const real_type ex = std::exp(-s);

    // sum the terms related to the sigmoidal functions that are needed for U, dU, and d^2U
    real_type sig_sum_0 = 0.0; // contributes to U(r)
    real_type sig_sum_1 = 0.0; // contributes to dU/dr
    real_type sig_sum_2 = 0.0; // contributes to d2U/dr2
    real_type n = 2.0;
    real_type sn = s * s; // s^{n  }
    real_type snm1 = s;   // s^{n-1}
    real_type snm2 = 1.0; // s^{n-2}
    for (int i = 0; i < B.size(); i++) {
      const real_type d = 1.0 + sn;
      const real_type d2 = d * d;
      const real_type soverd = s / d;
      const real_type noverd2 = n / d2;
      sig_sum_0 += B[i] * sn / d;
      sig_sum_1 += B[i] * ( sn * sn + sn - n * snm1 ) / d2;
      sig_sum_2 += B[i] * snm2 * ( d * noverd2 * noverd2 * ( sn - 1.0 ) + noverd2 * ( 1.0 + 2.0 * s ) - soverd * s );
      snm2 = snm1; // update s^{n-2}
      snm1 = sn;   // update s^{n-1}
      sn *= s;     // update s^{n  }
      n += 1.0;    // update n
    }

    // get some useful ratios ( ex / R0 and ex / R0^2 )
    const real_type exoverR0  = ex / R0;
    const real_type exoverR02 = exoverR0 / R0;

    // evaluate dU/dr
    dudr = A * ex + exoverR0 * sig_sum_1;

    // evaluate d^2U/dr^2
    d2udr2 = -A * exoverR0 + exoverR02 * sig_sum_2;

    // return U(r)
    return -1.0 * ex * ( A * R0 + sig_sum_0 );

  }

  ///compute U(r), dU/dr, d^2U/dr^2, and d^3U/dr^3
  inline real_type evaluate(real_type r, real_type& dudr, real_type& d2udr2, real_type& d3udr3) const
  {
    ReportEngine PRE("ShortRangeCuspFunctor", "evaluate(r, dudr, d2udr2, d3udr3)");
    PRE.error("evaluate(r, dudr, d2udr2, d3udr3) not implemented for ShortRangeCuspFunctor", true);
    return 0.0;
  }

  ///compute U(r) at multiple distances and sum the results
  inline real_type evaluateV(const int iat,
                             const int iStart,
                             const int iEnd,
                             const T* restrict _distArray,
                             T* restrict distArrayCompressed) const
  {
    real_type sum(0);
    for (int idx = iStart; idx < iEnd; idx++)
      if (idx != iat)
        sum += evaluate(_distArray[idx]);
    return sum;
  }

  ///compute U(r), dU/dr, and d^2U/dr^2 at multiple values of r
  inline void evaluateVGL(const int iat,
                          const int iStart,
                          const int iEnd,
                          const T* distArray,
                          T* restrict valArray,
                          T* restrict gradArray,
                          T* restrict laplArray,
                          T* restrict distArrayCompressed,
                          int* restrict distIndices) const
  {
    for (int idx = iStart; idx < iEnd; idx++)
    {
      valArray[idx] = evaluate(distArray[idx], gradArray[idx], laplArray[idx]);
      gradArray[idx] /= distArray[idx];
    }
    if (iat >= iStart && iat < iEnd)
      valArray[iat] = gradArray[iat] = laplArray[iat] = T(0);
  }

  ///compute U(r) at a particular distance or return zero if beyond the cutoff
  inline real_type f(real_type r) override
  {
    if (r >= cutoff_radius)
      return 0.0;
    return evaluate(r);
  }

  ///compute dU/dr at a particular distance or return zero if beyond the cutoff
  inline real_type df(real_type r) override
  {
    if (r >= cutoff_radius)
      return 0.0;
    real_type du, d2u;
    evaluate(r, du, d2u);
    return du;
  }

  /// compute derivatives of U(r), dU/dr, and d^2U/dr^2 with respect to the variational parameters
  inline bool evaluateDerivatives(real_type r, std::vector<TinyVector<real_type, 3>>& derivs) override
  {

    // get the ratio of the distance and the soft cutoff distance
    const real_type s = r / R0;

    // get the exponential factor
    const real_type ex = std::exp(-s);

    // get some useful ratios
    const real_type exoverR0  = ex / R0;
    const real_type exoverR02 = exoverR0 / R0;
    const real_type exoverR03 = exoverR02 / R0;

    // initialize the index that will track where to put different parameters' derivatives
    int i = 0;

    // evaluate derivatives with respect to the cusp condition variable A
    if (Opt_A)
    {
      derivs[i][0] = -ex * R0;   //dU/da
      derivs[i][1] = ex;         //d(dU/da)/dr
      derivs[i][2] = -exoverR0;  //d^2 (dU/da)/dr^2
      ++i;
    }

    // evaluate derivatives with respect to the soft cutoff radius
    if (Opt_R0)
    {

      // sum up terms related to the sigmoidal exansion
      real_type n = 2.0;
      real_type sn = s * s; // s^{n  }
      real_type snm1 = s;   // s^{n-1}
      real_type snm2 = 1.0; // s^{n-2}
      real_type sig_sum_0 = 0.0; // contributes to dUdR0
      real_type sig_sum_1 = 0.0; // contributes to d(dU/dR0)/dr
      real_type sig_sum_2 = 0.0; // contributes to d^2(dU/dR0)/dr^2
      for (int j = 0; j < B.size(); j++) {
        const real_type d = 1.0 + sn;
        const real_type d2 = d * d;
        const real_type soverd = s / d;
        const real_type noverd2 = n / d2;
        sig_sum_0 += B[j] * sn * ( noverd2 - soverd );
        sig_sum_1 += B[j] * snm1 * ( d * noverd2 * noverd2 * ( 1.0 - sn ) - 2.0 * noverd2 * s + soverd * ( s - 1.0 ) );
        sig_sum_2 += B[j] * snm2 * (   s * ( 3.0 * s - 1.0 ) * noverd2 - s * ( s - 2.0 ) * soverd
                                     + ( n + n * sn * ( sn - 4.0 ) + ( 3.0 * s + 1.0 ) * ( sn * sn - 1.0 ) ) * noverd2 * noverd2 );
        snm2 = snm1; // update s^{n-2}
        snm1 = sn;   // update s^{n-1}
        sn *= s;     // update s^{n  }
        n += 1.0;    // update n
      }

      // compute terms related to the cusp-inducing term
      const real_type cusp_part_0 = -A * ex        * ( s + 1.0 ); // contributes to dUdR0
      const real_type cusp_part_1 =  A * exoverR0  * s;           // contributes to d(dU/dR0)/dr
      const real_type cusp_part_2 = -A * exoverR02 * ( s - 1.0 ); // contributes to d^2(dU/dR0)/dr^2

      // compute the desired derivatives
      derivs[i][0] = cusp_part_0 + exoverR0  * sig_sum_0; // dUdR0
      derivs[i][1] = cusp_part_1 + exoverR02 * sig_sum_1; // d(dU/dR0)/dr
      derivs[i][2] = cusp_part_2 + exoverR03 * sig_sum_2; // d^2(dU/dR0)/dr^2

      // increment the index tracking where to put derivatives
      ++i;

    }

    // evaluate derivatives with respect to the sigmoidal expansion coefficients
    if (Opt_B)
    {

      // loop over the terms in the expansion over sigmoidal functions
      real_type n = 2.0;
      real_type sn = s * s; // s^{n  }
      real_type snm1 = s;   // s^{n-1}
      real_type snm2 = 1.0; // s^{n-2}
      for (int j = 0; j < B.size(); j++) {

        // compute some useful values
        const real_type d = 1.0 + sn;
        const real_type d2 = d * d;
        const real_type soverd = s / d;
        const real_type noverd2 = n / d2;

        // dU/dB_j
        derivs[i][0] = -ex * sn / d;

        // d(dU/dB_j)/dr
        derivs[i][1] = exoverR0 * ( sn * sn + sn - n * snm1 ) / d2;

        // d^2(dU/dB_j)/dr^2
        derivs[i][2] = exoverR02 * snm2 * ( d * noverd2 * noverd2 * ( sn - 1.0 ) + noverd2 * ( 1.0 + 2.0 * s ) - soverd * s );

        snm2 = snm1; // update s^{n-2}
        snm1 = sn;   // update s^{n-1}
        sn *= s;     // update s^{n  }
        n += 1.0;    // update n

        // increment the index tracking where to put derivatives
        ++i;

      }

    }

    return true;
  }

  /// compute derivatives of U(r) with respect to the variational parameters
  inline bool evaluateDerivatives(real_type r, std::vector<real_type>& derivs) override
  {

    // get the ratio of the distance and the soft cutoff distance
    const real_type s = r / R0;

    // get the exponential factor
    const real_type ex = std::exp(-s);

    // get some useful ratios
    const real_type exoverR0 = ex / R0;

    // initialize the index that will track where to put different parameters' derivatives
    int i = 0;

    // evaluate derivatives with respect to the cusp condition variable A
    if (Opt_A)
    {
      derivs[i] = -ex * R0;   //dU/da
      ++i;
    }

    // evaluate derivatives with respect to the soft cutoff radius
    if (Opt_R0)
    {

      // sum up terms related to the sigmoidal exansion
      real_type n = 2.0;
      real_type sn = s * s; // s^n
      real_type sig_sum = 0.0;
      for (int j = 0; j < B.size(); j++) {
        const real_type d = 1.0 + sn;
        const real_type d2 = d * d;
        const real_type soverd = s / d;
        const real_type noverd2 = n / d2;
        sig_sum += B[j] * sn * ( noverd2 - soverd );
        sn *= s;  // update s^n
        n += 1.0; // update n
      }

      // compute term related to the cusp-inducing term
      const real_type cusp_part = -A * ex * ( s + 1.0 );

      // compute dU/dR0
      derivs[i] = cusp_part + exoverR0 * sig_sum;

      // increment the index tracking where to put derivatives
      ++i;

    }

    // evaluate derivatives with respect to the sigmoidal expansion coefficients
    if (Opt_B)
    {
      real_type sn = s * s; // s^n
      for (int j = 0; j < B.size(); j++) {
        derivs[i] = -ex * sn / ( 1.0 + sn ); // dU/dB_j
        sn *= s;  // update s^n
        ++i;      // increment the index tracking where to put derivatives
      }
    }

    return true;
  }

  /// set up a variational parameter using the supplied xml node
  template <class U>
  void set_variable_from_xml(ReportEngine & PRE, xmlNodePtr cur, U & variable_to_set, std::string & id_to_set, bool & opt_to_set)
  {

    // get id, name, and optimize info
    std::string       id("string_not_set");
    std::string     name("string_not_set");
    std::string optimize("string_not_set");
    OhmmsAttributeSet rAttrib;
    rAttrib.add(      id,       "id");
    rAttrib.add(    name,     "name");
    rAttrib.add(optimize, "optimize");
    rAttrib.put(cur);

    // read in the variable
    putContent(variable_to_set, cur);

    // set the id if we have it
    if ( id != "string_not_set" )
      id_to_set = id;

    // if we are to optimize the variable, add it to the optimizable variable set
    optimize = lowerCase(optimize);
    if ( optimize == "yes" ) {
      if ( id == "string_not_set" )
        PRE.error("\"id\" must be set if we are going to optimize variable " + name, true);
      opt_to_set = true;
      myVars.insert(id, variable_to_set, opt_to_set, optimize::OTHER_P);
    } else if ( optimize == "no" ) {
      opt_to_set = false;
    } else if ( optimize != "string_not_set" ) {
      PRE.error("Unrecognized value for \"optimize\". Should be either yes or no", true);
    }

  }

  ///read in information about the functor from an xml node
  bool put(xmlNodePtr cur) override
  {

    // set up an object for info / warning / error reporting
    ReportEngine PRE("ShortRangeCuspFunctor", "put(xmlNodePtr)");

    // create some variables to read information in to
    real_type radius = -1.0;    // hard cutoff radius
    real_type cusp_in = -1.0;   // cusp condition

    // read from the xml node
    OhmmsAttributeSet rAttrib;
    rAttrib.add(radius, "rcut");
    rAttrib.add(radius, "cutoff");
    rAttrib.add(cusp_in, "cusp");
    rAttrib.put(cur);

    // set the hard cutoff radius if we have it
    if (radius > 0.0)
      cutoff_radius = radius;

    // set the cusp if we have it
    if (cusp_in > 0.0)
      setCusp(cusp_in);

    // loop over this node's children to read in variational parameter information
    xmlNodePtr childPtr = cur->xmlChildrenNode;
    while (childPtr != NULL)
    {

      // skip blank nodes
      if ( xmlIsBlankNode(childPtr) ) {
        childPtr = childPtr->next;
        continue;
      }

      // get the name of the child node
      std::string cname(lowerCase(castXMLCharToChar(childPtr->name)));
      // read in a variable
      if (cname == "var")
      {

        // read in the name of the variable
        std::string v_name("string_not_set");
        OhmmsAttributeSet att;
        att.add(v_name, "name");
        att.put(childPtr);

        // read in the variable's info
        if (v_name == "A")
          set_variable_from_xml(PRE, childPtr, A, ID_A, Opt_A);
        else if (v_name == "R0")
          set_variable_from_xml(PRE, childPtr, R0, ID_R0, Opt_R0);
        else if (v_name == "string_not_set")
          PRE.error("variable name not set", true);
        else
          PRE.error("unrecognized variable name: " + v_name, true);

      }

      // read in the B coefficients
      else if (cname == "coefficients")
      {

        // read in the id, whether to optimize, and the node type
        std::string       id("string_not_set");
        std::string     type("string_not_set");
        std::string optimize("string_not_set");
        OhmmsAttributeSet att;
        att.add(      id,       "id");
        att.add(    type,     "type");
        att.add(optimize, "optimize");
        att.put(childPtr);

        // sanity check that this node has been specified as an array type
        if (type != "Array")
          PRE.error("ShortRangeCuspFunctor expected coefficients parameter array type to be \"Array\"", true);

        // read in the vector of coefficients
        putContent(B, childPtr);

        // set the id if we have it
        if ( id != "string_not_set" )
          ID_B = id;

        // if the coefficients are to be optimized, add them to the variable set
        optimize = lowerCase(optimize);
        if ( optimize == "yes" ) {
          if ( id == "string_not_set" )
            PRE.error("\"id\" must be set if we are going to optimize B coefficients", true);
          Opt_B = true;
          for (int i = 0; i < B.size(); i++) {
            std::stringstream sstr;
            sstr << ID_B << "_" << i;
            myVars.insert(sstr.str(), B.at(i), Opt_B, optimize::OTHER_P);
          }
        } else if ( optimize == "no" ) {
          Opt_B = false;
        } else if ( optimize != "string_not_set" ) {
          PRE.error("Unrecognized value for \"optimize\". Should be either yes or no", true);
        }

      }

      // error if cname is not recognized
      else {
        PRE.error("\"" + cname + "\" is not a recognized value for \"cname\". Allowed values are \"var\" and \"coefficients\"", true);
      }

      // go to the next node
      childPtr = childPtr->next;
    }

    // summarize what was read in
    app_summary() << "    ---------------------------------------------------------" << std::endl;
    app_summary() << "                       cusp variable A: " << A << std::endl;
    app_summary() << "                        soft cutoff R0: " << R0 << std::endl;
    app_summary() << "    sigmoidal expansion coefficients B:";
    for (int i = 0; i < B.size(); i++)
      app_summary() << " " << B.at(i);
    app_summary() << std::endl;
    app_summary() << "                    hard cutoff radius: " << cutoff_radius << std::endl;
    app_summary() << "                                 Opt_A: " << ( Opt_A  ? "yes" : "no" ) << std::endl;
    app_summary() << "                                Opt_R0: " << ( Opt_R0 ? "yes" : "no" ) << std::endl;
    app_summary() << "                                 Opt_B: " << ( Opt_B  ? "yes" : "no" ) << std::endl;

    // reset the functor now that we've read its info in (currently does nothing)
    reset();

    return true;
  }

  void checkInVariables(opt_variables_type& active) override
  {
    active.insertFrom(myVars);
    //myVars.print(std::cout);
  }

  void checkOutVariables(const opt_variables_type& active) override
  {
    myVars.getIndex(active);
    //myVars.print(std::cout);
  }

  void resetParameters(const opt_variables_type& active) override
  {
    if (myVars.size())
    {
      int ia = myVars.where(0);
      if (ia > -1)
      {
        int i = 0;
        if (Opt_A)
          A = std::real( myVars[i++] = active[ia++] );
        if (Opt_R0)
          R0 = std::real( myVars[i++] = active[ia++] );
        for (int j = 0; Opt_B && j < B.size(); j++)
          B.at(j) = std::real( myVars[i++] = active[ia++] );
      }
      reset();
    }
  }

};

} // namespace qmcplusplus
#endif

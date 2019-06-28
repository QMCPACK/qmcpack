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
#include "Numerics/OptimizableFunctorBase.h"
#include "OhmmsData/AttributeSet.h"
#include <cmath>
#include <stdexcept>
#include "OhmmsPETE/TinyVector.h"

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
  ShortRangeCuspFunctor() : Opt_A(true), Opt_R0(true), Opt_B(true), A(1.0), R0(0.1), ID_A("A"), ID_R0("R0"), ID_B("B") { reset(); }

  ///sets the cusp condition and disables optimization of the cusp-determining parameter
  void setCusp(real_type cusp)
  {
    //throw std::runtime_error("ShortRangeCuspFunctor::setCusp was called");
    A     = cusp;
    Opt_A = false;
    reset();
  }

  ///clone the functor
  OptimizableFunctorBase* makeClone() const { return new ShortRangeCuspFunctor(*this); }

  ///Implement the reset function, which was pure virtual in OptimizableFunctorBase, even though we don't need it
  void reset()
  {
    //cutoff_radius = 1.0e4; //some big range
  }

  ///compute U(r) at a particular value of r
  inline real_type evaluate(real_type r) const
  {

    // get the ratio of the distance and the soft cutoff distance
    const real_type s = r / R0;

    // sum up the sigmoidal function expansion
    real_type sig_sum = 0.0;
    real_type n = 2.0;
    real_type sn = s * s; // s^n
    for (int i = 0; i < B.size(); i++) {
      sig_sum += B[i] * sn / ( 1.0 + sn );
      sn *= s;  // update s^n
      n += 1.0; // update n
    }

    // return U(r)
    return -1.0 * std::exp(-s) * ( A * R0 + sig_sum );

  }

  ///compute U(r), dU/dr, and d^2U/dr^2 at a particular value of r
  inline real_type evaluate(real_type r, real_type& dudr, real_type& d2udr2) const
  {

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
    throw std::runtime_error("evaluate d3udr3 not implemented for ShortRangeCuspFunctor");
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
  inline real_type f(real_type r)
  {
    if (r >= cutoff_radius)
      return 0.0;
    return evaluate(r);
  }

  ///compute dU/dr at a particular distance or return zero if beyond the cutoff
  inline real_type df(real_type r)
  {
    if (r >= cutoff_radius)
      return 0.0;
    real_type du, d2u;
    evaluate(r, du, d2u);
    return du;
  }

  /// compute derivatives of U(r), dU/dr, and d^2U/dr^2 with respect to the variational parameters
  inline bool evaluateDerivatives(real_type r, std::vector<TinyVector<real_type, 3>>& derivs)
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
  inline bool evaluateDerivatives(real_type r, std::vector<real_type>& derivs)
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
      real_type n = 2.0;
      real_type sn = s * s; // s^n
      for (int j = 0; j < B.size(); j++) {
        derivs[i] = -ex * sn / ( 1.0 + sn ); // dU/dB_j
        sn *= s;  // update s^n
        n += 1.0; // update n
        ++i;      // increment the index tracking where to put derivatives
      }
    }

    return true;
  }

  ///read in information about the functor from an xml node
  bool put(xmlNodePtr cur)
  {

    // set up an object for info / warning / error reporting
    ReportEngine PRE("ShortRangeCuspFunctor", "put(xmlNodePtr)");

    // create some variables to read information in to
    int nB = -1;                // number of terms in the expansion over sigmoidal functions
    real_type radius = -1.0;    // hard cutoff radius
    std::string name("0");      // a name used in naming the variational parameters
    std::string optA("no");     // whether to optimize the cusp condition (off by default)
    std::string optR0("yes");   // whether to optimize the soft cutoff radius (on by default)

    // read from the xml node
    OhmmsAttributeSet rAttrib;
    rAttrib.add(nB, "size");
    rAttrib.add(radius, "rcut");
    rAttrib.add(radius, "cutoff");
    rAttrib.add(A, "cusp");
    rAttrib.add(R0, "r0");
    rAttrib.add(optA, "optA");
    rAttrib.add(optR0, "optR0");
    rAttrib.add(name, "elementType");
    rAttrib.put(cur);

    // set the hard cutoff radius
    if (radius > 0.0)
      cutoff_radius = radius;

    // sanity check for the length of the sigmoidal expansion
    if (nB < 0)
      PRE.error("Number of sigmoidal expansion terms must be non-negative.  Did you forget to specify \"size\"?", true);

    // sanity check for whether we found a name
    if (name == "0")
      PRE.error("ShortRangeCuspFunctor did not find an acceptable name.", true);

    // set the strings that act as IDs for the different variational parameters
    ID_A  = name + "_SRC_A";
    ID_B  = name + "_SRC_B";
    ID_R0 = name + "_SRC_R0";

    // set which of A and R0 are optimizable and add whichever are to the variable set
    tolower(optA);
    tolower(optR0);
    Opt_A  = ( optA  == "yes" );
    Opt_R0 = ( optR0 == "yes" );
    if ( Opt_A )
      myVars.insert(ID_A, A, Opt_A, optimize::OTHER_P);
    if ( Opt_R0 )
      myVars.insert(ID_R0, R0, Opt_R0, optimize::OTHER_P);

    // loop over this node's children to find the the B coefficents for the sigmoidal expansion
    xmlNodePtr xmlCoefs = cur->xmlChildrenNode;
    while (xmlCoefs != NULL)
    {

      // read from the node named "coefficients"
      std::string cname((const char*)xmlCoefs->name);
      if (cname == "coefficients")
      {

        // strings to read in to
        std::string type("0");
        std::string optB("yes");

        // read attributes of this node
        OhmmsAttributeSet cAttrib;
        cAttrib.add(type, "type");
        cAttrib.add(optB, "optimize");
        cAttrib.put(xmlCoefs);

        // sanity check that this node has been specified as an array type
        if (type != "Array")
          PRE.error("ShortRangeCuspFunctor expected B paramter array type to be \"Array\"", true);

        // read in the vector of coefficients
        std::vector<real_type> params;
        putContent(params, xmlCoefs);

        // sanity check that the vector was the expected length
        if (params.size() != nB)
          PRE.error("ShortRangeCuspFunctor encountered a B parameter array that is the wrong length.", true);

        // store the coefficients, and, if they are to be optimized, add them to the variable set
        B = params;
        tolower(optB);
        Opt_B = ( optB == "yes" );
        if ( Opt_B ) {
          for (int i = 0; i < B.size(); i++)
          {
            std::stringstream sstr;
            sstr << ID_B << "_" << i;
            myVars.insert(sstr.str(), B.at(i), Opt_B, optimize::OTHER_P);
          }
        }

        // stop looping if we've found the coefficients
        break;
      }

      // go to the next node
      xmlCoefs = xmlCoefs->next;
    }

    // summarize what was read in
    app_summary() << "                   cusp variable A: " << A << std::endl;
    app_summary() << "                    soft cutoff R0: " << R0 << std::endl;
    app_summary() << "sigmoidal expansion coefficients B:";
    for (int i = 0; i < B.size(); i++)
      app_summary() << " " << B.at(i);
    app_summary() << std::endl;
    app_summary() << "                hard cutoff radius: " << cutoff_radius << std::endl;
    app_summary() << "                             Opt_A: " << ( Opt_A  ? "yes" : "no" ) << std::endl;
    app_summary() << "                            Opt_R0: " << ( Opt_R0 ? "yes" : "no" ) << std::endl;
    app_summary() << "                             Opt_B: " << ( Opt_B  ? "yes" : "no" ) << std::endl;

    // reset the functor now that we've read its info in (currently does nothing)
    reset();

    return true;
  }

  void checkInVariables(opt_variables_type& active)
  {
    active.insertFrom(myVars);
    //myVars.print(std::cout);
  }

  void checkOutVariables(const opt_variables_type& active)
  {
    myVars.getIndex(active);
    //myVars.print(std::cout);
  }

  void resetParameters(const opt_variables_type& active)
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

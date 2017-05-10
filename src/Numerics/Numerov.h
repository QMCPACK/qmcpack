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
    
    



#ifndef NUMEROVCLASS_H
#define NUMEROVCLASS_H

/**@defgroup NumerovTransform Numerov Transformation
 *@brief Transform functors which provide interfaces to Numerov<ST,FT>.
 *
 *The prototype: "template<class SourceType> ATransformFunctor",
 *where SourceType encapsulates the form of an external potential.
 *
 *For examples, the RegularLinearTransform is implemented to solve
 *Radial Schrodinger equation
 *\f[ \frac{d^2R_{nl}}{dr^2}+k^2(r) R_{nl} = 0, \f] where
 *\f[k^2(r)=2[\varepsilon-\frac{L(L+1)}{2r^2}-V(r)], \f] and
 *\f$ V(r) \f$ is a radial potential on a linear grid.
 *
 */
/** Generic class to solve eigen problems using the Numerov
 *Integration Algorithm.
 *
 *Solve second order differential equations of the form
 \f[\frac{d^2y}{dx^2} + k^2(x)y = S(x)\f]
 using the Numerov algorithm
 \f[
 (1+\frac{h^2}{12}k^2_{n+2})y_{n+2} -
 2(1-\frac{5h^2}{12}k^2_{n+1})y_{n+1} + (1+\frac{h^2}{12}k^2_n)y_n =
 \frac{h^2}{12}(S_{n+2}+10S_{n+1}+S_n) + \mathcal{O}(h^6), \f] by
 solving for \f$y_{n+2}\f$ and recursively integrating forward in x.
 The boundary conditons being \f$y_0 = y(x_0)\f$ and \f$y_1 = y(x_1).\f$
 *
 * Two template parameters are required:
 - ST: source functor type,  e.g., the transform functors in
 \ref NumerovTransform.
 - FT: grid functor type for the solution \f$y(x)\f$ and \f$V(x)\f$
 for \f$k^2(\epsilon,V(x))\f$.
 *
 *See \ref numerov_sec "Numerov algorihm".
 *
 *Specifically, the source functor should implement
 - value_type ST::setCusp(int i, value_type& z0, value_type& z1)
 - int ST::first() the first valid index for the grid
 - value_type ST::k2(i)
 - value_type ST::convert(value_type Z(x), x)
 *
 */
template<class ST, class FT>
class Numerov
{

public:

  enum {lower_energy=-1, nothing, raise_energy};

  typedef typename FT::value_type value_type;
  typedef typename FT::data_type data_type;

  ///constructor with source V(x) and target solution f(x) and TransForm y
  Numerov(ST& y, FT& f): Y(y), Target(f)
  {
    Z.resize(Target.size());
    //save the derivate of the first grid
    Y.setCusp(Y.first(), Target0, Target1);
  }

  int evaluate(value_type e);

  /**
   *\param lowerin lower bound for the eigen energy
   *\param upperin upper bound for the eigen energy
   *\param etol the energy tolerance
   *\return the eigen energy
   *\brief Solve for the eigen energy and eigen function self
   consistently.
  */
  inline value_type solve(value_type lowerin, value_type upperin,
                          value_type etol)
  {
    value_type lower = lowerin;
    value_type upper = upperin;
    value_type trialenergy=0.5*(lower+upper);
    value_type oldenergy = upper;
    int modifier = nothing;
    do
    {
      modifier = evaluate(trialenergy);
      if(modifier == raise_energy)
        lower = trialenergy;
      if(modifier == lower_energy)
        upper = trialenergy;
      oldenergy = trialenergy;
      trialenergy=0.5*(lower+upper);
    }
    while((std::abs(trialenergy-oldenergy))>etol);
    return trialenergy;
  }

  /**
   *\param e reference eigenvalue
   *\brief Reset the reference eigen energy for \f$k^2(\epsilon, V(x))\f$
   *and the classical TurningPoint \f$i\f$, where
   *\f$k^2(\epsilon,V(x_{i})) < 0\f$ and \f$k^2(\epsilon,V(x_{i-1})) > 0\f$
   */
  inline void reset(value_type e)
  {
    Target.m_Y = 0.0;
    TurningPoint = 0;
    Y.reset(e);
    int i=Target.size()-2;
    while(i > 1 && Y.k2(i)<0)
    {
      i--;
    }
    TurningPoint = i;
  }

private:

  ///engine that takes care of the transformation of variables
  ST& Y;

  ///real function that holds solution
  FT& Target;

  ///solution by Numerov  Target <-> TransForm::convert(Z)
  data_type Z;

  ///cusp condition
  value_type Target0, Target1;

  ///the maximum size of physical data
  int TurningPoint;

  ///index reserved to truncate the data
  int Current;
};


/**
 *\param e reference eigen energy for \f$k^2(\epsilon, V(x))\f$
 *\return integer flag for the change of the reference eigen energy
 */
template<class TransForm, class FT>
inline
int Numerov<TransForm,FT>::evaluate(value_type e)
{
  reset(e);
  const value_type MAX_VALUE = 1000;
  value_type dh2 = Target.dh()*Target.dh();
  value_type tentwelfth = 0.8333333333333333333333333333333333333333333333*dh2;
  value_type onetwelfth = 0.0833333333333333333333333333333333333333333333*dh2;
  //get the first index
  int first = Y.first();
  int second = first+1;
  Z[first] = Target0;
  Z[second] = Target1;
  value_type r0 = Target.r(first);
  value_type r1 = Target.r(second);
  value_type k2_m2 = Y.k2(first);
  value_type k2_m  = Y.k2(second);
  value_type y_m2  = Z[first];
  value_type y_m   = Z[second];
  //set the wave function
  Target(first) = Y.convert(y_m2,r0);
  Target(second) = Y.convert(y_m,r1);
  int num_nodes = 0;
  int nodes_to_find = Y.nodes();
  //not meaningful turning point, raise the energy
  if(TurningPoint == 1)
    return raise_energy;
  int e_mode =nothing;
  int i=first+2;
  while(i<=TurningPoint+2)
  {
    //get the radius of the current index
    r0 = Target.r(i);
    //get k^2
    value_type k2 = Y.k2(i);
    value_type y =
      ((2.0-tentwelfth*k2_m)*y_m -
       (1.0+onetwelfth*k2_m2)*y_m2)/(1.0+onetwelfth*k2);
    Z[i] = y;
    //avoid exponential solution
    if(std::abs(y)>MAX_VALUE)
    {
      value_type yinv = 1.0/y;
      for(int j=0; j<=i; j++)
        Z[j]*=yinv;
      y_m2*=yinv;
      y_m*=yinv;
      y = 1.0;
    }
    //check the node
    if(y_m*y < 0.0)
    {
      num_nodes++;
      if( num_nodes > nodes_to_find)
      {
        return lower_energy;
      }
    }
    k2_m2 = k2_m;
    k2_m = k2;
    y_m2 = y_m;
    y_m = y;
    Target(i) = Y.convert(Z[i],r0);
    i++;
  }
  if(e_mode == nothing)
  {
    i = TurningPoint+3;
    while(i<=Target.size()-2)
    {
      r0 = Target.r(i);
      //get k^2
      value_type k2 = Y.k2(i);
      value_type y
      = ((2.0-tentwelfth*k2_m)*y_m - (1.0+onetwelfth*k2_m2)*y_m2)
        /(1.0+onetwelfth*k2);
      Z[i] = y;
      //growing beyond the turning point
      if(y/y_m > 1)
      {
        return raise_energy;
      }
      if(std::abs(y)>MAX_VALUE)
      {
        value_type yinv = 1.0/y;
        for(int j=0; j<=i; j++)
          Z[j]*=yinv;
        y_m2*=yinv;
        y_m*=yinv;
        y = 1.0;
      }
      if(y_m*y < 0.0)
      {
        num_nodes++;
        if(num_nodes > nodes_to_find)
          return lower_energy;
      }
      //assign numerov value to the real function
      k2_m2 = k2_m;
      k2_m = k2;
      y_m2 = y_m;
      y_m = y;
      Target(i) = Y.convert(Z[i],r0);
      i++;
    }
  }
  if( num_nodes < nodes_to_find )
    return raise_energy;
  else
    return nothing;
}

#endif

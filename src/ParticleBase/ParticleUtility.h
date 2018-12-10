//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#ifndef OHMMS_PARTICLEUTILITY_H
#define OHMMS_PARTICLEUTILITY_H

namespace qmcplusplus
{
//////////////////////////////////////////////////////////
// functors to generate vectors with each element [-0.5,0.5)
//////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// specialized of RandomVector in Utilities/RandomGenerator.h for TinyVector<T,D>
//////////////////////////////////////////////////////////

// template<class T, unsigned D, class RNG>
// struct
// RandomVector<TinyVector<T,D>, RNG> {

//   typedef TinyVector<T,D> Return_t;

//   RandomVector(RNG& gen): RG(gen){ }

//   inline Return_t operator()() {
//     Return_t res(0);
//     for(int i=0; i<D; i++) res[i] = 2*RG()-1;
//     return res;
//   }

//   inline Return_t norm() {
//     Return_t res;
//     for(int i=0; i<D; i++) res[i] = RG();
//     T norm = 1.0/sqrt(dot(res,res));
//     return norm*res;
//   }

//   inline void get(Return_t& a, Return_t& b) {
//     for(int i=0; i<D; i++) RG.bivariate(a[i], b[i]);
//   }

//   RNG& RG;
// };


// // specialized for TinyVector<T,D>
// template<class RNG>
// struct RandomVector<TinyVector<double,3>, RNG> {
//   typedef TinyVector<double,3> Return_t;

//   RandomVector(RNG& gen): RG(gen){ }

//   inline Return_t operator()() {
//     return Return_t(2*RG()-1, 2*RG()-1, 2*RG()-1);
//   }

//   inline Return_t norm() {
//     double x,y,s;
//     s = 2.0e0;
//     while(s > 1.0e0) {
//       x = 2.0e0*RG() - 1.0e0;  // [-0.5,0.5)
//       y = 2.0e0*RG() - 1.0e0;  // [-0.5,0.5)
//      s = x*x + y*y;
//     }
//     double z = 2.0e0*sqrt(1.0e0-s);
//     return Return_t(z*x, z*y, 1.0e0-2.0e0*s);
//   }

//   inline void get(Return_t& a, Return_t& b) {
//     RG.bivariate(a[0], b[0]);
//     RG.bivariate(a[1], b[1]);
//     RG.bivariate(a[2], b[2]);
//   }


//   RNG& RG;
// };

// //////////////////////////////////////////////////////////
// // Templated function to generate random Vector<T>
// //////////////////////////////////////////////////////////
// // forward declaration
// template<class T, unsigned D> class CrystalLattice;
// template<class T> class ParticleAttrib;

// ///////////////////////////////////////////////////////////
// // specialization of RandomSqence with ParticleAttrib<T>
// ///////////////////////////////////////////////////////////
// template<class T, class RNG>
// struct RandomSequence<ParticleAttrib<T>, RNG> {

//   static void apply(ParticleAttrib<T>& v, RandomVector<T,RNG>& rnd) {
//     for(int i=0; i<v.size(); i++) v[i] = rnd();
//   }
// };


// ///////////////////////////////////////////////////////////
// // Functor to generate a sequence of Normalized vectors
// ///////////////////////////////////////////////////////////

// template<class RV, class RNG>
// struct RandomNormSequence { };

// template<class T, class RNG>
// struct RandomNormSequence<ParticleAttrib<T>, RNG> {

//   static void apply(ParticleAttrib<T>& v, RandomVector<T,RNG>& rnd) {
//     for(int i=0; i<v.size(); i++) v[i] = rnd.norm();
//   }
// };

////////////////////////////////////////////////////////////////
// Iterator is exposed. Parallel Implementation requires special care
////////////////////////////////////////////////////////////////
template<class PL, class PV>
void convert(const PL& lat, const PV& pin, PV& pout)
{
  if(pin.InUnit == pout.InUnit)
  {
    pout = pin;
    return;
  }
  if(pin.InUnit)
  {
    for(int i=0; i<pin.size(); i++)
      pout[i] = lat.toCart(pin[i]);
    return;
  }
  else
  {
    for(int i=0; i<pin.size(); i++)
      pout[i] = lat.toUnit(pin[i]);
    return;
  }
}

////////////////////////////////////////////////////////////////
// Iterator is exposed. Parallel Implementation requires special care
////////////////////////////////////////////////////////////////
template<class PL, class PV>
void convert2Cart(const PL& lat, PV& pin)
{
  if(pin.InUnit)
  {
    PV tmp(pin.size());
    tmp = pin;
    pin.InUnit = false;
    for(int i=0; i<pin.size(); i++)
      pin[i] = lat.toCart(pin[i]);
  }
}

template<class PL, class PV>
void convert2Unit(const PL& lat, PV& pin)
{
  if(!pin.InUnit)
  {
    PV tmp(pin.size());
    tmp = pin;
    pin.InUnit = true;
    for(int i=0; i<pin.size(); i++)
      pin[i] = lat.toUnit(pin[i]);
  }
}

////////////////////////////////////////////////////////////////
// Apply BC conditions to put the position type in lattice box [0,1)
////////////////////////////////////////////////////////////////
template<class PL, class PV>
void wrapAroundBox(const PL& lat, const PV& pin, PV& pout)
{
  if(pin.InUnit)
  {
    if(pout.InUnit)
    {
      for(int i=0; i<pin.size(); i++)
      {
        pout[i] = lat.BConds.wrap(pin[i]);            //unit -> unit
      }
    }
    else
    {
      for(int i=0; i<pin.size(); i++)
        pout[i] = lat.toCart(lat.BConds.wrap(pin[i]));//unit -> cart
    }
  }
  else
  {
    if(pout.InUnit)
    {
      for(int i=0; i<pin.size(); i++)
        pout[i] = lat.BConds.wrap(lat.toUnit(pin[i]));//cart -> unit
    }
    else
    {
      for(int i=0; i<pin.size(); i++)
        pout[i] = lat.toCart(lat.BConds.wrap(lat.toUnit(pin[i])));//cart -> cart
    }
  }
}

/////////////////////////////////////////////////////////////////
/*\fn template<class T, unsigned D>
 *    T Dot(const ParticleAttrib<TinyVector<T, D> >& pa,
 *      const ParticleAttrib<TinyVector<T, D> >& pb)
 * \return a dot product of an array
 */
/////////////////////////////////////////////////////////////////
template<typename T, unsigned D>
inline T Dot(const ParticleAttrib<TinyVector<T, D> >& pa,
                  const ParticleAttrib<TinyVector<T, D> >& pb)
{
  T sum = 0;
  for(int i=0; i<pa.size(); i++)
  {
    sum += dot(pa[i],pb[i]);
  }
  return sum;
}

template<typename T>
inline T Sum (const ParticleAttrib<T>& pa)
{
  T sum = 0;
  for(int i=0; i<pa.size(); i++)
  {
    sum += pa[i];
  }
  return sum;
}

template<class T, unsigned D>
void  normalize(ParticleAttrib<TinyVector<T, D> >& pa)
{
  T factor = Dot(pa,pa);
  factor = 1.0/sqrt(factor);
  pa *= factor;
}

}
#endif

